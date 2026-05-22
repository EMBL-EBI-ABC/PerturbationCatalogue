import argparse
import csv
import json
import traceback
import instructor
from pathlib import Path
import os
from datetime import datetime
from threading import Lock
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Type
from pydantic import BaseModel
from controlled_vocab_model import MavedbMetadataExtractionSchema, MavedbMetadataSchema
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]

PROMPT_TEMPLATE_FILE = SCRIPT_DIR / "metadata_extraction_prompt_template.md"
MAVEDB_METADATA_DIR = SCRIPT_DIR / "mavedb_metadata"
MAVEDB_URN_TO_DOIS_FILE = SCRIPT_DIR / "mavedb_urn_to_dois.json"

DEFAULT_LLM_MODEL_NAME = os.getenv("LLM_MODEL_NAME", "google/gemini-2.5-flash")
EXTRACTED_METADATA_OUTPUT_DIR = SCRIPT_DIR / "extracted_metadata"
DEFAULT_CURATED_METADATA_CSV_PATH = EXTRACTED_METADATA_OUTPUT_DIR / "clean_metadata.csv"
PUBLICATION_FULL_TEXT_MD_DIR = SCRIPT_DIR / "pub_full_text_md"
METADATA_EXTRACTION_LOG_FILE = SCRIPT_DIR / 'metadata_extraction.log'
METADATA_EXTRACTION_LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
DEFAULT_DOWNLOAD_MAX_WORKERS = min(32, (os.cpu_count() or 1) * 4)
DEFAULT_BULK_EXCLUDED_PUBLICATION_FILES = {
    "10_1101_2024_04_26_591310.md",
}
PRINT_LOCK = Lock()
LOG_SEPARATOR_WIDTH = 80
JSON_INDENT = 2

_MAVEDB_URN_TO_DOIS_CACHE: dict[str, list[str]] | None = None

EXAMPLE_PUBLICATION_FILE = PUBLICATION_FULL_TEXT_MD_DIR / "10_1016_j_ajhg_2024_02_002.md"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--llm-model",
        default=DEFAULT_LLM_MODEL_NAME,
        help="LLM model ID to use for metadata extraction.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_DOWNLOAD_MAX_WORKERS,
        help="Number of worker threads to use for bulk extraction.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing extracted metadata outputs.",
    )
    parser.add_argument(
        "--create-csv",
        action="store_true",
        help="Create a single CSV from the curated clean metadata JSON outputs.",
    )
    return parser

def _append_log_block(title: str, *lines: str) -> None:
    timestamp = datetime.now().isoformat(timespec='seconds')
    separator = "=" * LOG_SEPARATOR_WIDTH
    with METADATA_EXTRACTION_LOG_FILE.open('a', encoding='utf-8') as log_file:
        log_file.write(f"\n[{timestamp}] {separator}\n")
        log_file.write(f"[{timestamp}] {title}\n")
        for line in lines:
            log_file.write(f"[{timestamp}] {line}\n")
        log_file.write(f"[{timestamp}] {separator}\n")


def _append_log_line(message: str) -> None:
    timestamp = datetime.now().isoformat(timespec='seconds')
    with METADATA_EXTRACTION_LOG_FILE.open('a', encoding='utf-8') as log_file:
        log_file.write(f"[{timestamp}] {message}\n")

def print_status_block(title: str, *lines: str, echo: bool = True) -> None:
    separator = "=" * LOG_SEPARATOR_WIDTH
    with PRINT_LOCK:
        _append_log_block(title, *lines)
        if not echo:
            return
        print(f"\n{separator}")
        print(title)
        for line in lines:
            print(line)
        print(separator)


def _normalize_prompt_text(value: str | None) -> str | None:
    if not value:
        return None
    normalized_value = " ".join(str(value).split())
    return normalized_value or None


def _format_identifier_for_lookup(identifier: str) -> str:
    return str(identifier).replace(".", "_").replace("/", "_")


def _format_urn_for_filename(urn: str) -> str:
    return str(urn).replace(":", "_")


def _load_mavedb_urn_to_dois() -> dict[str, list[str]]:
    global _MAVEDB_URN_TO_DOIS_CACHE
    if _MAVEDB_URN_TO_DOIS_CACHE is None:
        _MAVEDB_URN_TO_DOIS_CACHE = json.loads(
            MAVEDB_URN_TO_DOIS_FILE.read_text(encoding="utf-8")
        )
    return _MAVEDB_URN_TO_DOIS_CACHE


def _find_matching_mavedb_entry_paths(
    publication_full_text_path: str | Path,
) -> list[tuple[str, Path]]:
    publication_full_text_path = Path(publication_full_text_path).resolve()
    publication_identifier = publication_full_text_path.stem
    urn_to_dois = _load_mavedb_urn_to_dois()
    matching_entries: list[tuple[str, Path]] = []
    for urn, identifiers in sorted(urn_to_dois.items()):
        if any(_format_identifier_for_lookup(identifier) == publication_identifier for identifier in identifiers):
            entry_path = MAVEDB_METADATA_DIR / f"{urn.replace(':', '_')}.json"
            if entry_path.is_file():
                matching_entries.append((urn, entry_path))
    return matching_entries


def _extract_curated_mavedb_prompt_metadata(entry_payload: dict) -> dict[str, object]:
    experiment_payload = entry_payload.get("experiment") or {}
    target_genes = sorted(
        {
            gene.get("name")
            for gene in entry_payload.get("targetGenes", [])
            if gene.get("name")
        }
    )
    publication_titles = sorted(
        {
            publication.get("title")
            for publication in entry_payload.get("primaryPublicationIdentifiers", [])
            if publication.get("title")
        }
    )

    curated_metadata: dict[str, object] = {
        "score_set_title": _normalize_prompt_text(entry_payload.get("title")),
        "score_set_short_description": _normalize_prompt_text(entry_payload.get("shortDescription")),
        "score_set_abstract": _normalize_prompt_text(entry_payload.get("abstractText")),
        "score_set_method": _normalize_prompt_text(entry_payload.get("methodText")),
        "experiment_title": _normalize_prompt_text(experiment_payload.get("title")),
        "experiment_short_description": _normalize_prompt_text(experiment_payload.get("shortDescription")),
        "experiment_abstract": _normalize_prompt_text(experiment_payload.get("abstractText")),
        "experiment_method": _normalize_prompt_text(experiment_payload.get("methodText")),
        "target_genes": target_genes or None,
        "score_columns": entry_payload.get("datasetColumns", {}).get("scoreColumns") or None,
        "primary_publication_titles": publication_titles or None,
    }
    return {
        field_name: field_value
        for field_name, field_value in curated_metadata.items()
        if field_value not in (None, [], {})
    }


def _build_mavedb_context_signature(curated_metadata: dict[str, object]) -> str:
    signature_payload = {
        "score_set_short_description": curated_metadata.get("score_set_short_description"),
        "score_set_method": curated_metadata.get("score_set_method"),
        "experiment_short_description": curated_metadata.get("experiment_short_description"),
        "experiment_method": curated_metadata.get("experiment_method"),
    }
    if not any(signature_payload.values()):
        signature_payload = {
            "score_set_title": curated_metadata.get("score_set_title"),
            "experiment_title": curated_metadata.get("experiment_title"),
        }
    return json.dumps(signature_payload, sort_keys=True, ensure_ascii=True)


def _merge_prompt_metadata_value(existing_value: object, new_value: object) -> object:
    if existing_value == new_value or new_value in (None, [], {}):
        return existing_value
    if existing_value in (None, [], {}):
        return new_value

    if isinstance(existing_value, list):
        merged_values = list(existing_value)
    else:
        merged_values = [existing_value]

    if isinstance(new_value, list):
        for value in new_value:
            if value not in merged_values:
                merged_values.append(value)
    elif new_value not in merged_values:
        merged_values.append(new_value)
    return merged_values


def _build_mavedb_prompt_contexts(
    publication_full_text_path: str | Path,
) -> list[dict[str, object]]:
    matching_entries = _find_matching_mavedb_entry_paths(publication_full_text_path)
    if not matching_entries:
        return []

    grouped_contexts: dict[str, dict[str, object]] = {}
    for urn, entry_path in matching_entries:
        entry_payload = json.loads(entry_path.read_text(encoding="utf-8"))
        curated_metadata = _extract_curated_mavedb_prompt_metadata(entry_payload)
        context_signature = _build_mavedb_context_signature(curated_metadata)
        if context_signature not in grouped_contexts:
            grouped_contexts[context_signature] = {
                "source_urns": [urn],
                "source_files": [entry_path.name],
                "metadata": curated_metadata,
            }
            continue

        grouped_contexts[context_signature]["source_urns"].append(urn)
        grouped_contexts[context_signature]["source_files"].append(entry_path.name)
        merged_metadata = grouped_contexts[context_signature]["metadata"]
        for field_name, field_value in curated_metadata.items():
            merged_metadata[field_name] = _merge_prompt_metadata_value(
                merged_metadata.get(field_name),
                field_value,
            )

    return sorted(
        grouped_contexts.values(),
        key=lambda context: tuple(context["source_urns"]),
    )


def _format_supplementary_mavedb_metadata(
    prompt_context: dict[str, object] | None,
) -> str:
    if not prompt_context:
        return "No supplementary MaveDB metadata was available for this publication."

    return json.dumps(
        {
            "source_urns": prompt_context["source_urns"],
            "source_files": prompt_context["source_files"],
            "curated_mavedb_metadata": prompt_context["metadata"],
        },
        indent=JSON_INDENT,
        ensure_ascii=True,
    )


def _build_context_output_suffix(
    prompt_context: dict[str, object] | None,
    context_index: int,
    total_contexts: int,
) -> str:
    base_suffix = "" if total_contexts == 1 else f"__ctx_{context_index:02d}"
    if not prompt_context or not prompt_context.get("source_urns"):
        return base_suffix

    urn_suffix = "__" + "__".join(
        _format_urn_for_filename(urn) for urn in prompt_context["source_urns"]
    )
    return f"{base_suffix}{urn_suffix}" if base_suffix else urn_suffix


def _build_output_metadata(prompt_context: dict[str, object] | None) -> dict[str, object]:
    if not prompt_context:
        return {}

    return {
        "__source_urns": list(prompt_context.get("source_urns", [])),
        "__source_files": list(prompt_context.get("source_files", [])),
    }


def _filter_bulk_publication_paths(
    publication_full_text_paths: list[str | Path],
) -> tuple[list[Path], list[Path]]:
    publication_paths = [Path(path).resolve() for path in publication_full_text_paths]
    excluded_paths: list[Path] = []
    included_paths: list[Path] = []
    for publication_path in publication_paths:
        if publication_path.name in DEFAULT_BULK_EXCLUDED_PUBLICATION_FILES:
            excluded_paths.append(publication_path)
            continue
        included_paths.append(publication_path)
    return included_paths, excluded_paths


def get_metadata_output_paths(
    publication_full_text_path: str | Path,
    output_dir: str | Path = EXTRACTED_METADATA_OUTPUT_DIR,
    suffix: str = "",
) -> tuple[Path, Path]:
    publication_full_text_path = Path(publication_full_text_path).resolve()
    publication_stem = publication_full_text_path.stem
    output_dir = Path(output_dir).resolve()
    with_evidence_path = output_dir / "with_evidence" / f"{publication_stem}{suffix}.json"
    clean_path = output_dir / "clean" / f"{publication_stem}{suffix}.json"
    return with_evidence_path, clean_path


def _format_metadata_csv_cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=True, sort_keys=True)
    return str(value)


def create_csv_from_curated_metadata_json(
    input_dir: str | Path = EXTRACTED_METADATA_OUTPUT_DIR / "clean",
    output_csv_path: str | Path = DEFAULT_CURATED_METADATA_CSV_PATH,
) -> Path:
    input_dir = Path(input_dir).resolve()
    output_csv_path = Path(output_csv_path).resolve()

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Curated metadata directory not found: {input_dir}")

    json_paths = sorted(input_dir.glob("*.json"))
    if not json_paths:
        raise FileNotFoundError(f"No curated metadata JSON files found in: {input_dir}")

    rows: list[dict[str, str]] = []
    fieldnames: list[str] = ["source_json_file"]
    seen_fieldnames = set(fieldnames)

    for json_path in json_paths:
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        row = {"source_json_file": json_path.name}
        for field_name, field_value in payload.items():
            row[field_name] = _format_metadata_csv_cell(field_value)
            if field_name not in seen_fieldnames:
                fieldnames.append(field_name)
                seen_fieldnames.add(field_name)
        rows.append(row)

    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    with output_csv_path.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print_status_block(
        "Curated metadata CSV created",
        f"Input directory: {input_dir}",
        f"JSON files processed: {len(json_paths)}",
        f"CSV output: {output_csv_path}",
    )
    return output_csv_path
        

def save_metadata_outputs(
    publication_full_text_path: str | Path,
    extraction_result_with_evidence: BaseModel,
    clean_result: MavedbMetadataSchema,
    output_dir: str | Path = EXTRACTED_METADATA_OUTPUT_DIR,
    suffix: str = "",
    prompt_context: dict[str, object] | None = None,
) -> tuple[Path, Path]:
    publication_full_text_path = Path(publication_full_text_path).resolve()
    output_dir = Path(output_dir).resolve()
    with_evidence_path, clean_path = get_metadata_output_paths(
        publication_full_text_path=publication_full_text_path,
        output_dir=output_dir,
        suffix=suffix,
    )
    with_evidence_dir = with_evidence_path.parent
    clean_dir = clean_path.parent
    with_evidence_dir.mkdir(parents=True, exist_ok=True)
    clean_dir.mkdir(parents=True, exist_ok=True)

    output_metadata = _build_output_metadata(prompt_context)
    with_evidence_payload = extraction_result_with_evidence.model_dump()
    clean_payload = clean_result.model_dump()
    with_evidence_payload.update(output_metadata)
    clean_payload.update(output_metadata)

    with_evidence_path.write_text(json.dumps(with_evidence_payload, indent=JSON_INDENT))
    clean_path.write_text(json.dumps(clean_payload, indent=JSON_INDENT))
    print_status_block(
        "Metadata outputs saved",
        f"Publication text: {publication_full_text_path}",
        f"With evidence: {with_evidence_path}",
        f"Clean output: {clean_path}",
    )
    return with_evidence_path, clean_path


def to_final_metadata_schema(
    extraction_result_with_evidence: BaseModel,
) -> MavedbMetadataSchema:
    final_payload = {
        field_name: field_value
        for field_name, field_value in extraction_result_with_evidence.model_dump().items()
        if not field_name.endswith("_evidence")
    }
    return MavedbMetadataSchema.model_validate(final_payload)

def build_metadata_extraction_prompt(
    publication_full_text_path: str | Path,
    prompt_context: dict[str, object] | None = None,
) -> str:
    publication_full_text_path = Path(publication_full_text_path).resolve()
    publication_full_text = publication_full_text_path.read_text(encoding='utf-8')
    _append_log_line(
        f"Loaded publication text from {publication_full_text_path}; characters: {len(publication_full_text)}"
    )
    prompt_template = PROMPT_TEMPLATE_FILE.read_text(encoding='utf-8')
    supplementary_mavedb_metadata = _format_supplementary_mavedb_metadata(prompt_context)
    return prompt_template.format(
        supplementary_mavedb_metadata=supplementary_mavedb_metadata,
        publication_full_text=publication_full_text,
    )


def _extract_metadata_for_prompt_context(
    publication_full_text_path: Path,
    output_dir: Path,
    overwrite: bool,
    model_name: str,
    extraction_schema: Type[BaseModel],
    prompt_context: dict[str, object] | None,
    output_suffix: str,
) -> MavedbMetadataSchema:
    with_evidence_output_path, clean_output_path = get_metadata_output_paths(
        publication_full_text_path=publication_full_text_path,
        output_dir=output_dir,
        suffix=output_suffix,
    )

    if not overwrite and clean_output_path.is_file():
        print_status_block(
            "Metadata extraction skipped - clean output already exists",
            f"Publication text: {publication_full_text_path}",
            f"Clean output: {clean_output_path}",
            f"Output with evidence: {with_evidence_output_path}",
        )
        return MavedbMetadataSchema.model_validate_json(clean_output_path.read_text(encoding='utf-8'))

    prompt = build_metadata_extraction_prompt(
        publication_full_text_path=publication_full_text_path,
        prompt_context=prompt_context,
    )
    _append_log_line(
        f"Built metadata extraction prompt for {publication_full_text_path}; characters: {len(prompt)}; output suffix: '{output_suffix or '[default]'}'"
    )

    client = instructor.from_provider(
        model_name,
        location='global',
        vertexai=True,
    )
    extraction_response = client.create(
        response_model=extraction_schema,
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
    )
    response = to_final_metadata_schema(extraction_response)
    with_evidence_output_path, clean_output_path = save_metadata_outputs(
        publication_full_text_path=publication_full_text_path,
        extraction_result_with_evidence=extraction_response,
        clean_result=response,
        output_dir=output_dir,
        suffix=output_suffix,
        prompt_context=prompt_context,
    )
    extracted_field_count = sum(
        field_value is not None for field_value in response.model_dump().values()
    )
    print_status_block(
        "Metadata extraction complete",
        f"Publication text: {publication_full_text_path}",
        f"Extracted non-null fields: {extracted_field_count}",
        f"Output with evidence: {with_evidence_output_path}",
        f"Clean output: {clean_output_path}",
        f"Log file: {METADATA_EXTRACTION_LOG_FILE}",
    )
    return response

def extract_metadata_from_publication(
    publication_full_text_path: str | Path,
    extraction_schema: Type[BaseModel],
    output_dir: str | Path = EXTRACTED_METADATA_OUTPUT_DIR,
    overwrite: bool = False,
    model_name: str = DEFAULT_LLM_MODEL_NAME,
) -> None:
    """Main function to extract metadata from a publication given the path to its full text.
    The extracted metadata will be saved in two formats: one with evidence quotes and one clean version without evidence. Both will be saved as JSON files in the specified output directory.
    If multiple distinct supplementary MaveDB metadata contexts map to the same publication,
    a separate prompt and output pair will be generated for each distinct context.
    Parameters:
    - publication_full_text_path: Path to the full text of the publication (e.g., a markdown file).
    - extraction_schema: Pydantic schema passed to instructor as the response model.
    - output_dir: Directory where the extracted metadata JSON files will be saved.
    - overwrite: Whether to overwrite existing extracted metadata outputs.
    - model_name: LLM model ID to use for metadata extraction.
    """
    publication_full_text_path = Path(publication_full_text_path).resolve()
    output_dir = Path(output_dir).resolve()

    prompt_contexts = _build_mavedb_prompt_contexts(publication_full_text_path)
    context_count = len(prompt_contexts) if prompt_contexts else 1

    print_status_block(
        "Starting metadata extraction",
        f"Publication text: {publication_full_text_path}",
        f"Output directory: {output_dir}",
        f"Model: {model_name}",
        f"Extraction schema: {extraction_schema.__name__}",
        f"Overwrite: {overwrite}",
        f"Matched MaveDB prompt contexts: {context_count}",
    )

    try:
        if not prompt_contexts:
            _extract_metadata_for_prompt_context(
                publication_full_text_path=publication_full_text_path,
                output_dir=output_dir,
                overwrite=overwrite,
                model_name=model_name,
                extraction_schema=extraction_schema,
                prompt_context=None,
                output_suffix="",
            )
            return

        for context_index, prompt_context in enumerate(prompt_contexts, start=1):
            output_suffix = _build_context_output_suffix(
                prompt_context=prompt_context,
                context_index=context_index,
                total_contexts=len(prompt_contexts),
            )
            _append_log_line(
                "Resolved supplementary MaveDB metadata context"
                f" for {publication_full_text_path}; context {context_index}/{len(prompt_contexts)};"
                f" source URNs: {', '.join(prompt_context['source_urns'])};"
                f" output suffix: '{output_suffix or '[default]'}'"
            )
            _extract_metadata_for_prompt_context(
                publication_full_text_path=publication_full_text_path,
                output_dir=output_dir,
                overwrite=overwrite,
                model_name=model_name,
                extraction_schema=extraction_schema,
                prompt_context=prompt_context,
                output_suffix=output_suffix,
            )

        if len(prompt_contexts) > 1:
            print_status_block(
                "Multiple metadata extraction contexts processed",
                f"Publication text: {publication_full_text_path}",
                f"Distinct MaveDB prompt contexts: {len(prompt_contexts)}",
                "Separate output files were written with context suffixes.",
            )
        return
    except Exception as exc:
        print_status_block(
            "Metadata extraction failed",
            f"Publication text: {publication_full_text_path}",
            f"Output directory: {output_dir}",
            f"Error: {exc}",
            f"Traceback: {traceback.format_exc()}"
        )
        raise


def bulk_extract_metadata_from_publications(
    publication_full_text_paths: list[str | Path],
    output_dir: str | Path = EXTRACTED_METADATA_OUTPUT_DIR,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
    overwrite: bool = False,
    model_name: str = DEFAULT_LLM_MODEL_NAME,
    extraction_schema: Type[BaseModel] = MavedbMetadataExtractionSchema,
    create_csv: bool = False,
) -> list[Path]:
    """Extract metadata from multiple publication full-text files in parallel.

    Parameters:
    - publication_full_text_paths: Paths to publication full-text markdown files.
    - output_dir: Directory where extracted metadata JSON files will be saved.
    - max_workers: Number of worker threads used to process publications concurrently.
    - overwrite: Whether to overwrite existing extracted metadata outputs.
    - model_name: LLM model ID to use for metadata extraction.
    - extraction_schema: Pydantic schema passed to instructor as the response model.
    - create_csv: Whether to create a single CSV from the curated clean metadata JSON outputs.

    Returns:
        - A list of publication paths successfully processed, in the same order as the input
            list. Publications that fail to process are omitted from the result.
    """
    output_dir = Path(output_dir).resolve()
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")

    publication_paths, excluded_paths = _filter_bulk_publication_paths(publication_full_text_paths)
    print_status_block(
        "Starting bulk metadata extraction",
        f"Publications queued: {len(publication_paths)}",
        f"Output directory: {output_dir}",
        f"Model: {model_name}",
        f"Max workers: {max_workers}",
        f"Overwrite: {overwrite}",
        f"Excluded publications: {len(excluded_paths)}",
        f"Log file: {METADATA_EXTRACTION_LOG_FILE}",
    )
    for excluded_path in excluded_paths:
        _append_log_line(
            f"Bulk metadata extraction excluded publication: {excluded_path}"
        )

    extracted_by_index: dict[int, Path] = {}
    completed_publications = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {
            executor.submit(
                extract_metadata_from_publication,
                publication_full_text_path=publication_path,
                output_dir=output_dir,
                overwrite=overwrite,
                model_name=model_name,
                extraction_schema=extraction_schema,
            ): index
            for index, publication_path in enumerate(publication_paths)
        }

        for future in tqdm(
            as_completed(future_to_index),
            total=len(future_to_index),
            desc="Extracting publication metadata",
            unit="publication",
        ):
            index = future_to_index[future]
            publication_path = publication_paths[index]
            try:
                future.result()
                extracted_by_index[index] = publication_path
                completed_publications += 1
                _append_log_line(
                    f"Bulk metadata extraction progress: {completed_publications}/{len(publication_paths)} publications processed; latest file: {publication_path}; status: ok"
                )
            except Exception as exc:
                completed_publications += 1
                _append_log_line(
                    f"Bulk metadata extraction progress: {completed_publications}/{len(publication_paths)} publications processed; latest file: {publication_path}; status: error; error: {exc}"
                )
                print_status_block(
                    "Bulk metadata extraction item failed",
                    f"Publication text: {publication_path}",
                    f"Error: {exc}",
                )

    extracted_publication_paths = [
        extracted_by_index[index] for index in sorted(extracted_by_index)
    ]
    print_status_block(
        "Bulk metadata extraction complete",
        f"Publications queued: {len(publication_paths)}",
        f"Publications extracted: {len(extracted_publication_paths)}",
        f"Publications failed: {len(publication_paths) - len(extracted_publication_paths)}",
        f"Output directory: {output_dir}",
        f"Log file: {METADATA_EXTRACTION_LOG_FILE}",
    )

    if create_csv:
        create_csv_from_curated_metadata_json(
            input_dir=output_dir / "clean",
            output_csv_path=output_dir / "clean_metadata.csv",
        )

    return extracted_publication_paths
    
    
    
if __name__ == "__main__":
    args = build_parser().parse_args()
    if args.max_workers < 1:
        raise ValueError("--max-workers must be at least 1")

    bulk_extract_metadata_from_publications(
        publication_full_text_paths=[
            str(path) for path in PUBLICATION_FULL_TEXT_MD_DIR.glob("*.md")
        ],
        output_dir=EXTRACTED_METADATA_OUTPUT_DIR,
        max_workers=args.max_workers,
        overwrite=args.overwrite,
        model_name=args.llm_model,
        create_csv=args.create_csv,
    )