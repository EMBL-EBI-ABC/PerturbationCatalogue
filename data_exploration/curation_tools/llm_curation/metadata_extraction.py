import csv
import json
import os
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Type

import instructor
from pydantic import BaseModel, create_model
from tqdm import tqdm

from curation_tools.llm_curation.logging_utils import (
    _ensure_log_file,
    append_log_line,
    print_status_block,
)


DEFAULT_LLM_MODEL_NAME = os.getenv("LLM_MODEL_NAME", "google/gemini-3.5-flash")
DEFAULT_DOWNLOAD_MAX_WORKERS = min(32, (os.cpu_count() or 1) * 4)
JSON_INDENT = 2

PromptContext = dict[str, object]


def format_prompt_context_as_json(prompt_context: PromptContext | None) -> str:
    """Render prompt context as JSON with a generic fallback message."""
    if not prompt_context:
        return "No supplementary metadata was provided for this publication."
    return json.dumps(prompt_context, indent=JSON_INDENT, ensure_ascii=True)


def build_default_context_output_suffix(
    prompt_context: PromptContext | None,
    context_index: int,
    total_contexts: int,
) -> str:
    """Build a context suffix without assuming a domain-specific identifier."""
    del prompt_context
    return "" if total_contexts == 1 else f"__ctx_{context_index:02d}"


def build_default_output_metadata(prompt_context: PromptContext | None) -> dict[str, object]:
    """Return no extra metadata unless the caller provides a custom builder."""
    del prompt_context
    return {}


def _filter_bulk_publication_paths(
    publication_full_text_paths: list[str | Path],
    excluded_publication_files: set[str] | frozenset[str] | None = None,
) -> tuple[list[Path], list[Path]]:
    """Split publication paths into included and excluded sets for bulk processing."""
    publication_paths = [Path(path).resolve() for path in publication_full_text_paths]
    excluded_paths: list[Path] = []
    included_paths: list[Path] = []
    excluded_publication_files = set(excluded_publication_files or ())
    for publication_path in publication_paths:
        if publication_path.name in excluded_publication_files:
            excluded_paths.append(publication_path)
            continue
        included_paths.append(publication_path)
    return included_paths, excluded_paths


def get_metadata_output_paths(
    publication_full_text_path: str | Path,
    output_dir: str | Path,
    suffix: str = "",
) -> tuple[Path, Path]:
    """Return the with-evidence and clean JSON output paths for a publication."""
    publication_full_text_path = Path(publication_full_text_path).resolve()
    publication_stem = publication_full_text_path.stem
    output_dir = Path(output_dir).resolve()
    with_evidence_path = output_dir / "with_evidence" / f"{publication_stem}{suffix}.json"
    clean_path = output_dir / "clean" / f"{publication_stem}{suffix}.json"
    return with_evidence_path, clean_path


def _format_metadata_csv_cell(value: object) -> str:
    """Serialize a metadata value into a CSV-safe string representation."""
    if value is None:
        return ""
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=True, sort_keys=True)
    return str(value)


def create_csv_from_curated_metadata_json(
    input_dir: str | Path,
    output_csv_path: str | Path,
    log_file: str | Path | None = None,
) -> Path:
    """Flatten curated clean-metadata JSON outputs into a single CSV file."""
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

    if log_file is not None:
        print_status_block(
            log_file,
            "Curated metadata CSV created",
            f"Input directory: {input_dir}",
            f"JSON files processed: {len(json_paths)}",
            f"CSV output: {output_csv_path}",
        )
    return output_csv_path


def save_metadata_outputs(
    publication_full_text_path: str | Path,
    extraction_result_with_evidence: BaseModel,
    clean_result: BaseModel,
    output_dir: str | Path,
    log_file: str | Path,
    suffix: str = "",
    prompt_context: PromptContext | None = None,
    output_metadata_builder=build_default_output_metadata,
) -> tuple[Path, Path]:
    """Persist metadata extraction results in both evidence-rich and clean forms."""
    publication_full_text_path = Path(publication_full_text_path).resolve()
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    with_evidence_path, clean_path = get_metadata_output_paths(
        publication_full_text_path=publication_full_text_path,
        output_dir=output_dir,
        suffix=suffix,
    )
    with_evidence_path.parent.mkdir(parents=True, exist_ok=True)
    clean_path.parent.mkdir(parents=True, exist_ok=True)

    output_metadata = output_metadata_builder(prompt_context)
    with_evidence_payload = extraction_result_with_evidence.model_dump()
    clean_payload = clean_result.model_dump()
    with_evidence_payload.update(output_metadata)
    clean_payload.update(output_metadata)

    with_evidence_path.write_text(
        json.dumps(with_evidence_payload, indent=JSON_INDENT),
        encoding="utf-8",
    )
    clean_path.write_text(
        json.dumps(clean_payload, indent=JSON_INDENT),
        encoding="utf-8",
    )
    print_status_block(
        log_file,
        "Metadata outputs saved",
        f"Publication text: {publication_full_text_path}",
        f"With evidence: {with_evidence_path}",
        f"Clean output: {clean_path}",
    )
    return with_evidence_path, clean_path


def to_final_metadata_schema(
    extraction_result_with_evidence: BaseModel,
    extraction_schema: Type[BaseModel],
) -> BaseModel:
    """Drop evidence fields and validate the remaining payload against the final schema."""
    final_payload = {
        field_name: field_value
        for field_name, field_value in extraction_result_with_evidence.model_dump().items()
        if not field_name.endswith("_evidence")
    }
    if not any(
        field_name.endswith("_evidence")
        for field_name in extraction_schema.model_fields
    ):
        return extraction_schema.model_validate(final_payload)

    no_evidence_schema = create_model(
        f"{extraction_schema.__name__}NoEvidence",
        **{
            field_name: (field_info.annotation, field_info)
            for field_name, field_info in extraction_schema.model_fields.items()
            if not field_name.endswith("_evidence")
        },
    )
    return no_evidence_schema.model_validate(final_payload)


def build_metadata_extraction_prompt(
    publication_full_text_path: str | Path,
    prompt_template_file: str | Path,
    log_file: str | Path,
    prompt_context: PromptContext | None = None,
    prompt_context_formatter=format_prompt_context_as_json,
) -> str:
    """Render the metadata extraction prompt for a publication and optional context."""
    publication_full_text_path = Path(publication_full_text_path).resolve()
    prompt_template_file = Path(prompt_template_file).resolve()
    log_file = _ensure_log_file(log_file)
    publication_full_text = publication_full_text_path.read_text(encoding="utf-8")
    append_log_line(
        log_file,
        f"Loaded publication text from {publication_full_text_path}; characters: {len(publication_full_text)}",
    )
    prompt_template = prompt_template_file.read_text(encoding="utf-8")
    supplementary_metadata = prompt_context_formatter(prompt_context)
    return prompt_template.format(
        supplementary_metadata=supplementary_metadata,
        supplementary_mavedb_metadata=supplementary_metadata,
        publication_full_text=publication_full_text,
    )


def _extract_metadata_for_prompt_context(
    publication_full_text_path: Path,
    output_dir: str | Path,
    log_file: str | Path,
    prompt_template_file: str | Path,
    overwrite: bool,
    model_name: str,
    extraction_schema: Type[BaseModel],
    prompt_context: PromptContext | None,
    output_suffix: str,
    prompt_context_formatter,
    output_metadata_builder,
) -> BaseModel:
    """Extract metadata for one publication under a single prompt context."""
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    with_evidence_output_path, clean_output_path = get_metadata_output_paths(
        publication_full_text_path=publication_full_text_path,
        output_dir=output_dir,
        suffix=output_suffix,
    )

    if not overwrite and clean_output_path.is_file():
        print_status_block(
            log_file,
            "Metadata extraction skipped - clean output already exists",
            f"Publication text: {publication_full_text_path}",
            f"Clean output: {clean_output_path}",
            f"Output with evidence: {with_evidence_output_path}",
        )
        return extraction_schema.model_validate_json(
            clean_output_path.read_text(encoding="utf-8")
        )

    prompt = build_metadata_extraction_prompt(
        publication_full_text_path=publication_full_text_path,
        prompt_template_file=prompt_template_file,
        log_file=log_file,
        prompt_context=prompt_context,
        prompt_context_formatter=prompt_context_formatter,
    )
    append_log_line(
        log_file,
        f"Built metadata extraction prompt for {publication_full_text_path}; characters: {len(prompt)}; output suffix: '{output_suffix or '[default]'}'",
    )

    client = instructor.from_provider(
        model_name,
        location="global",
        vertexai=True,
    )
    extraction_response = client.create(
        response_model=extraction_schema,
        messages=[{"role": "user", "content": prompt}],
    )
    response = to_final_metadata_schema(extraction_response, extraction_schema)
    with_evidence_output_path, clean_output_path = save_metadata_outputs(
        publication_full_text_path=publication_full_text_path,
        extraction_result_with_evidence=extraction_response,
        clean_result=response,
        output_dir=output_dir,
        log_file=log_file,
        suffix=output_suffix,
        prompt_context=prompt_context,
        output_metadata_builder=output_metadata_builder,
    )
    extracted_field_count = sum(
        field_value is not None for field_value in response.model_dump().values()
    )
    other_field_count = sum(
        field_value == "Other" for field_value in response.model_dump().values()
    )
    print_status_block(
        log_file,
        "Metadata extraction complete",
        f"Publication text: {publication_full_text_path}",
        f"Extracted non-null fields: {extracted_field_count}",
        f"Extracted 'Other' fields: {other_field_count}",
        f"Output with evidence: {with_evidence_output_path}",
        f"Clean output: {clean_output_path}",
        f"Log file: {log_file}",
    )
    return response


def extract_metadata_from_publication(
    publication_full_text_path: str | Path,
    extraction_schema: Type[BaseModel],
    output_dir: str | Path,
    log_file: str | Path,
    prompt_template_file: str | Path,
    overwrite: bool = False,
    model_name: str = DEFAULT_LLM_MODEL_NAME,
    prompt_context_builder=None,
    prompt_context_formatter=format_prompt_context_as_json,
    context_output_suffix_builder=build_default_context_output_suffix,
    output_metadata_builder=build_default_output_metadata,
) -> None:
    """Extract metadata for one publication and write the resulting JSON outputs."""
    publication_full_text_path = Path(publication_full_text_path).resolve()
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    prompt_template_file = Path(prompt_template_file).resolve()
    prompt_contexts = (
        prompt_context_builder(publication_full_text_path)
        if prompt_context_builder is not None
        else []
    )
    context_count = len(prompt_contexts) if prompt_contexts else 1

    print_status_block(
        log_file,
        "Starting metadata extraction",
        f"Publication text: {publication_full_text_path}",
        f"Output directory: {output_dir}",
        f"Model: {model_name}",
        f"Extraction schema: {extraction_schema.__name__}",
        f"Overwrite: {overwrite}",
        f"Matched prompt contexts: {context_count}",
    )

    try:
        if not prompt_contexts:
            _extract_metadata_for_prompt_context(
                publication_full_text_path=publication_full_text_path,
                output_dir=output_dir,
                log_file=log_file,
                prompt_template_file=prompt_template_file,
                overwrite=overwrite,
                model_name=model_name,
                extraction_schema=extraction_schema,
                prompt_context=None,
                output_suffix="",
                prompt_context_formatter=prompt_context_formatter,
                output_metadata_builder=output_metadata_builder,
            )
            return

        for context_index, prompt_context in enumerate(prompt_contexts, start=1):
            output_suffix = context_output_suffix_builder(
                prompt_context,
                context_index,
                len(prompt_contexts),
            )
            append_log_line(
                log_file,
                "Resolved supplementary metadata context"
                f" for {publication_full_text_path}; context {context_index}/{len(prompt_contexts)};"
                f" output suffix: '{output_suffix or '[default]'}'",
            )
            _extract_metadata_for_prompt_context(
                publication_full_text_path=publication_full_text_path,
                output_dir=output_dir,
                log_file=log_file,
                prompt_template_file=prompt_template_file,
                overwrite=overwrite,
                model_name=model_name,
                extraction_schema=extraction_schema,
                prompt_context=prompt_context,
                output_suffix=output_suffix,
                prompt_context_formatter=prompt_context_formatter,
                output_metadata_builder=output_metadata_builder,
            )

        if len(prompt_contexts) > 1:
            print_status_block(
                log_file,
                "Multiple metadata extraction contexts processed",
                f"Publication text: {publication_full_text_path}",
                f"Distinct prompt contexts: {len(prompt_contexts)}",
                "Separate output files were written with context suffixes.",
            )
    except Exception as exc:
        print_status_block(
            log_file,
            "Metadata extraction failed",
            f"Publication text: {publication_full_text_path}",
            f"Output directory: {output_dir}",
            f"Error: {exc}",
            f"Traceback: {traceback.format_exc()}",
        )
        raise


def bulk_extract_metadata_from_publications(
    publication_full_text_paths: list[str | Path],
    extraction_schema: Type[BaseModel],
    output_dir: str | Path,
    log_file: str | Path,
    prompt_template_file: str | Path,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
    overwrite: bool = False,
    model_name: str = DEFAULT_LLM_MODEL_NAME,
    create_csv: bool = False,
    excluded_publication_files: set[str] | frozenset[str] | None = None,
    prompt_context_builder=None,
    prompt_context_formatter=format_prompt_context_as_json,
    context_output_suffix_builder=build_default_context_output_suffix,
    output_metadata_builder=build_default_output_metadata,
) -> list[Path]:
    """Extract metadata for many publication files in parallel."""
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    prompt_template_file = Path(prompt_template_file).resolve()

    publication_paths, excluded_paths = _filter_bulk_publication_paths(
        publication_full_text_paths,
        excluded_publication_files=excluded_publication_files,
    )
    print_status_block(
        log_file,
        "Starting bulk metadata extraction",
        f"Publications queued: {len(publication_paths)}",
        f"Output directory: {output_dir}",
        f"Model: {model_name}",
        f"Max workers: {max_workers}",
        f"Overwrite: {overwrite}",
        f"Excluded publications: {len(excluded_paths)}",
        f"Log file: {log_file}",
    )
    for excluded_path in excluded_paths:
        append_log_line(
            log_file,
            f"Bulk metadata extraction excluded publication: {excluded_path}",
        )

    extracted_by_index: dict[int, Path] = {}
    completed_publications = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {
            executor.submit(
                extract_metadata_from_publication,
                publication_full_text_path=publication_path,
                output_dir=output_dir,
                log_file=log_file,
                prompt_template_file=prompt_template_file,
                overwrite=overwrite,
                model_name=model_name,
                extraction_schema=extraction_schema,
                prompt_context_builder=prompt_context_builder,
                prompt_context_formatter=prompt_context_formatter,
                context_output_suffix_builder=context_output_suffix_builder,
                output_metadata_builder=output_metadata_builder,
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
                append_log_line(
                    log_file,
                    f"Bulk metadata extraction progress: {completed_publications}/{len(publication_paths)} publications processed; latest file: {publication_path}; status: ok",
                )
            except Exception as exc:
                completed_publications += 1
                append_log_line(
                    log_file,
                    f"Bulk metadata extraction progress: {completed_publications}/{len(publication_paths)} publications processed; latest file: {publication_path}; status: error; error: {exc}",
                )
                print_status_block(
                    log_file,
                    "Bulk metadata extraction item failed",
                    f"Publication text: {publication_path}",
                    f"Error: {exc}",
                )

    extracted_publication_paths = [
        extracted_by_index[index] for index in sorted(extracted_by_index)
    ]
    print_status_block(
        log_file,
        "Bulk metadata extraction complete",
        f"Publications queued: {len(publication_paths)}",
        f"Publications extracted: {len(extracted_publication_paths)}",
        f"Publications failed: {len(publication_paths) - len(extracted_publication_paths)}",
        f"Output directory: {output_dir}",
        f"Log file: {log_file}",
    )

    if create_csv:
        create_csv_from_curated_metadata_json(
            input_dir=output_dir / "clean",
            output_csv_path=output_dir / "clean_metadata.csv",
            log_file=log_file,
        )

    return extracted_publication_paths