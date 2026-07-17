import csv
import json
import os
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Type

import instructor
from pydantic import BaseModel
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


def get_evidence_output_path(
    publication_full_text_path: str | Path,
    output_dir: str | Path,
    suffix: str = "",
) -> Path:
    """Return the Step 1 evidence JSON output path for a publication."""
    publication_full_text_path = Path(publication_full_text_path).resolve()
    publication_stem = publication_full_text_path.stem
    output_dir = Path(output_dir).resolve()
    return output_dir / "step1_evidence" / f"{publication_stem}{suffix}.json"


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
    """
    Flatten curated evidence JSON outputs into a single CSV file.
    ---
    Parameters:
        input_dir: Directory containing evidence JSON files.
        output_csv_path: Path to write the resulting CSV file.
        log_file: Optional path to a log file for status messages.
    Returns:
        Path to the created CSV file.
    """
    input_dir = Path(input_dir).resolve()
    output_csv_path = Path(output_csv_path).resolve()

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Evidence directory not found: {input_dir}")

    json_paths = sorted(input_dir.glob("*.json"))
    if not json_paths:
        raise FileNotFoundError(f"No evidence JSON files found in: {input_dir}")
    if not output_csv_path:
        raise ValueError("Output CSV path must be specified.")
    if output_csv_path.is_dir():
        raise ValueError(f"Output CSV path must be a file, not a directory: {output_csv_path}")

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
            "Step 1 Evidence CSV created",
            f"Input directory: {input_dir}",
            f"JSON files processed: {len(json_paths)}",
            f"CSV output: {output_csv_path}",
        )
    return output_csv_path


def save_evidence_outputs(
    publication_full_text_path: str | Path,
    evidence_result: BaseModel,
    output_dir: str | Path,
    log_file: str | Path,
    suffix: str = "",
    prompt_context: PromptContext | None = None,
    output_metadata_builder=build_default_output_metadata,
) -> Path:
    """Persist Step 1 evidence extraction results."""
    publication_full_text_path = Path(publication_full_text_path).resolve()
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    
    output_path = get_evidence_output_path(
        publication_full_text_path=publication_full_text_path,
        output_dir=output_dir,
        suffix=suffix,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    output_metadata = output_metadata_builder(prompt_context)
    evidence_payload = evidence_result.model_dump()
    evidence_payload.update(output_metadata)
    
    output_path.write_text(
        json.dumps(evidence_payload, indent=JSON_INDENT),
        encoding="utf-8",
    )
    print_status_block(
        log_file,
        "Step 1 Evidence outputs saved",
        f"Publication text: {publication_full_text_path}",
        f"Step 1 Evidence output: {output_path}",
    )
    return output_path


def build_metadata_extraction_prompt(
    publication_full_text_path: str | Path,
    prompt_template_file: str | Path,
    log_file: str | Path,
    prompt_context: PromptContext | None = None,
    prompt_context_formatter=format_prompt_context_as_json,
) -> str:
    """Render the metadata/evidence extraction prompt for a publication and optional context."""
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


def _extract_evidence_for_prompt_context(
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
    """Extract evidence for one publication under a single prompt context (Step 1)."""
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    output_path = get_evidence_output_path(
        publication_full_text_path=publication_full_text_path,
        output_dir=output_dir,
        suffix=output_suffix,
    )

    if not overwrite and output_path.is_file():
        print_status_block(
            log_file,
            "Step 1 evidence extraction skipped - output already exists",
            f"Publication text: {publication_full_text_path}",
            f"Evidence output: {output_path}",
        )
        return extraction_schema.model_validate_json(
            output_path.read_text(encoding="utf-8")
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
        f"Built Step 1 evidence extraction prompt for {publication_full_text_path}; characters: {len(prompt)}; output suffix: '{output_suffix or '[default]'}'",
    )

    client = instructor.from_provider(
        model_name,
        location="global",
        vertexai=True,
    )
    extraction_response = client.create(
        response_model=extraction_schema,
        messages=[{"role": "user", "content": prompt}],
        thinking_config={
            "thinking_level": "high",
        },
    )
    
    save_evidence_outputs(
        publication_full_text_path=publication_full_text_path,
        evidence_result=extraction_response,
        output_dir=output_dir,
        log_file=log_file,
        suffix=output_suffix,
        prompt_context=prompt_context,
        output_metadata_builder=output_metadata_builder,
    )
    extracted_field_count = sum(
        field_value is not None for field_value in extraction_response.model_dump().values()
    )
    print_status_block(
        log_file,
        "Step 1 evidence extraction complete",
        f"Publication text: {publication_full_text_path}",
        f"Extracted evidence fields: {extracted_field_count}",
        f"Output: {output_path}",
        f"Log file: {log_file}",
    )
    return extraction_response


def extract_evidence_from_publication(
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
    """Extract evidence for one publication and write the resulting JSON output (Step 1)."""
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
        "Starting Step 1 evidence extraction",
        f"Publication text: {publication_full_text_path}",
        f"Output directory: {output_dir}",
        f"Model: {model_name}",
        f"Extraction schema: {extraction_schema.__name__}",
        f"Overwrite: {overwrite}",
        f"Matched prompt contexts: {context_count}",
    )

    try:
        if not prompt_contexts:
            _extract_evidence_for_prompt_context(
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
            _extract_evidence_for_prompt_context(
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
                "Multiple evidence extraction contexts processed",
                f"Publication text: {publication_full_text_path}",
                f"Distinct prompt contexts: {len(prompt_contexts)}",
                "Separate output files were written with context suffixes.",
            )
    except Exception as exc:
        print_status_block(
            log_file,
            "Step 1 evidence extraction failed",
            f"Publication text: {publication_full_text_path}",
            f"Output directory: {output_dir}",
            f"Error: {exc}",
            f"Traceback: {traceback.format_exc()}",
        )
        raise


run_step1_evidence_extraction = extract_evidence_from_publication


def bulk_extract_evidence_from_publications(
    publication_full_text_paths: list[str | Path],
    extraction_schema: Type[BaseModel],
    output_dir: str | Path,
    log_file: str | Path,
    prompt_template_file: str | Path,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
    overwrite: bool = False,
    model_name: str = DEFAULT_LLM_MODEL_NAME,
    create_csv: bool = True,
    excluded_publication_files: set[str] | frozenset[str] | None = None,
    prompt_context_builder=None,
    prompt_context_formatter=format_prompt_context_as_json,
    context_output_suffix_builder=build_default_context_output_suffix,
    output_metadata_builder=build_default_output_metadata,
) -> list[Path]:
    """Extract evidence for many publication files in parallel."""
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
        "Starting bulk evidence extraction",
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
            f"Bulk evidence extraction excluded publication: {excluded_path}",
        )

    extracted_by_index: dict[int, Path] = {}
    completed_publications = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {
            executor.submit(
                extract_evidence_from_publication,
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
            desc="Extracting publication evidence",
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
                    f"Bulk evidence extraction progress: {completed_publications}/{len(publication_paths)} publications processed; latest file: {publication_path}; status: ok",
                )
            except Exception as exc:
                completed_publications += 1
                append_log_line(
                    log_file,
                    f"Bulk evidence extraction progress: {completed_publications}/{len(publication_paths)} publications processed; latest file: {publication_path}; status: error; error: {exc}",
                )
                print_status_block(
                    log_file,
                    "Bulk evidence extraction item failed",
                    f"Publication text: {publication_path}",
                    f"Error: {exc}",
                )

    extracted_publication_paths = [
        extracted_by_index[index] for index in sorted(extracted_by_index)
    ]
    print_status_block(
        log_file,
        "Bulk evidence extraction complete",
        f"Publications queued: {len(publication_paths)}",
        f"Publications extracted: {len(extracted_publication_paths)}",
        f"Publications failed: {len(publication_paths) - len(extracted_publication_paths)}",
        f"Output directory: {output_dir}",
        f"Log file: {log_file}",
    )

    if create_csv:
        create_csv_from_curated_metadata_json(
            input_dir=output_dir / "step1_evidence",
            output_csv_path=output_dir / "step1_evidence.csv",
            log_file=log_file,
        )

    return extracted_publication_paths

