"""Pipeline for mapping verbatim Step 1 evidence to controlled ontology vocabularies (Step 2)."""

import argparse
import json
import os
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Type

import instructor
from pydantic import BaseModel
from tqdm import tqdm

from curation_tools.llm_curation.logging_utils import (
    _ensure_log_file,
    append_log_line,
    print_status_block,
)
from curation_tools.llm_curation.llm_curation_schema import SpecificTermExtractionSchema
from curation_tools.llm_curation.metadata_extraction import (
    create_csv_from_curated_metadata_json,
)

DEFAULT_LLM_MODEL_NAME = os.getenv("LLM_MODEL_NAME", "google/gemini-3.5-flash")
DEFAULT_CONCURRENCY_WORKERS = min(64, (os.cpu_count() or 1) * 8)
JSON_INDENT = 2


def get_normalized_output_path(
    step1_evidence_path: str | Path,
    output_dir: str | Path,
) -> Path:
    """Return the Step 2 normalized JSON output path for a Step 1 evidence artifact."""
    step1_evidence_path = Path(step1_evidence_path).resolve()
    filename_stem = step1_evidence_path.stem
    output_dir = Path(output_dir).resolve()
    return output_dir / "step2_normalized" / f"{filename_stem}.json"


def load_and_format_mave_metadata(
    evidence_payload: dict[str, Any],
    mavedb_metadata_dir: Path | str | None = None,
) -> str:
    """Format matching MaveDB metadata from Step 1 payload source urns/files."""
    if not evidence_payload or not evidence_payload.get("__source_urns"):
        return "No supplementary MaveDB metadata was available for this publication."

    if not mavedb_metadata_dir:
        try:
            from curation_tools.llm_curation.mavedb.processing import (
                MAVEDB_METADATA_OUTPUT_DIR,
            )

            mavedb_metadata_dir = MAVEDB_METADATA_OUTPUT_DIR
        except ImportError:
            pass

    if not mavedb_metadata_dir:
        return "No supplementary MaveDB metadata was available for this publication."

    mavedb_metadata_dir = Path(mavedb_metadata_dir).resolve()
    source_files = evidence_payload.get("__source_files", [])
    source_urns = evidence_payload.get("__source_urns", [])

    merged_metadata: dict[str, Any] = {}
    for filename in source_files:
        filepath = mavedb_metadata_dir / filename
        if filepath.is_file():
            try:
                entry_payload = json.loads(filepath.read_text(encoding="utf-8"))
                from curation_tools.llm_curation.mavedb.processing import (
                    extract_curated_mavedb_prompt_metadata,
                    merge_prompt_metadata_value,
                )

                curated_metadata = extract_curated_mavedb_prompt_metadata(entry_payload)
                for field_name, field_value in curated_metadata.items():
                    merged_metadata[field_name] = merge_prompt_metadata_value(
                        merged_metadata.get(field_name), field_value
                    )
            except Exception:
                pass

    if not merged_metadata:
        return "No supplementary MaveDB metadata was available for this publication."

    return json.dumps(
        {
            "source_urns": source_urns,
            "source_files": source_files,
            "curated_mavedb_metadata": merged_metadata,
        },
        indent=JSON_INDENT,
        ensure_ascii=True,
    )


def normalize_single_evidence_artifact(
    step1_path: Path,
    output_dir: Path,
    log_file: Path,
    prompt_template_file: Path,
    mavedb_metadata_dir: Path | None,
    overwrite: bool,
    model_name: str,
    normalization_schema: Type[BaseModel] = SpecificTermExtractionSchema,
    verbose: bool = False,
) -> Path:
    """Normalize a single Step 1 evidence JSON artifact to controlled vocabularies."""
    output_path = get_normalized_output_path(step1_path, output_dir)

    if not overwrite and output_path.is_file():
        print_status_block(
            log_file,
            "Step 2 evidence normalization skipped - output already exists",
            f"Step 1 artifact: {step1_path}",
            f"Normalized output: {output_path}",
        )
        return output_path

    evidence_payload = json.loads(step1_path.read_text(encoding="utf-8"))

    # Prepare supplementary MaveDB metadata context
    supplementary_mavedb_metadata = load_and_format_mave_metadata(
        evidence_payload=evidence_payload,
        mavedb_metadata_dir=mavedb_metadata_dir,
    )

    # Clean the input evidence JSON to contain only relevant fields for prompt clarity
    cleaned_evidence = {
        k: v
        for k, v in evidence_payload.items()
        if k.endswith("_evidence") or k in ("dataset_id", "data_modality")
    }
    step1_evidence_json_str = json.dumps(cleaned_evidence, indent=JSON_INDENT)

    # Read and render Step 2 prompt template
    prompt_template = prompt_template_file.read_text(encoding="utf-8")
    prompt = prompt_template.format(
        step1_evidence=step1_evidence_json_str,
        supplementary_mavedb_metadata=supplementary_mavedb_metadata,
    )

    append_log_line(
        log_file,
        f"Built Step 2 normalization prompt for {step1_path}; characters: {len(prompt)}",
    )

    if verbose:
        print_status_block(
            log_file,
            "[VERBOSE] Step 2 full prompt",
            f"Step 1 artifact: {step1_path}",
            "----- PROMPT START -----",
            prompt,
            "----- PROMPT END -----",
        )

    client = instructor.from_provider(
        model_name,
        location="global",
        vertexai=True,
    )
    normalized_response = client.create(
        response_model=normalization_schema,
        messages=[{"role": "user", "content": prompt}],
        thinking_config={
            "thinking_level": "high",
        },
        generation_config={
            "temperature": 0.2,
        }
    )

    # Persist the output
    normalized_payload = normalized_response.model_dump()

    # Retain critical tracking fields from Step 1
    for tracking_field in ["__source_urns", "__source_files", "dataset_id"]:
        if tracking_field in evidence_payload:
            normalized_payload[tracking_field] = evidence_payload[tracking_field]

    # Add Curation Agent Metadata
    normalized_payload["curation_agent_type"] = "LLM"
    normalized_payload["curation_agent_name"] = model_name

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(normalized_payload, indent=JSON_INDENT),
        encoding="utf-8",
    )

    print_status_block(
        log_file,
        "Step 2 evidence normalization complete",
        f"Step 1 artifact: {step1_path}",
        f"Normalized output: {output_path}",
    )
    return output_path


def normalize_evidence_artifacts(
    step1_dir: str | Path,
    output_dir: str | Path,
    log_file: str | Path,
    prompt_template_file: str | Path,
    mavedb_metadata_dir: str | Path | None = None,
    max_workers: int = DEFAULT_CONCURRENCY_WORKERS,
    overwrite: bool = False,
    model_name: str = DEFAULT_LLM_MODEL_NAME,
    create_csv: bool = True,
    verbose: bool = False,
    normalization_schema: Type[BaseModel] = SpecificTermExtractionSchema,
) -> list[Path]:
    """Runner to concurrently map and normalize a directory of Step 1 JSON evidence files."""
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")

    step1_dir = Path(step1_dir).resolve()
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    prompt_template_file = Path(prompt_template_file).resolve()
    resolved_mavedb_metadata_dir = (
        Path(mavedb_metadata_dir).resolve() if mavedb_metadata_dir else None
    )

    if not step1_dir.is_dir():
        raise FileNotFoundError(f"Step 1 evidence directory not found: {step1_dir}")

    step1_paths = sorted(step1_dir.glob("*.json"))
    if not step1_paths:
        print_status_block(
            log_file,
            "Step 2 normalization warning - no files found",
            f"Directory searched: {step1_dir}",
        )
        return []

    print_status_block(
        log_file,
        "Starting Step 2 evidence normalization",
        f"Step 1 artifacts queued: {len(step1_paths)}",
        f"Output directory: {output_dir}",
        f"Model: {model_name}",
        f"Max workers: {max_workers}",
        f"Overwrite: {overwrite}",
        f"Log file: {log_file}",
    )

    normalized_outputs: list[Path] = []
    completed_files = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_path = {
            executor.submit(
                normalize_single_evidence_artifact,
                step1_path=step1_path,
                output_dir=output_dir,
                log_file=log_file,
                prompt_template_file=prompt_template_file,
                mavedb_metadata_dir=resolved_mavedb_metadata_dir,
                overwrite=overwrite,
                model_name=model_name,
                normalization_schema=normalization_schema,
                verbose=verbose,
            ): step1_path
            for step1_path in step1_paths
        }

        for future in tqdm(
            as_completed(future_to_path),
            total=len(future_to_path),
            desc="Normalizing evidence to ontology terms",
            unit="artifact",
        ):
            step1_path = future_to_path[future]
            try:
                output_path = future.result()
                normalized_outputs.append(output_path)
                completed_files += 1
                append_log_line(
                    log_file,
                    f"Step 2 normalization progress: {completed_files}/{len(step1_paths)} processed; latest file: {step1_path}; status: ok",
                )
            except Exception as exc:
                completed_files += 1
                append_log_line(
                    log_file,
                    f"Step 2 normalization progress: {completed_files}/{len(step1_paths)} processed; latest file: {step1_path}; status: error; error: {exc}",
                )
                print_status_block(
                    log_file,
                    "Step 2 normalization item failed",
                    f"Step 1 artifact: {step1_path}",
                    f"Error: {exc}",
                    f"Traceback: {traceback.format_exc()}",
                )

    print_status_block(
        log_file,
        "Step 2 evidence normalization batch complete",
        f"Artifacts queued: {len(step1_paths)}",
        f"Artifacts successfully normalized: {len(normalized_outputs)}",
        f"Artifacts failed: {len(step1_paths) - len(normalized_outputs)}",
        f"Output directory: {output_dir}",
        f"Log file: {log_file}",
    )

    if create_csv and normalized_outputs:
        try:
            create_csv_from_curated_metadata_json(
                input_dir=output_dir / "step2_normalized",
                output_csv_path=output_dir / "clean_metadata.csv",
                log_file=log_file,
            )
        except Exception as csv_exc:
            print_status_block(
                log_file,
                "Step 2 CSV creation failed",
                f"Error: {csv_exc}",
                f"Traceback: {traceback.format_exc()}",
            )

    return normalized_outputs


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for Step 2 specific term extraction."""
    parser = argparse.ArgumentParser(
        description="Normalize verbatim Step 1 evidence JSON files to controlled vocabularies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--step1-dir",
        type=Path,
        required=True,
        help="Directory containing Step 1 evidence JSON files to process.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory used to write Step 2 normalized metadata outputs.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        required=True,
        help="Log file used for progress and errors.",
    )
    parser.add_argument(
        "--prompt-template-file",
        type=Path,
        required=True,
        help="Prompt template file used to build normalization prompts.",
    )
    parser.add_argument(
        "--mavedb-metadata-dir",
        type=Path,
        default=None,
        help="Optional directory containing MaveDB metadata entries for lookup.",
    )
    parser.add_argument(
        "--llm-model",
        default=DEFAULT_LLM_MODEL_NAME,
        help="LLM model ID to use for metadata normalization.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_CONCURRENCY_WORKERS,
        help="Number of worker threads to use for bulk normalization.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing extracted metadata outputs.",
    )
    parser.add_argument(
        "--no-csv",
        action="store_false",
        dest="create_csv",
        help="Disable creation of a single CSV from the curated metadata JSON outputs.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print and log the full prompt, MaveDB metadata, and Step 1 evidence text for each item.",
    )
    return parser


def main() -> None:
    """Main CLI entry point for Step 2 evidence normalization."""
    args = build_parser().parse_args()
    normalize_evidence_artifacts(
        step1_dir=args.step1_dir,
        output_dir=args.output_dir,
        log_file=args.log_file,
        prompt_template_file=args.prompt_template_file,
        mavedb_metadata_dir=args.mavedb_metadata_dir,
        max_workers=args.max_workers,
        overwrite=args.overwrite,
        model_name=args.llm_model,
        create_csv=args.create_csv,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
