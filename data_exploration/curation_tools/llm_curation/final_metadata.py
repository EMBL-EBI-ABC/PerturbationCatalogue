"""Step 5: combine normalized LLM metadata with MaveDB metadata."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from curation_tools.llm_curation.logging_utils import (
    _ensure_log_file,
    append_log_line,
    print_status_block,
)
from curation_tools.llm_curation.mavedb.processing import (
    MAVEDB_METADATA_OUTPUT_DIR,
    extract_curated_mavedb_prompt_metadata,
    format_urn_for_filename,
    merge_prompt_metadata_value,
)
from curation_tools.llm_curation.metadata_extraction import (
    create_csv_from_curated_metadata_json,
)
from curation_tools.perturbseq_anndata_schema import ObsSchema

JSON_INDENT = 2
OBS_SCHEMA_FIELDS = tuple(ObsSchema.to_schema().columns.keys())
PROVENANCE_FIELDS = ("__source_urns", "__source_files")
EXCLUDED_INPUT_JSON_NAMES = {
    "step4_backfill_audit.json",
    "step3_ontology_candidates.json",
    "approved_ontology_terms.json",
}


def get_final_json_dir(output_dir: str | Path) -> Path:
    """Return the directory containing Step 5 JSON and CSV outputs."""
    output_dir = Path(output_dir).resolve()
    return (
        output_dir if output_dir.name == "step5_final" else output_dir / "step5_final"
    )


def get_final_csv_path(output_dir: str | Path) -> Path:
    """Return the editable Step 5 CSV path for an output directory."""
    return get_final_json_dir(output_dir) / "step5_final_metadata.csv"


def load_final_csv(csv_path: str | Path) -> pd.DataFrame:
    """Load a final metadata CSV without coercing empty cells to NaN."""
    return pd.read_csv(Path(csv_path).resolve(), dtype=str, keep_default_na=False)


def build_csv_cell_changes(
    original: pd.DataFrame, edited: pd.DataFrame
) -> list[dict[str, object]]:
    """Return cell-level differences between two editor snapshots."""
    if list(original.columns) != list(edited.columns):
        raise ValueError("CSV column names/order cannot be changed in the editor")
    if len(original) != len(edited):
        raise ValueError("CSV rows cannot be added or removed in the editor")

    changes: list[dict[str, object]] = []
    for row_index in range(len(original)):
        dataset_id = edited.iloc[row_index].get("dataset_id", "")
        for column in original.columns:
            old_value = str(original.iloc[row_index][column])
            new_value = str(edited.iloc[row_index][column])
            if old_value != new_value:
                changes.append(
                    {
                        "row_index": row_index,
                        "dataset_id": dataset_id,
                        "field_name": column,
                        "previous_value": old_value,
                        "new_value": new_value,
                    }
                )
    return changes


def save_final_csv_edits(
    csv_path: str | Path,
    original: pd.DataFrame,
    edited: pd.DataFrame,
) -> list[dict[str, object]]:
    """Persist cell edits and a compact audit trail beside the CSV."""
    csv_path = Path(csv_path).resolve()
    changes = build_csv_cell_changes(original, edited)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    edited.to_csv(csv_path, index=False)
    audit_path = csv_path.with_name("step5_csv_edit_audit.json")
    audit_path.write_text(json.dumps(changes, indent=JSON_INDENT), encoding="utf-8")
    return changes


def get_final_output_path(
    normalized_metadata_path: str | Path,
    output_dir: str | Path,
) -> Path:
    """Return the Step 5 JSON path for one normalized metadata artifact."""
    normalized_metadata_path = Path(normalized_metadata_path).resolve()
    return get_final_json_dir(output_dir) / normalized_metadata_path.name


def _is_missing(value: object) -> bool:
    return value is None or value == "" or value == [] or value == {}


def _unique_strings(values: list[object]) -> list[str]:
    result: list[str] = []
    for value in values:
        if isinstance(value, str) and value and value not in result:
            result.append(value)
    return result


def _join_values(values: list[object]) -> str | None:
    unique_values = _unique_strings(values)
    return "|".join(unique_values) if unique_values else None


def _load_source_entries(
    normalized_payload: dict[str, Any],
    mavedb_metadata_dir: Path,
) -> list[dict[str, Any]]:
    """Load the raw MaveDB entries referenced by a normalized artifact."""
    source_files = normalized_payload.get("__source_files") or []
    if not source_files:
        source_files = [
            f"{format_urn_for_filename(urn)}.json"
            for urn in normalized_payload.get("__source_urns", [])
        ]

    entries: list[dict[str, Any]] = []
    seen_paths: set[Path] = set()
    for source_file in source_files:
        source_path = mavedb_metadata_dir / Path(str(source_file)).name
        if source_path in seen_paths or not source_path.is_file():
            continue
        seen_paths.add(source_path)
        try:
            payload = json.loads(source_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(payload, dict):
            entries.append(payload)
    return entries


def _extract_mavedb_fields(entries: list[dict[str, Any]]) -> dict[str, object]:
    """Map useful MaveDB-only metadata into fields defined by ``ObsSchema``."""
    if not entries:
        return {}

    target_records: list[dict[str, Any]] = []
    merged_prompt_metadata: dict[str, object] = {}
    for entry in entries:
        for target in entry.get("targetGenes") or []:
            if isinstance(target, dict):
                target_records.append(target)
        prompt_metadata = extract_curated_mavedb_prompt_metadata(entry)
        for field_name, field_value in prompt_metadata.items():
            merged_prompt_metadata[field_name] = merge_prompt_metadata_value(
                merged_prompt_metadata.get(field_name), field_value
            )

    target_names = _unique_strings([target.get("name") for target in target_records])
    taxonomies: list[object] = []
    for target in target_records:
        sequence = target.get("targetSequence") or {}
        taxonomy = sequence.get("taxonomy") or {}
        taxonomies.append(taxonomy.get("organismName"))

    final_fields: dict[str, object] = {}
    joined_target_names = _join_values(target_names)
    if joined_target_names:
        final_fields["perturbed_target_symbol"] = joined_target_names
    if target_names:
        final_fields["perturbed_target_number"] = len(target_names)
        final_fields["number_of_perturbed_targets"] = str(len(target_names))

    total_variants = merged_prompt_metadata.get("total_variants")
    if total_variants is not None:
        final_fields["library_total_variants"] = total_variants

    publication_titles = merged_prompt_metadata.get("primary_publication_titles")
    publication_dois = merged_prompt_metadata.get("primary_publication_dois")
    publication_years = merged_prompt_metadata.get("primary_publication_years")
    if publication_titles:
        final_fields["study_title"] = _join_values(
            publication_titles
            if isinstance(publication_titles, list)
            else [publication_titles]
        )
    if publication_dois:
        final_fields["study_uri"] = _join_values(
            publication_dois
            if isinstance(publication_dois, list)
            else [publication_dois]
        )
    if isinstance(publication_years, list) and publication_years:
        final_fields["study_year"] = publication_years[0]

    experiment_title = merged_prompt_metadata.get("experiment_title")
    if experiment_title:
        final_fields["experiment_title"] = experiment_title
    experiment_summary = _join_values(
        [
            merged_prompt_metadata.get("experiment_short_description"),
            merged_prompt_metadata.get("score_set_short_description"),
        ]
    )
    if experiment_summary:
        final_fields["experiment_summary"] = experiment_summary

    license_label = merged_prompt_metadata.get("license")
    if license_label:
        final_fields["license_label"] = license_label
    if any(taxonomy == "Homo sapiens" for taxonomy in taxonomies):
        final_fields["species"] = "Homo sapiens"
    if merged_prompt_metadata:
        final_fields["data_modality"] = "MAVE"
    return final_fields


def merge_final_metadata(
    normalized_payload: dict[str, Any],
    mavedb_metadata_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
) -> dict[str, object]:
    """Merge one Step 2/4 payload with MaveDB fields and project to ``ObsSchema``."""
    mavedb_metadata_dir = Path(mavedb_metadata_dir).resolve()
    source_entries = _load_source_entries(normalized_payload, mavedb_metadata_dir)
    mavedb_fields = _extract_mavedb_fields(source_entries)

    final_payload = {
        field_name: normalized_payload.get(field_name)
        for field_name in OBS_SCHEMA_FIELDS
    }
    for field_name, field_value in mavedb_fields.items():
        if field_name in final_payload and _is_missing(final_payload[field_name]):
            final_payload[field_name] = field_value

    # Preserve provenance for auditability; these are operational fields rather
    # than AnnData ``obs`` columns and are intentionally kept outside the schema
    # projection above.
    for field_name in PROVENANCE_FIELDS:
        if field_name in normalized_payload:
            final_payload[field_name] = normalized_payload[field_name]
    return final_payload


def _resolve_input_json_dir(input_dir: str | Path) -> Path:
    input_dir = Path(input_dir).resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Normalized metadata directory not found: {input_dir}")
    if any(
        path.is_file()
        and path.suffix == ".json"
        and path.name not in EXCLUDED_INPUT_JSON_NAMES
        for path in input_dir.iterdir()
    ):
        return input_dir
    for subdirectory_name in ("step4_backfilled", "step2_normalized"):
        candidate = input_dir / subdirectory_name
        if candidate.is_dir():
            return _resolve_input_json_dir(candidate)
    raise FileNotFoundError(f"No normalized metadata JSON files found in: {input_dir}")


def finalize_metadata(
    normalized_metadata_dir: str | Path,
    output_dir: str | Path,
    log_file: str | Path,
    mavedb_metadata_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
    overwrite: bool = False,
    create_csv: bool = True,
) -> list[Path]:
    """Create one schema-projected Step 5 JSON object per normalized study."""
    input_json_dir = _resolve_input_json_dir(normalized_metadata_dir)
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    output_json_dir = get_final_json_dir(output_dir)
    output_json_dir.mkdir(parents=True, exist_ok=True)

    output_paths: list[Path] = []
    for input_path in sorted(input_json_dir.glob("*.json")):
        if input_path.name in EXCLUDED_INPUT_JSON_NAMES:
            continue
        output_path = get_final_output_path(input_path, output_dir)
        if output_path.is_file() and not overwrite:
            output_paths.append(output_path)
            continue
        try:
            normalized_payload = json.loads(input_path.read_text(encoding="utf-8"))
            if not isinstance(normalized_payload, dict):
                raise ValueError("normalized metadata payload must be a JSON object")
            final_payload = merge_final_metadata(
                normalized_payload=normalized_payload,
                mavedb_metadata_dir=mavedb_metadata_dir,
            )
            output_path.write_text(
                json.dumps(final_payload, indent=JSON_INDENT, ensure_ascii=True),
                encoding="utf-8",
            )
            output_paths.append(output_path)
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            append_log_line(log_file, f"Step 5 Error: {input_path}: {exc}")

    csv_path = get_final_csv_path(output_dir)
    should_write_csv = (
        create_csv and output_paths and (overwrite or not csv_path.is_file())
    )
    if should_write_csv:
        try:
            create_csv_from_curated_metadata_json(
                input_dir=output_json_dir,
                output_csv_path=csv_path,
                log_file=log_file,
            )
        except (OSError, ValueError) as exc:
            append_log_line(log_file, f"Step 5 CSV creation failed: {exc}")

    print_status_block(
        log_file,
        "Step 5 final metadata complete",
        f"Normalized metadata directory: {input_json_dir}",
        f"Output directory: {output_json_dir}",
        f"Studies finalized: {len(output_paths)}",
    )
    return output_paths


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for Step 5 final metadata assembly."""
    parser = argparse.ArgumentParser(
        description="Combine normalized LLM metadata with MaveDB metadata into final ObsSchema-shaped JSON objects.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--normalized-metadata-dir",
        type=Path,
        required=True,
        help="Directory containing Step 4 backfilled or Step 2 normalized JSON files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory used to write Step 5 final outputs.",
    )
    parser.add_argument(
        "--mavedb-metadata-dir",
        type=Path,
        default=MAVEDB_METADATA_OUTPUT_DIR,
        help="Directory containing raw MaveDB metadata JSON entries.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        required=True,
        help="Log file used for progress and errors.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing Step 5 output files.",
    )
    parser.add_argument(
        "--no-csv",
        action="store_false",
        dest="create_csv",
        help="Disable creation of a final metadata CSV.",
    )
    return parser


def main() -> None:
    """Run Step 5 final metadata assembly from the command line."""
    args = build_parser().parse_args()
    finalize_metadata(
        normalized_metadata_dir=args.normalized_metadata_dir,
        output_dir=args.output_dir,
        log_file=args.log_file,
        mavedb_metadata_dir=args.mavedb_metadata_dir,
        overwrite=args.overwrite,
        create_csv=args.create_csv,
    )


if __name__ == "__main__":
    main()
