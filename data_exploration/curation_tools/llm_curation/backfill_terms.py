"""Module and CLI runner for Step 4: Backfilling approved terms into Step 2 normalized outputs."""

import argparse
import json
import shutil
import traceback
from pathlib import Path
from typing import Any

from curation_tools.llm_curation.logging_utils import (
    append_log_line,
    print_status_block,
)
from curation_tools.llm_curation.metadata_extraction import (
    create_csv_from_curated_metadata_json,
)


def preview_backfill_changes(
    step2_dir: str | Path,
    decisions_data: dict[str, Any] | str | Path,
) -> list[dict[str, Any]]:
    """Preview field replacement records ('Other' -> approved_term) without modifying files on disk.

    Returns a list of preview audit record dicts:
    [
        {
            "source_file": "...",
            "dataset_id": "...",
            "field_name": "...",
            "previous_value": "Other",
            "new_value": "approved_term"
        },
        ...
    ]
    """
    step2_dir = Path(step2_dir).resolve()
    if not step2_dir.is_dir():
        return []

    if not list(step2_dir.glob("*.json")) and (step2_dir / "step2_normalized").is_dir():
        step2_dir = step2_dir / "step2_normalized"

    if isinstance(decisions_data, (str, Path)):
        decisions_file = Path(decisions_data).resolve()
        if not decisions_file.is_file():
            return []
        decisions = json.loads(decisions_file.read_text(encoding="utf-8"))
    elif isinstance(decisions_data, dict):
        decisions = decisions_data
    else:
        return []

    file_field_updates: dict[str, dict[str, str]] = {}

    for field_name, candidate_list in decisions.items():
        if not isinstance(candidate_list, list):
            continue
        for candidate in candidate_list:
            if candidate.get("status") == "Approved":
                approved_term = candidate.get("term")
                if not approved_term:
                    continue
                evidence_list = candidate.get("supporting_evidence", [])
                for ev in evidence_list:
                    source_file = None
                    if isinstance(ev, dict):
                        source_file = ev.get("source_file")
                    elif isinstance(ev, str):
                        if ev.endswith(".json"):
                            source_file = ev

                    if source_file:
                        source_filename = Path(source_file).name
                        if source_filename not in file_field_updates:
                            file_field_updates[source_filename] = {}
                        file_field_updates[source_filename][field_name] = approved_term

    preview_records: list[dict[str, Any]] = []

    for filename, field_updates in file_field_updates.items():
        target_path = step2_dir / filename
        if not target_path.exists():
            continue

        try:
            data = json.loads(target_path.read_text(encoding="utf-8"))
            dataset_id = data.get("dataset_id") or "unknown"

            for field_name, approved_term in field_updates.items():
                if field_name in data and data[field_name] == "Other":
                    preview_records.append(
                        {
                            "source_file": filename,
                            "dataset_id": dataset_id,
                            "field_name": field_name,
                            "previous_value": "Other",
                            "new_value": approved_term,
                        }
                    )
        except Exception:
            pass

    return preview_records


def backfill_approved_terms(
    step2_dir: str | Path,
    output_dir: str | Path,
    decisions_file: str | Path,
    log_file: str | Path,
    create_csv: bool = True,
) -> tuple[int, int, list[dict[str, Any]]]:
    """Copy Step 2 normalized JSON files to output_dir and replace 'Other' values with approved terms.

    Returns tuple: (total_files_copied, total_fields_updated, audit_records)
    """
    step2_dir = Path(step2_dir).resolve()
    output_dir = Path(output_dir).resolve()
    decisions_file = Path(decisions_file).resolve()
    log_file = Path(log_file).resolve()

    if not step2_dir.is_dir():
        raise FileNotFoundError(f"Step 2 directory not found: {step2_dir}")
    if not list(step2_dir.glob("*.json")) and (step2_dir / "step2_normalized").is_dir():
        step2_dir = step2_dir / "step2_normalized"
    if not decisions_file.is_file():
        raise FileNotFoundError(f"Decisions file not found: {decisions_file}")

    output_dir.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    # 1. Copy all Step 2 files into output_dir
    step2_files = sorted(step2_dir.glob("*.json"))
    for s2_file in step2_files:
        dest_file = output_dir / s2_file.name
        shutil.copy2(s2_file, dest_file)

    append_log_line(
        log_file,
        f"Step 4: Copied {len(step2_files)} files from {step2_dir} to {output_dir}",
    )

    # 2. Load decisions file
    decisions_data = json.loads(decisions_file.read_text(encoding="utf-8"))

    # 3. Process approved decisions and map target files -> field -> approved_term
    file_field_updates: dict[str, dict[str, str]] = {}

    for field_name, candidate_list in decisions_data.items():
        if not isinstance(candidate_list, list):
            continue
        for candidate in candidate_list:
            if candidate.get("status") == "Approved":
                approved_term = candidate.get("term")
                if not approved_term:
                    continue
                evidence_list = candidate.get("supporting_evidence", [])
                for ev in evidence_list:
                    source_file = None
                    if isinstance(ev, dict):
                        source_file = ev.get("source_file")
                    elif isinstance(ev, str):
                        if ev.endswith(".json"):
                            source_file = ev

                    if source_file:
                        source_filename = Path(source_file).name
                        if source_filename not in file_field_updates:
                            file_field_updates[source_filename] = {}
                        file_field_updates[source_filename][field_name] = approved_term

    # 4. Perform replacements on copied Step 4 files
    updated_fields_count = 0
    updated_files_count = 0
    audit_records: list[dict[str, Any]] = []

    for filename, field_updates in file_field_updates.items():
        target_path = output_dir / filename
        if not target_path.exists():
            append_log_line(
                log_file,
                f"Step 4 Warning: Target file {target_path} referenced in decisions does not exist in Step 2.",
            )
            continue

        try:
            data = json.loads(target_path.read_text(encoding="utf-8"))
            dataset_id = data.get("dataset_id") or "unknown"
            file_modified = False

            for field_name, approved_term in field_updates.items():
                if field_name in data and data[field_name] == "Other":
                    prev_val = data[field_name]
                    data[field_name] = approved_term
                    file_modified = True
                    updated_fields_count += 1
                    audit_records.append(
                        {
                            "source_file": filename,
                            "dataset_id": dataset_id,
                            "field_name": field_name,
                            "previous_value": prev_val,
                            "new_value": approved_term,
                        }
                    )
                    append_log_line(
                        log_file,
                        f"Step 4 Backfill: {filename}.{field_name}: '{prev_val}' -> '{approved_term}'",
                    )

            if file_modified:
                target_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
                updated_files_count += 1

        except Exception as exc:
            append_log_line(
                log_file,
                f"Step 4 Error: Failed updating {target_path}: {exc}\n{traceback.format_exc()}",
            )

    # Save audit records JSON
    audit_file = output_dir / "step4_backfill_audit.json"
    audit_file.write_text(json.dumps(audit_records, indent=2), encoding="utf-8")

    print_status_block(
        log_file,
        "Step 4 Backfill Complete",
        f"Step 2 Source Dir: {step2_dir}",
        f"Step 4 Output Dir: {output_dir}",
        f"Decisions File: {decisions_file}",
        f"Total Files Copied: {len(step2_files)}",
        f"Files Backfilled: {updated_files_count}",
        f"Fields Replaced ('Other' -> Approved): {updated_fields_count}",
        f"Audit Report Saved: {audit_file}",
    )

    # 5. Optionally create merged CSV
    if create_csv and len(step2_files) > 0:
        csv_path = output_dir / "step4_backfilled_metadata.csv"
        try:
            create_csv_from_curated_metadata_json(
                input_dir=output_dir,
                output_csv_path=csv_path,
                log_file=log_file,
            )
        except Exception as csv_exc:
            append_log_line(log_file, f"Step 4 CSV creation failed: {csv_exc}")

    return len(step2_files), updated_fields_count, audit_records


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for Step 4 backfilling."""
    parser = argparse.ArgumentParser(
        description="Copy Step 2 normalized JSON files into a Step 4 directory and backfill approved terms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--step2-dir",
        type=Path,
        required=True,
        help="Directory containing Step 2 normalized JSON files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory used to write Step 4 backfilled outputs.",
    )
    parser.add_argument(
        "--decisions-file",
        type=Path,
        required=True,
        help="Path to approved_ontology_terms.json containing approved decisions.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        required=True,
        help="Log file used for progress and errors.",
    )
    parser.add_argument(
        "--no-csv",
        action="store_false",
        dest="create_csv",
        help="Disable creation of a single CSV from the backfilled metadata JSON outputs.",
    )
    return parser


def main() -> None:
    """Main CLI entry point for Step 4 backfill."""
    args = build_parser().parse_args()
    backfill_approved_terms(
        step2_dir=args.step2_dir,
        output_dir=args.output_dir,
        decisions_file=args.decisions_file,
        log_file=args.log_file,
        create_csv=args.create_csv,
    )


if __name__ == "__main__":
    main()
