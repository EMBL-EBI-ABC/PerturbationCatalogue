"""Helper functions for Streamlit curation GUI status tables and file tracking."""

import json
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any


def get_step1_file_status(input_dir: Path, output_dir: Path) -> list[dict[str, Any]]:
    """Scan Step 1 input markdown files and check status against output JSON files."""
    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir).resolve()

    if not input_dir.is_dir():
        return []

    md_files = sorted(input_dir.glob("*.md"))
    records = []

    for md_file in md_files:
        stem = md_file.stem
        size_kb = md_file.stat().st_size / 1024.0

        # Check matching output JSON files in step1 output_dir
        matching_outputs = []
        if output_dir.is_dir():
            matching_outputs = [
                out_p.name
                for out_p in output_dir.glob("*.json")
                if out_p.name.startswith(stem)
            ]

        status = "Completed" if matching_outputs else "Pending"
        records.append(
            {
                "file_name": md_file.name,
                "size_kb": f"{size_kb:.1f} KB",
                "status": status,
                "output_count": len(matching_outputs),
                "output_files": (
                    ", ".join(matching_outputs) if matching_outputs else "None"
                ),
            }
        )

    return records


def get_mavedb_urn_status(
    mapping_file: Path,
    output_dir: Path,
    metadata_dir: Path,
    target_urns: set[str] | None = None,
) -> list[dict[str, Any]]:
    """Scan MaveDB URNs from mapping file and check status against step 1 output directory."""
    from curation_tools.llm_curation.mavedb.processing import (
        format_urn_for_filename,
        load_mavedb_urn_to_dois,
    )

    mapping_file = Path(mapping_file).resolve()
    output_dir = Path(output_dir).resolve()
    metadata_dir = Path(metadata_dir).resolve()

    if not mapping_file.is_file():
        return []

    try:
        urn_to_dois = load_mavedb_urn_to_dois(mapping_file)
    except Exception:
        return []

    records = []
    for urn, dois in sorted(urn_to_dois.items()):
        if target_urns and urn not in target_urns:
            continue

        urn_stem = format_urn_for_filename(urn)
        out_target = output_dir / f"{urn_stem}.json"
        if (
            not out_target.exists()
            and (output_dir / "step1_evidence" / f"{urn_stem}.json").exists()
        ):
            out_target = output_dir / "step1_evidence" / f"{urn_stem}.json"

        status = "Completed" if out_target.exists() else "Pending"

        title = "-"
        meta_file = metadata_dir / f"{urn_stem}.json"
        if meta_file.exists():
            try:
                data = json.loads(meta_file.read_text(encoding="utf-8"))
                title = (
                    data.get("title") or data.get("experiment", {}).get("title") or "-"
                )
            except Exception:
                pass

        records.append(
            {
                "urn": urn,
                "title": title,
                "primary_dois": ", ".join(dois) if dois else "None",
                "status": status,
                "output_json": out_target.name if out_target.exists() else "-",
            }
        )

    return records


def get_step2_file_status(step1_dir: Path, step2_dir: Path) -> list[dict[str, Any]]:
    """Scan Step 1 evidence files and check status & 'Other' count in Step 2 outputs."""
    step1_dir = Path(step1_dir).resolve()
    step2_dir = Path(step2_dir).resolve()

    if not step1_dir.is_dir():
        if (step1_dir / "step1_evidence").is_dir():
            step1_dir = step1_dir / "step1_evidence"
        else:
            return []

    s1_files = [
        f
        for f in sorted(step1_dir.glob("*.json"))
        if not f.name.endswith("_audit.json")
    ]
    if not s1_files and (step1_dir / "step1_evidence").is_dir():
        s1_files = [
            f
            for f in sorted((step1_dir / "step1_evidence").glob("*.json"))
            if not f.name.endswith("_audit.json")
        ]

    records = []

    for s1_file in s1_files:
        s2_target = step2_dir / s1_file.name
        if (
            not s2_target.exists()
            and (step2_dir / "step2_normalized" / s1_file.name).exists()
        ):
            s2_target = step2_dir / "step2_normalized" / s1_file.name

        status = "Pending"
        other_count = 0

        if s2_target.exists():
            try:
                data = json.loads(s2_target.read_text(encoding="utf-8"))
                status = "Completed"
                if isinstance(data, dict):
                    other_count = sum(1 for v in data.values() if v == "Other")
            except Exception:
                status = "Error"

        records.append(
            {
                "file_name": s1_file.name,
                "status": status,
                "other_fields_count": other_count if status == "Completed" else "-",
                "output_path": str(s2_target) if s2_target.exists() else "-",
            }
        )

    return records


def get_step3_other_corpus_summary(
    step1_dir: Path,
    step2_dir: Path,
    selected_files: list[str | Path] | None = None,
) -> list[dict[str, Any]]:
    """Aggregate 'Other' fields across Step 2 files for Step 3a pre-run analysis."""
    from curation_tools.llm_curation.candidate_discovery import (
        aggregate_unmapped_evidence,
    )

    step1_dir = Path(step1_dir).resolve()
    step2_dir = Path(step2_dir).resolve()

    if not step1_dir.is_dir() or not step2_dir.is_dir():
        return []

    try:
        unmapped = aggregate_unmapped_evidence(
            step1_dir=step1_dir,
            step2_dir=step2_dir,
            log_file=step2_dir.parent / "step3_summary.log",
            selected_files=selected_files,
        )
        records = []
        for field_name, evidence_list in unmapped.items():
            records.append(
                {
                    "field_name": field_name,
                    "other_instances_count": len(evidence_list),
                    "sample_evidence": (
                        (
                            evidence_list[0].get("evidence")
                            if isinstance(evidence_list[0], dict)
                            else str(evidence_list[0])
                        )
                        if evidence_list
                        else ""
                    ),
                }
            )
        records.sort(key=lambda x: x["other_instances_count"], reverse=True)
        return records
    except Exception:
        return []


def read_last_log_lines(log_file: Path, num_lines: int = 25) -> str:
    """Return the last num_lines lines from a log file."""
    log_file = Path(log_file)
    if not log_file.exists():
        return "Log file not created yet."
    try:
        lines = log_file.read_text(encoding="utf-8").splitlines()
        return "\n".join(lines[-num_lines:])
    except Exception as e:
        return f"Error reading log file: {e}"
