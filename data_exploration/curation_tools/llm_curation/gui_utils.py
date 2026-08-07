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


def _normalized_json_files(directory: Path) -> list[Path]:
    """Return normalized JSON artifacts from a directory or its nested output folder."""
    directory = Path(directory).resolve()
    files = [
        path
        for path in sorted(directory.glob("*.json"))
        if not path.name.endswith("_audit.json")
    ]
    if not files and (directory / "step2_normalized").is_dir():
        files = [
            path
            for path in sorted((directory / "step2_normalized").glob("*.json"))
            if not path.name.endswith("_audit.json")
        ]
        directory = directory / "step2_normalized"
    return files


def _pipeline_artifact_files(directory: Path) -> list[Path]:
    """Return files that participate in pipeline freshness signatures."""
    directory = Path(directory).resolve()
    if not directory.is_dir():
        return []

    files = [
        path
        for path in sorted(directory.iterdir())
        if path.is_file()
        and path.suffix in {".json", ".md"}
        and not path.name.endswith("_audit.json")
        and path.name != "pipeline_manifest.json"
    ]
    if not files and (directory / "step2_normalized").is_dir():
        return _pipeline_artifact_files(directory / "step2_normalized")
    return files


def calculate_file_signature(file_path: Path) -> str | None:
    """Return a SHA-256 signature for a file, or None when it is unavailable."""
    file_path = Path(file_path).resolve()
    if not file_path.is_file():
        return None
    return sha256(file_path.read_bytes()).hexdigest()


def get_publication_dois_for_source_file(
    source_file: str,
    normalized_dir: Path,
    mapping_file: Path,
) -> list[str]:
    """Resolve publication DOIs from a Step 2/4 source file's MaveDB URNs."""
    from curation_tools.llm_curation.mavedb.processing import format_urn_for_filename

    normalized_dir = Path(normalized_dir).resolve()
    source_path = normalized_dir / source_file
    if not source_path.is_file() and (normalized_dir / "step2_normalized").is_dir():
        source_path = normalized_dir / "step2_normalized" / source_file
    if not source_path.is_file():
        matches = list(normalized_dir.rglob(source_file))
        source_path = matches[0] if matches else source_path

    try:
        urn_to_dois = json.loads(
            Path(mapping_file).resolve().read_text(encoding="utf-8")
        )
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return []

    source_urns: list[str] = []
    if source_path.is_file():
        try:
            source_payload = json.loads(source_path.read_text(encoding="utf-8"))
            source_urns = source_payload.get("__source_urns", [])
        except (json.JSONDecodeError, OSError):
            source_urns = []

    dois: set[str] = set()
    for source_urn in source_urns:
        urn_variants = {
            str(source_urn),
            str(source_urn).removeprefix("urn:"),
        }
        for urn_variant in urn_variants:
            dois.update(urn_to_dois.get(urn_variant, []))

    if not dois:
        for urn, urn_dois in urn_to_dois.items():
            if f"{format_urn_for_filename(urn)}.json" == source_file:
                dois.update(urn_dois)

    return sorted(dois)


def calculate_directory_signature(directory: Path) -> str | None:
    """Return a stable SHA-256 signature for pipeline artifacts in a directory."""
    files = _pipeline_artifact_files(directory)
    if not files:
        return None

    digest = sha256()
    for file_path in files:
        digest.update(file_path.name.encode("utf-8"))
        digest.update(file_path.read_bytes())
    return digest.hexdigest()


def load_pipeline_manifest(manifest_path: Path) -> dict[str, Any]:
    """Load the pipeline manifest, returning an empty structure when absent or invalid."""
    manifest_path = Path(manifest_path).resolve()
    if not manifest_path.is_file():
        return {"steps": {}}
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {"steps": {}}
    return manifest if isinstance(manifest, dict) else {"steps": {}}


def record_pipeline_step(
    manifest_path: Path,
    step: str,
    input_dir: Path | None = None,
    output_dir: Path | None = None,
    selected_items: list[str] | None = None,
    decisions_file: Path | None = None,
    source_step: str | None = None,
) -> dict[str, Any]:
    """Record a successful pipeline step and its freshness metadata."""
    manifest_path = Path(manifest_path).resolve()
    manifest = load_pipeline_manifest(manifest_path)
    manifest.setdefault("steps", {})
    completed_at = datetime.now(timezone.utc).isoformat()
    step_record: dict[str, Any] = {
        "status": "completed",
        "run_id": f"{step}-{completed_at}",
        "completed_at": completed_at,
    }

    if input_dir is not None:
        input_dir = Path(input_dir).resolve()
        step_record["input_dir"] = str(input_dir)
        step_record["input_signature"] = calculate_directory_signature(input_dir)
    if output_dir is not None:
        output_dir = Path(output_dir).resolve()
        step_record["output_dir"] = str(output_dir)
        step_record["output_signature"] = calculate_directory_signature(output_dir)
    if selected_items is not None:
        step_record["selected_items"] = sorted(selected_items)
    if decisions_file is not None:
        step_record["decisions_file"] = str(Path(decisions_file).resolve())
        step_record["decisions_signature"] = calculate_file_signature(decisions_file)

    if source_step:
        source_record = manifest["steps"].get(source_step, {})
        step_record["source_step"] = source_step
        step_record["source_run_id"] = source_record.get("run_id")
        step_record["source_completed_at"] = source_record.get("completed_at")
        step_record["source_output_signature"] = source_record.get("output_signature")

    manifest["steps"][step] = step_record
    timestamp_key = {
        "step1": "timestamp_evidence_extraction",
        "step2": "timestamp_term_normalization",
        "step3a": "timestamp_candidate_discovery",
        "step4": "timestamp_backfill",
        "step5": "timestamp_final_metadata",
    }.get(step)
    if timestamp_key:
        manifest[timestamp_key] = completed_at

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return step_record


def resolve_effective_normalized_dir(
    step2_dir: Path,
    step4_dir: Path,
    manifest_path: Path | None = None,
    decisions_file: Path | None = None,
) -> Path:
    """Use a manifest-verified Step 4 output, otherwise use Step 2 artifacts."""
    step2_dir = Path(step2_dir).resolve()
    step4_dir = Path(step4_dir).resolve()
    step2_files = _normalized_json_files(step2_dir)
    step4_files = _normalized_json_files(step4_dir)

    if not step4_files:
        return step2_dir

    if manifest_path is not None:
        manifest = load_pipeline_manifest(manifest_path)
        step2_record = manifest.get("steps", {}).get("step2", {})
        step4_record = manifest.get("steps", {}).get("step4", {})
        current_step2_signature = calculate_directory_signature(step2_dir)
        current_step4_signature = calculate_directory_signature(step4_dir)
        current_decisions_signature = (
            calculate_file_signature(decisions_file) if decisions_file else None
        )
        source_run_matches = not step4_record.get("source_run_id") or (
            step4_record.get("source_run_id") == step2_record.get("run_id")
        )
        source_timestamp_matches = not step4_record.get("source_completed_at") or (
            step4_record.get("source_completed_at") == step2_record.get("completed_at")
        )

        step4_is_current = (
            step4_record.get("status") == "completed"
            and step4_record.get("input_dir") == str(step2_dir)
            and step4_record.get("output_dir") == str(step4_dir)
            and step4_record.get("source_step") == "step2"
            and source_run_matches
            and source_timestamp_matches
            and step4_record.get("input_signature") == current_step2_signature
            and step4_record.get("output_signature") == current_step4_signature
            and step4_record.get("decisions_signature") == current_decisions_signature
        )
        return step4_dir if step4_is_current else step2_dir

    if not step2_files:
        return step4_dir

    latest_step2_mtime = max(path.stat().st_mtime for path in step2_files)
    latest_step4_mtime = max(path.stat().st_mtime for path in step4_files)
    return step4_dir if latest_step4_mtime >= latest_step2_mtime else step2_dir


def is_candidate_discovery_current(
    candidates_file: Path,
    effective_input_dir: Path,
    manifest_path: Path,
) -> bool:
    """Return whether the Step 3a output matches the current input artifacts."""
    candidates_file = Path(candidates_file).resolve()
    effective_input_dir = Path(effective_input_dir).resolve()
    manifest = load_pipeline_manifest(manifest_path)
    step3_record = manifest.get("steps", {}).get("step3a", {})
    output_dir = Path(step3_record.get("output_dir", "")).resolve()

    return (
        step3_record.get("status") == "completed"
        and step3_record.get("input_dir") == str(effective_input_dir)
        and step3_record.get("input_signature")
        == calculate_directory_signature(effective_input_dir)
        and output_dir == candidates_file.parent
        and candidates_file.is_file()
        and step3_record.get("output_signature")
        == calculate_directory_signature(output_dir)
    )


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
