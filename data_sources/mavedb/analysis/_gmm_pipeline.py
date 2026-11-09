"""Portable plans and checked collection for independently executed GMM batches."""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any


def file_identity(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def configuration(config: Any) -> dict[str, Any]:
    return {
        name: getattr(config, name)
        for name in (
            "profile",
            "model",
            "input_format",
            "study_type",
            "variants_per_batch",
            "max_studies",
            "max_study_rows",
        )
    } | {"study_ids": list(config.study_ids), "settings": asdict(config.settings)}


def implementation_identity() -> dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    paths = [root / "scripts/fit_gmm_batch.py"]
    paths.extend(
        root / "analysis" / name
        for name in (
            "gmm.py",
            "_gmm_batch.py",
            "_gmm_inputs.py",
            "_gmm_stan.py",
            "_gmm_pipeline.py",
        )
    )
    paths.extend(sorted((root / "analysis/stan").glob("*.stan")))
    return {str(path.relative_to(root)): file_identity(path) for path in paths}


def make_plan(config: Any, summary: Any) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "status": "batch_plan",
        "configuration": configuration(config),
        "inputs": {
            "scores": file_identity(config.scores_path),
            "metadata": file_identity(config.meta_path),
        },
        "implementation": implementation_identity(),
        "batch_count": summary.batch_count,
        "assignments": [list(row) for row in summary.assignments],
    }


def read_plan(path: Path) -> dict[str, Any]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != 1 or plan.get("status") != "batch_plan":
        raise ValueError("Unsupported GMM batch plan")
    count = plan.get("batch_count")
    if not isinstance(count, int) or isinstance(count, bool) or count < 0:
        raise ValueError("GMM plan requires a nonnegative integer batch count")
    for row in plan["assignments"]:
        if (
            len(row) != 3
            or not isinstance(row[0], int)
            or isinstance(row[0], bool)
            or not isinstance(row[1], str)
            or not row[1]
            or not isinstance(row[2], int)
            or isinstance(row[2], bool)
            or row[2] < 0
        ):
            raise ValueError("Invalid GMM batch assignment")
    indices = {row[0] for row in plan["assignments"]}
    if indices != set(range(1, plan["batch_count"] + 1)):
        raise ValueError("GMM plan has inconsistent batch indices")
    studies = [row[1] for row in plan["assignments"]]
    if len(studies) != len(set(studies)):
        raise ValueError("GMM plan contains duplicate studies")
    return plan


def verify_plan(config: Any, path: Path) -> dict[str, Any]:
    plan = read_plan(path)
    if config.index < 1 or config.index > plan["batch_count"]:
        raise ValueError("A planned fit requires a positive, in-range --index")
    if configuration(config) != plan["configuration"]:
        raise ValueError("GMM settings differ from the frozen batch plan")
    if implementation_identity() != plan["implementation"]:
        raise ValueError("GMM implementation differs from the frozen batch plan")
    for name, source in (
        ("scores", config.scores_path),
        ("metadata", config.meta_path),
    ):
        if file_identity(source) != plan["inputs"][name]:
            raise ValueError(f"GMM {name} input differs from the frozen batch plan")
    return plan


def batch_record(config: Any, summary: Any, plan_path: Path) -> dict[str, Any]:
    plan = read_plan(plan_path)
    expected = [row[1] for row in plan["assignments"] if row[0] == config.index]
    if set(summary.selected_study_ids) != set(expected):
        raise ValueError("Fitted studies differ from the frozen batch assignment")
    if [list(row) for row in summary.assignments] != plan["assignments"]:
        raise ValueError("Recomputed GMM batch assignments differ from the plan")
    outputs = {}
    for name, path in summary.output_paths.items():
        relative = path.resolve().relative_to(Path(config.out).resolve()).as_posix()
        outputs[name] = {"path": relative, **file_identity(path)}
    return {
        "schema_version": 1,
        "status": "completed" if summary.failed_study_count == 0 else "study_failures",
        "plan": file_identity(plan_path),
        "batch_index": config.index,
        "study_ids": list(summary.selected_study_ids),
        "processed_study_count": summary.processed_study_count,
        "successful_study_count": summary.successful_study_count,
        "failed_study_count": summary.failed_study_count,
        "score_rows_written": summary.score_rows_written,
        "failures": [
            asdict(item) if is_dataclass(item) else str(item)
            for item in summary.failures
        ],
        "outputs": outputs,
    }


def collect_batches(plan_path: Path, batches: list[Path], output: Path) -> None:
    """Publish only a complete, checksum-verified set of successful native batches."""
    plan = read_plan(plan_path)
    if not plan["batch_count"]:
        raise ValueError("GMM plan contains no eligible studies")
    plan_identity = file_identity(plan_path)
    seen = set()
    records = []
    if output.exists():
        raise FileExistsError(f"GMM collection output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".gmm-collect-", dir=output.parent
    ) as temporary:
        staged = Path(temporary) / "gmm"
        staged.mkdir()
        for batch in batches:
            record = json.loads((batch / "summary.json").read_text(encoding="utf-8"))
            if record.get("schema_version") != 1:
                raise ValueError("Unsupported GMM batch result")
            index = record["batch_index"]
            expected = {row[1] for row in plan["assignments"] if row[0] == index}
            if not expected or index in seen:
                raise ValueError("Unexpected or duplicate GMM batch result")
            if record["plan"] != plan_identity:
                raise ValueError("GMM batch belongs to a different plan")
            if (
                record["status"] != "completed"
                or record["failed_study_count"] != 0
                or record["processed_study_count"] != len(expected)
                or record["successful_study_count"] != len(expected)
                or len(record["study_ids"]) != len(expected)
                or set(record["study_ids"]) != expected
            ):
                raise ValueError(
                    f"GMM batch {index} did not successfully fit its assigned studies"
                )
            required = {
                "scores",
                "densities",
                "params",
                "posterior",
                "fit_summary",
                "loo",
                "study_diagnostics",
                "run_metadata",
                "eligibility",
                "exclusions",
            }
            if not required.issubset(record["outputs"]):
                raise ValueError("GMM batch is missing required native outputs or QC")
            for identity in record["outputs"].values():
                relative = Path(identity["path"])
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError("GMM batch output escapes its bundle")
                source = batch / relative
                if file_identity(source) != {
                    key: identity[key] for key in ("bytes", "sha256")
                }:
                    raise ValueError(f"GMM output identity mismatch: {relative}")
                destination = staged / relative
                if destination.exists():
                    raise ValueError(f"GMM output collision: {relative}")
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
            write_json(staged / "batches" / f"batch{index}.json", record)
            seen.add(index)
            records.append(record)
        if seen != set(range(1, plan["batch_count"] + 1)):
            raise ValueError("Missing GMM batch results")
        shutil.copyfile(plan_path, staged / "batch_plan.json")
        write_json(
            staged / "completion.json",
            {
                "schema_version": 1,
                "status": "completed",
                "plan": plan_identity,
                "batch_count": len(seen),
                "successful_study_count": sum(
                    row["successful_study_count"] for row in records
                ),
                "failed_study_count": 0,
                "score_rows_written": sum(row["score_rows_written"] for row in records),
            },
        )
        staged.rename(output)
