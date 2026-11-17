"""Validate existing published outputs without invoking their producers."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import tempfile
from typing import Any

from data_sources.mavedb.analysis._gmm_pipeline import (
    collect_batches,
    file_identity,
    read_plan,
    write_json,
)
from data_sources.mavedb.scripts.package_gmm_outputs import (
    collection_paths,
    merge_tsv_files,
    model_output_root,
    REQUIRED_OUTPUTS,
)


def _record(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _path(root: Path, relative: str) -> Path:
    path = Path(relative)
    if not relative or path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Output path escapes its bundle: {relative!r}")
    source = root / path
    if not source.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Output path escapes its bundle: {relative!r}")
    return source


def _verify(path: Path, identity: dict[str, Any]) -> None:
    if not isinstance(identity, dict) or file_identity(path) != {
        key: identity.get(key) for key in ("bytes", "sha256")
    }:
        raise ValueError(f"Output identity mismatch: {path}")


def unpack_gmm(bundle: Path, output: Path) -> dict[str, Any]:
    """Check a packaged result and restore a portable collected GMM bundle.

    Receipts are checked against their saved plan, never against today's model
    implementation. No model imports or sampler executions are needed. The
    reconstructed bundle can be passed directly to ``package_outputs``.
    """
    bundle, output = Path(bundle), Path(output)
    if output.exists():
        raise FileExistsError(f"Adoption output already exists: {output}")
    plan_path = bundle / "batch_plan.json"
    plan = read_plan(plan_path)
    model_root = model_output_root(plan)
    completion = _record(bundle / "completion.json")
    if completion.get("schema_version") != 1:
        raise ValueError("Unsupported GMM completion record")
    if (
        completion.get("plan") != file_identity(plan_path)
        or completion.get("batch_count") != plan["batch_count"]
        or plan["batch_count"] < 1
    ):
        raise ValueError("GMM completion differs from its saved plan")
    receipt_paths = sorted((bundle / "batches").glob("*.json"))
    expected_receipts = {
        f"batch{index}.json" for index in range(1, plan["batch_count"] + 1)
    }
    if {path.name for path in receipt_paths} != expected_receipts:
        raise ValueError("GMM bundle has missing or unexpected batch receipts")
    merged = completion.get("merged_outputs")
    if not isinstance(merged, dict):
        raise ValueError("GMM bundle is missing merged output identities")

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".gmm-adopt-", dir=output.parent
    ) as temporary:
        temporary_root = Path(temporary)
        batches = []
        merge_sources: dict[str, list[tuple[int, Path]]] = {}
        for receipt_path in receipt_paths:
            record = _record(receipt_path)
            index = record.get("batch_index")
            if receipt_path.name != f"batch{index}.json":
                raise ValueError("GMM receipt filename differs from its batch index")
            outputs = record.get("outputs")
            if not isinstance(outputs, dict) or not REQUIRED_OUTPUTS.issubset(outputs):
                raise ValueError("GMM batch receipt is missing required outputs")
            batch = temporary_root / f"batch-{index}"
            batch.mkdir()
            for name, identity in outputs.items():
                if not isinstance(identity, dict):
                    raise ValueError(f"Invalid GMM output identity: {name}")
                native = _path(batch, identity.get("path", ""))
                expected_published, merged_relative = collection_paths(
                    model_root, name, Path(identity.get("path", "")), index
                )
                if identity.get("published_path") != expected_published.as_posix():
                    raise ValueError(f"Unexpected published GMM output path: {name}")
                source = _path(bundle, identity["published_path"])
                _verify(source, identity)
                native.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, native)
                if merged_relative is not None:
                    merge_sources.setdefault(merged_relative.as_posix(), []).append(
                        (index, native)
                    )
                identity.pop("published_path")
            write_json(batch / "summary.json", record)
            batches.append(batch)

        if set(merged) != set(merge_sources):
            raise ValueError("GMM bundle has missing or unexpected merged outputs")
        for relative, sources in merge_sources.items():
            identity = merged[relative]
            if not isinstance(identity, dict) or identity.get("path") != relative:
                raise ValueError(f"Invalid merged output identity: {relative}")
            _verify(_path(bundle, relative), identity)
            # Also prove merged rows/headers agree with the retained batch data.
            reconstructed = temporary_root / "merged" / relative
            rows = merge_tsv_files(
                [path for _, path in sorted(sources)],
                reconstructed,
                drop_columns=(
                    {"dataset_id"}
                    if relative == (model_root / "merged/score.tsv").as_posix()
                    else None
                ),
            )
            if rows != identity.get("rows") or file_identity(reconstructed) != {
                key: identity.get(key) for key in ("bytes", "sha256")
            }:
                raise ValueError(f"Merged output differs from GMM batches: {relative}")

        collected = temporary_root / "collected"
        collect_batches(plan_path, batches, collected)
        regenerated = _record(collected / "completion.json")
        if regenerated != {
            key: value for key, value in completion.items() if key != "merged_outputs"
        }:
            raise ValueError("GMM completion differs from verified batch receipts")
        collected.rename(output)
    return plan


def adopt_processing(path: Path) -> dict[str, Any]:
    """Validate canonical tables and their original processing manifest."""
    from data_sources.mavedb.scripts.process_mavedb import validate_tables

    path = Path(path)
    manifest = _record(path / "processing_manifest.json")
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported processing manifest")
    outputs = manifest.get("outputs")
    required = {"scores.tsv", "metadata.tsv", "dataset_processing.tsv"}
    if manifest.get("settings", {}).get("compatibility_tables"):
        required.update(f"maves/{name}.tsv" for name in ("singles", "doubles", "wide"))
    if not isinstance(outputs, dict) or not required.issubset(outputs):
        raise ValueError("Processing manifest is missing required outputs")
    for name, identity in outputs.items():
        if not isinstance(identity, dict) or identity.get("path") != name:
            raise ValueError(f"Invalid processing output identity: {name}")
        _verify(_path(path, name), identity)
    if validate_tables(path) != manifest.get("counts"):
        raise ValueError("Processing counts differ from its manifest")
    return manifest
