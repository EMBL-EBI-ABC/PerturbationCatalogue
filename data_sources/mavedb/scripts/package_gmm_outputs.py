#!/usr/bin/env python3
"""Package verified GMM batch outputs and create merged TSV tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import tempfile
from pathlib import Path
from typing import Any

NATIVE_TABLE_FILENAMES = {
    "scores": "score.tsv",
    "densities": "density.tsv",
    "params": "fitted_params.tsv",
    "posterior": "posterior.tsv",
    "fit_summary": "fit_summary.tsv",
    "loo": "loo.tsv",
}
QC_TSV_OUTPUTS = {
    "study_diagnostics",
    "eligibility",
    "exclusions",
    "smoke_selected_rows",
}
REQUIRED_OUTPUTS = {
    *NATIVE_TABLE_FILENAMES,
    "study_diagnostics",
    "run_metadata",
    "eligibility",
    "exclusions",
}


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


def model_output_root(plan: dict[str, Any]) -> Path:
    configuration = plan.get("configuration", {})
    profile = configuration.get("profile")
    model = configuration.get("model")
    if profile not in {"production", "smoke"} or model not in {
        "modes",
        "flexible",
    }:
        raise ValueError("GMM plan has an unsupported profile or model")
    if profile == "smoke":
        return Path("smoke") / "maves" / model
    return Path("maves") / model


def collection_paths(
    model_root: Path,
    name: str,
    relative: Path,
    batch_index: int,
) -> tuple[Path, Path | None]:
    batch_root = model_root / "batches" / f"batch{batch_index}"
    if name in NATIVE_TABLE_FILENAMES:
        filename = NATIVE_TABLE_FILENAMES[name]
        source_filename = f"{filename.removesuffix('.tsv')}{batch_index}.tsv"
        if relative != model_root / source_filename:
            raise ValueError(f"Unexpected path for GMM output {name!r}: {relative}")
        return batch_root / filename, model_root / "merged" / filename

    qc_root = model_root / "qc" / f"batch{batch_index}"
    if relative.parent != qc_root:
        raise ValueError(f"Unexpected path for GMM output {name!r}: {relative}")
    if name == "run_metadata" and relative.name == "run_metadata.json":
        return batch_root / "qc" / relative.name, None
    if name in QC_TSV_OUTPUTS and relative.suffix == ".tsv":
        return (
            batch_root / "qc" / relative.name,
            model_root / "merged" / "qc" / relative.name,
        )
    raise ValueError(f"Unsupported GMM output {name!r}: {relative}")


def merge_tsv_files(source_paths: list[Path], destination: Path) -> int:
    if not source_paths:
        raise ValueError(f"No TSV inputs supplied for {destination}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    expected_header = None
    rows_written = 0
    with destination.open("w", encoding="utf-8", newline="") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        for source in source_paths:
            with source.open(encoding="utf-8", newline="") as input_handle:
                reader = csv.reader(input_handle, delimiter="\t")
                header = next(reader, None)
                if not header:
                    raise ValueError(f"GMM TSV is empty: {source}")
                if expected_header is None:
                    expected_header = header
                    writer.writerow(header)
                elif header != expected_header:
                    raise ValueError(
                        f"GMM TSV headers differ while merging {destination}: "
                        f"{source}"
                    )

                for row_number, row in enumerate(reader, start=2):
                    if len(row) != len(expected_header):
                        raise ValueError(
                            f"GMM TSV row {row_number} in {source} has "
                            f"{len(row)} fields; expected {len(expected_header)}"
                        )
                    writer.writerow(row)
                    rows_written += 1
    return rows_written


def package_outputs(bundle: Path, output: Path) -> None:
    plan_path = bundle / "batch_plan.json"
    completion_path = bundle / "completion.json"
    if not plan_path.is_file() or not completion_path.is_file():
        raise ValueError("Collected GMM bundle is missing its plan or completion")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    batch_count = plan.get("batch_count")
    if (
        not isinstance(batch_count, int)
        or isinstance(batch_count, bool)
        or batch_count < 1
    ):
        raise ValueError("GMM plan has an invalid batch count")
    plan_identity = file_identity(plan_path)
    if completion.get("plan") != plan_identity:
        raise ValueError("GMM completion record belongs to a different plan")
    if completion.get("batch_count") != batch_count:
        raise ValueError("GMM completion batch count differs from its plan")
    model_root = model_output_root(plan)

    receipt_paths = {
        int(path.stem.removeprefix("batch")): path
        for path in (bundle / "batches").glob("batch*.json")
        if path.stem.removeprefix("batch").isdigit()
    }
    if set(receipt_paths) != set(range(1, batch_count + 1)):
        raise ValueError(
            "Collected GMM bundle has missing or unexpected batch receipts"
        )

    if output.exists():
        raise FileExistsError(f"Packaged GMM output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".gmm-package-", dir=output.parent
    ) as temporary:
        staged = Path(temporary) / "gmm"
        staged.mkdir()
        merge_sources: dict[Path, list[tuple[int, Path]]] = {}

        for index in sorted(receipt_paths):
            record = json.loads(receipt_paths[index].read_text(encoding="utf-8"))
            if (
                record.get("schema_version") != 1
                or record.get("batch_index") != index
                or record.get("plan") != plan_identity
            ):
                raise ValueError(f"GMM batch receipt {index} does not match its plan")
            outputs = record.get("outputs")
            if not isinstance(outputs, dict) or not REQUIRED_OUTPUTS.issubset(outputs):
                raise ValueError(f"GMM batch {index} is missing required outputs")

            for name, identity in outputs.items():
                if not isinstance(identity, dict):
                    raise ValueError(f"Invalid identity for GMM output {name!r}")
                relative = Path(identity.get("path", ""))
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError("GMM batch output escapes its bundle")
                source = bundle / relative
                if file_identity(source) != {
                    key: identity.get(key) for key in ("bytes", "sha256")
                }:
                    raise ValueError(f"GMM output identity mismatch: {relative}")

                destination_relative, merged_relative = collection_paths(
                    model_root, name, relative, index
                )
                destination = staged / destination_relative
                if destination.exists():
                    raise ValueError(f"GMM output collision: {destination_relative}")
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
                identity["published_path"] = destination_relative.as_posix()
                if merged_relative is not None:
                    merge_sources.setdefault(merged_relative, []).append(
                        (index, destination_relative)
                    )

            write_json(staged / "batches" / f"batch{index}.json", record)

        merged_outputs = {}
        for merged_relative, sources in sorted(merge_sources.items()):
            ordered_sources = [
                staged / relative
                for _, relative in sorted(sources, key=lambda item: item[0])
            ]
            destination = staged / merged_relative
            row_count = merge_tsv_files(ordered_sources, destination)
            merged_outputs[merged_relative.as_posix()] = {
                "path": merged_relative.as_posix(),
                "rows": row_count,
                **file_identity(destination),
            }

        shutil.copyfile(plan_path, staged / "batch_plan.json")
        completion["merged_outputs"] = merged_outputs
        write_json(staged / "completion.json", completion)
        staged.rename(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    package_outputs(args.bundle, args.out)


if __name__ == "__main__":
    main()
