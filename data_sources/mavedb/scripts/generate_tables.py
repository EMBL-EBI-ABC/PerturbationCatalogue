#!/usr/bin/env python3
"""Generate derived MAVE tables without replacing canonical processing."""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import csv
from datetime import datetime, timezone
import hashlib
import json
import logging
import math
from pathlib import Path
import sys
import tempfile
from typing import Any, Mapping

import pandas as pd

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.scripts import process_mavedb as processing

LOGGER = logging.getLogger(__name__)
TABLE_COLUMNS = {
    "singles": processing.SINGLE_COLUMNS,
    "doubles": processing.DOUBLE_COLUMNS,
    "wide": processing.WIDE_COLUMNS,
}
SCORE_FIELDS = ("hgvs_nt", "hgvs_pro", "hgvs_splice", "score", "score_raw")


def _check_identity(path: Path, expected: Mapping[str, Any], label: str) -> None:
    if processing.file_identity(path) != {
        key: expected[key] for key in ("bytes", "sha256")
    }:
        raise ValueError(f"{label} differs from the canonical processing manifest")


def _row_bytes(row: Mapping[str, Any]) -> bytes:
    values: list[Any] = [int(row["source_row"])]
    for field in SCORE_FIELDS:
        value = row[field]
        if pd.isna(value) or value in ("", "NA"):
            value = None
        elif field == "score":
            numeric = float(value)
            value = None if math.isnan(numeric) else numeric.hex()
        else:
            value = str(value)
        values.append(value)
    return (json.dumps(values, ensure_ascii=False) + "\n").encode("utf-8")


def _canonical_scores(path: Path, selected: set[str]) -> dict[str, str]:
    """Hash analysis-relevant source rows in bounded memory, by dataset."""
    digests = {dataset_id: hashlib.sha256() for dataset_id in selected}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"dataset_id", "source_row", *SCORE_FIELDS}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Canonical scores are missing source analysis columns")
        for row in reader:
            if row["dataset_id"] not in selected:
                raise ValueError(
                    "Canonical scores contain a dataset absent from metadata"
                )
            digests[row["dataset_id"]].update(_row_bytes(row))
    return {dataset_id: digest.hexdigest() for dataset_id, digest in digests.items()}


def generate_tables(
    input_path: Path,
    processed: Path,
    taxonomy: Path,
    output: Path,
) -> dict[str, Any]:
    """Write source-compatible tables for exactly the completed canonical selection.

    A processing manifest is required to establish the target metadata/taxonomy
    snapshot. Source score rows are additionally compared with the canonical
    rows, because older ingest manifests do not hash extracted CSV files.
    """
    input_path, processed, taxonomy, output = map(
        Path, (input_path, processed, taxonomy, output)
    )
    if output.exists():
        raise FileExistsError(f"Derived-table output already exists: {output}")
    manifest_path = processed / "processing_manifest.json"
    if not manifest_path.is_file():
        raise ValueError("Derived tables require processing_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported canonical processing manifest")
    metadata_path = processing.resolve_metadata_path(input_path)
    _check_identity(metadata_path, manifest["ingest"]["metadata"], "Source metadata")
    _check_identity(taxonomy, manifest["taxonomy"], "Taxonomy")
    for name in ("scores.tsv", "metadata.tsv", "dataset_processing.tsv"):
        _check_identity(processed / name, manifest["outputs"][name], name)
    processing.validate_tables(processed)
    with (processed / "metadata.tsv").open(encoding="utf-8", newline="") as handle:
        selected = {row["dataset_id"] for row in csv.DictReader(handle, delimiter="\t")}
    expected_scores = _canonical_scores(processed / "scores.tsv", selected)
    datasets = processing.load_datasets(metadata_path)
    processing.add_kingdoms(
        {key: dataset["targets"] for key, dataset in datasets.items()}, taxonomy
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".mavedb-tables-", dir=output.parent
    ) as temporary:
        staged = Path(temporary) / "tables"
        (staged / "maves").mkdir(parents=True)
        seen = set()
        source_files = {}
        with ExitStack() as stack:
            tables = {}
            for name, columns in TABLE_COLUMNS.items():
                handle = stack.enter_context(
                    (staged / "maves" / f"{name}.tsv").open(
                        "w", encoding="utf-8", newline=""
                    )
                )
                processing.write_tsv(pd.DataFrame(columns=columns), handle)
                tables[name] = handle
            for score_file in sorted(
                (metadata_path.parent / "csv").glob("*.scores.csv")
            ):
                dataset_id = score_file.name.removesuffix(".scores.csv").replace(
                    "urn-mavedb-", "urn:mavedb:", 1
                )
                if dataset_id not in selected:
                    continue
                if dataset_id in seen:
                    raise ValueError(f"Duplicate score file for {dataset_id}")
                seen.add(dataset_id)
                dataset = datasets.get(dataset_id)
                if dataset is None or not dataset["targets"]:
                    raise ValueError(f"Missing source targets for {dataset_id}")
                scores = processing.read_scores(score_file)
                digest = hashlib.sha256()
                fields = ("source_row", *SCORE_FIELDS)
                for values in scores.loc[:, fields].itertuples(index=False, name=None):
                    digest.update(_row_bytes(dict(zip(fields, values))))
                if digest.hexdigest() != expected_scores[dataset_id]:
                    raise ValueError(
                        f"Source scores differ from canonical rows: {dataset_id}"
                    )
                categories = processing.category_scores(scores, dataset["targets"])
                if not categories:
                    raise ValueError(
                        f"Canonical dataset has no source analysis: {dataset_id}"
                    )
                scores["dataset_id"] = dataset_id
                scores["score_urn"] = dataset_id
                processing.dataset_analysis(dataset, scores, categories, tables)
                source_files[score_file.name] = processing.file_identity(score_file)
        if seen != selected:
            raise ValueError(
                f"Missing source score files for: {', '.join(sorted(selected - seen))}"
            )
        outputs = {}
        for name in TABLE_COLUMNS:
            relative = f"maves/{name}.tsv"
            path = staged / relative
            with path.open(encoding="utf-8", newline="") as handle:
                rows = sum(1 for _ in csv.reader(handle, delimiter="\t")) - 1
            outputs[relative] = {
                "path": relative,
                "rows": rows,
                **processing.file_identity(path),
            }
        completion = {
            "schema_version": 1,
            "status": "completed",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "processing_manifest": processing.file_identity(manifest_path),
            "source_metadata": processing.file_identity(metadata_path),
            "source_scores": source_files,
            "taxonomy": processing.file_identity(taxonomy),
            "implementation": {
                "generator": processing.file_identity(Path(__file__)),
                "processor": processing.file_identity(Path(processing.__file__)),
            },
            "dataset_ids": sorted(selected),
            "outputs": outputs,
        }
        (staged / "completion.json").write_text(
            json.dumps(completion, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        staged.rename(output)
    return completion


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--processed", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        generate_tables(args.input, args.processed, args.taxonomy, args.output)
    except (OSError, ValueError, KeyError, pd.errors.ParserError) as error:
        parser.error(str(error))
    LOGGER.info("Wrote derived tables to %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
