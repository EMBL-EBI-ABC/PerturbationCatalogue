#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tempfile
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd

from io_schemas import DEA_SCHEMA, GSEA_SCHEMA


def normalize_leading_edge(value: Any) -> list[str]:
    if value is None:
        return []
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, str):
        if not value:
            return []
        delimiter = ";" if ";" in value else ","
        return [item.strip() for item in value.split(delimiter) if item.strip()]
    return [str(item) for item in value if str(item)]


def frame_to_table(
    frame: pd.DataFrame,
    schema: pa.Schema,
    dataset_id: str,
    ingested_at: pd.Timestamp,
) -> pa.Table:
    arrays = []
    for field in schema:
        if field.name == "dataset_id":
            array = pa.repeat(pa.scalar(dataset_id, type=field.type), len(frame))
        elif field.name == "max_ingested_at":
            array = pa.repeat(pa.scalar(ingested_at, type=field.type), len(frame))
        elif field.name not in frame.columns:
            if not field.nullable:
                raise ValueError(f"Missing required output column: {field.name}")
            array = pa.nulls(len(frame), type=field.type)
        elif pa.types.is_list(field.type):
            values = [normalize_leading_edge(value) for value in frame[field.name]]
            array = pa.array(values, type=field.type)
        else:
            array = pa.array(frame[field.name], type=field.type, from_pandas=False)
        arrays.append(array)
    return pa.Table.from_arrays(arrays, schema=schema)


def merge_parquet_files(
    paths: list[str],
    schema: pa.Schema,
    sort_columns: list[str],
    output_path: str | Path,
    dataset_id: str,
    ingested_at: pd.Timestamp,
    unique_columns: tuple[str, ...] = (),
) -> dict[str, Any]:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    unique_values = {column: set() for column in unique_columns}
    target_groups: dict[tuple[bool, str], list[int]] = {}
    row_count = 0

    with tempfile.TemporaryDirectory(
        prefix=f".{output.name}.", dir=output.parent
    ) as temporary_directory:
        temporary_directory = Path(temporary_directory)
        sorted_path = temporary_directory / "sorted_groups.parquet"
        temporary_output = temporary_directory / "merged.parquet"
        row_group_index = 0

        with pq.ParquetWriter(sorted_path, schema, compression="zstd") as writer:
            for path in paths:
                if not path:
                    continue
                # ponytail: one batch in RAM; stream row groups if batch_size grows.
                frame = pd.read_parquet(path)
                if frame.empty:
                    continue
                row_count += len(frame)
                for column, values in unique_values.items():
                    values.update(frame[column].dropna().unique())
                if "leading_edge" in schema.names and "leading_edge" not in frame:
                    raise KeyError("leading_edge")
                frame = frame.sort_values(sort_columns, kind="stable").reset_index(
                    drop=True
                )

                target_column = sort_columns[0]
                for target, group in frame.groupby(
                    target_column, sort=False, dropna=False
                ):
                    target_key = (
                        bool(pd.isna(target)),
                        "" if pd.isna(target) else str(target),
                    )
                    table = frame_to_table(group, schema, dataset_id, ingested_at)
                    writer.write_table(table, row_group_size=max(table.num_rows, 1))
                    target_groups.setdefault(target_key, []).append(row_group_index)
                    row_group_index += 1

        sorted_file = pq.ParquetFile(sorted_path)
        with pq.ParquetWriter(temporary_output, schema, compression="zstd") as writer:
            for _, row_groups in sorted(target_groups.items()):
                tables = [sorted_file.read_row_group(index) for index in row_groups]
                if len(tables) == 1:
                    table = tables[0]
                else:
                    # Prepared batches assign each target once; only duplicate
                    # targets are materialized together to restore their sort.
                    frame = pd.concat(
                        [table.to_pandas() for table in tables], ignore_index=True
                    )
                    frame = frame.sort_values(sort_columns, kind="stable").reset_index(
                        drop=True
                    )
                    table = frame_to_table(frame, schema, dataset_id, ingested_at)
                writer.write_table(table, row_group_size=max(table.num_rows, 1))
        os.replace(temporary_output, output)

    return {
        "rows": row_count,
        "unique_counts": {
            column: len(values) for column, values in unique_values.items()
        },
    }


def read_json(path: str | None) -> dict[str, Any]:
    if not path:
        return {}
    with open(path) as handle:
        return json.load(handle)


def read_metrics(paths: list[str]) -> list[dict[str, Any]]:
    metrics = []
    for path in paths:
        with open(path) as handle:
            metrics.append(json.load(handle))
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge batch DEA/GSEA Parquet outputs."
    )
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dea-files", nargs="*", default=[])
    parser.add_argument("--gsea-files", nargs="*", default=[])
    parser.add_argument("--metrics-files", nargs="*", default=[])
    parser.add_argument("--manifest", default="")
    parser.add_argument("--dea-output", default="")
    parser.add_argument("--gsea-output", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    max_ingested_at = pd.Timestamp.now(tz="UTC").floor("s")

    dea_output = outdir / (args.dea_output or f"{args.dataset_id}.dea.parquet")
    gsea_output = outdir / (args.gsea_output or f"{args.dataset_id}.gsea.parquet")
    dea_stats = merge_parquet_files(
        args.dea_files,
        DEA_SCHEMA,
        ["perturbed_target_symbol", "effect_gene_symbol"],
        dea_output,
        args.dataset_id,
        max_ingested_at,
        unique_columns=("perturbed_target_symbol",),
    )
    gsea_stats = merge_parquet_files(
        args.gsea_files,
        GSEA_SCHEMA,
        ["perturbed_target_symbol", "term"],
        gsea_output,
        args.dataset_id,
        max_ingested_at,
        unique_columns=("term",),
    )

    manifest = read_json(args.manifest)
    batch_metrics = read_metrics(args.metrics_files)
    summary = {
        "dataset_id": args.dataset_id,
        "max_ingested_at": max_ingested_at.isoformat(),
        "dea_output": str(dea_output),
        "gsea_output": str(gsea_output),
        "n_dea_rows": dea_stats["rows"],
        "n_gsea_rows": gsea_stats["rows"],
        "n_perturbations_in_dea": dea_stats["unique_counts"]["perturbed_target_symbol"],
        "n_terms_in_gsea": gsea_stats["unique_counts"]["term"],
        "manifest": manifest,
        "batch_metrics": batch_metrics,
    }
    with open(outdir / f"{args.dataset_id}.summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    print(f"DEA rows: {summary['n_dea_rows']:,} -> {dea_output}")
    print(f"GSEA rows: {summary['n_gsea_rows']:,} -> {gsea_output}")
    print(f"Shared max_ingested_at: {summary['max_ingested_at']}")


if __name__ == "__main__":
    main()
