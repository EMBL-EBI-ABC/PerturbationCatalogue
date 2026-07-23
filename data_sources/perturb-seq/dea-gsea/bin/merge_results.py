#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from io_schemas import DEA_SCHEMA, GSEA_SCHEMA, write_parquet


def empty_frame(schema) -> pd.DataFrame:
    return pd.DataFrame(columns=[field.name for field in schema])


def read_parquet_many(paths: list[str], schema) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path:
            continue
        df = pd.read_parquet(path)
        if not df.empty:
            frames.append(df)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return empty_frame(schema)


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

    dea = read_parquet_many(args.dea_files, DEA_SCHEMA)
    gsea = read_parquet_many(args.gsea_files, GSEA_SCHEMA)
    max_ingested_at = pd.Timestamp.now(tz="UTC").floor("s")

    if not dea.empty:
        dea["dataset_id"] = args.dataset_id
        dea["max_ingested_at"] = max_ingested_at
        dea = dea.sort_values(
            ["perturbed_target_symbol", "effect_gene_symbol"]
        ).reset_index(drop=True)
    else:
        dea = empty_frame(DEA_SCHEMA)

    if not gsea.empty:
        gsea["dataset_id"] = args.dataset_id
        gsea["max_ingested_at"] = max_ingested_at
        gsea["leading_edge"] = gsea["leading_edge"].apply(normalize_leading_edge)
        gsea = gsea.sort_values(["perturbed_target_symbol", "term"]).reset_index(
            drop=True
        )
    else:
        gsea = empty_frame(GSEA_SCHEMA)

    dea_output = outdir / (args.dea_output or f"{args.dataset_id}.dea.parquet")
    gsea_output = outdir / (args.gsea_output or f"{args.dataset_id}.gsea.parquet")
    write_parquet(dea, dea_output, DEA_SCHEMA)
    write_parquet(gsea, gsea_output, GSEA_SCHEMA)

    manifest = read_json(args.manifest)
    batch_metrics = read_metrics(args.metrics_files)
    summary = {
        "dataset_id": args.dataset_id,
        "max_ingested_at": max_ingested_at.isoformat(),
        "dea_output": str(dea_output),
        "gsea_output": str(gsea_output),
        "n_dea_rows": int(len(dea)),
        "n_gsea_rows": int(len(gsea)),
        "n_perturbations_in_dea": (
            int(dea["perturbed_target_symbol"].nunique()) if not dea.empty else 0
        ),
        "n_terms_in_gsea": int(gsea["term"].nunique()) if not gsea.empty else 0,
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
