#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


def decode_value(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def decode_array(values: np.ndarray) -> np.ndarray:
    return np.asarray([decode_value(value) for value in values], dtype=object)


def read_h5ad_column_full(handle: h5py.File, group: str, column: str) -> np.ndarray:
    key = f"{group}/{column}"
    if key not in handle:
        raise KeyError(f"Missing H5AD column: {key}")

    obj = handle[key]
    if isinstance(obj, h5py.Group):
        encoding = obj.attrs.get("encoding-type")
        if encoding != "categorical":
            raise ValueError(f"Unsupported group encoding for {key}: {encoding}")
        categories = decode_array(obj["categories"][:])
        codes = obj["codes"][:]
        out = np.empty(codes.shape[0], dtype=object)
        valid = codes >= 0
        out[~valid] = ""
        out[valid] = categories[codes[valid]]
        return out

    values = obj[:]
    if values.dtype.kind in {"O", "S", "U"}:
        return decode_array(values)
    return values


def read_h5ad_index(handle: h5py.File, group: str) -> np.ndarray:
    group_obj = handle[group]
    index_name = decode_value(group_obj.attrs["_index"])
    return read_h5ad_column_full(handle, group, index_name)


def h5ad_column_exists(handle: h5py.File, group: str, column: str) -> bool:
    return f"{group}/{column}" in handle


def parse_gtf_attributes(raw: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for item in raw.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        fields[key] = value.strip().strip('"')
    return fields


def load_gtf_gene_symbols(path: str | None) -> dict[str, str]:
    if not path:
        return {}

    opener = gzip.open if str(path).endswith(".gz") else open
    mapping: dict[str, str] = {}
    with opener(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name")
            if gene_id and gene_name:
                mapping[gene_id] = gene_name
                mapping.setdefault(gene_id.split(".", 1)[0], gene_name)
    return mapping


def load_tsv_gene_symbols(path: str | None) -> dict[str, str]:
    if not path:
        return {}

    df = pd.read_csv(path, sep="\t")
    normalized = {col.lower(): col for col in df.columns}
    id_col = next(
        (
            normalized[name]
            for name in [
                "gene_id",
                "ensembl_gene_id",
                "ensembl_id",
                "feature_id",
                "id",
            ]
            if name in normalized
        ),
        None,
    )
    symbol_col = next(
        (
            normalized[name]
            for name in ["gene_symbol", "symbol", "gene_name", "name"]
            if name in normalized
        ),
        None,
    )
    if id_col is None or symbol_col is None:
        raise ValueError(
            "Gene map TSV must contain gene id and symbol columns. "
            f"Observed columns: {list(df.columns)}"
        )

    mapping: dict[str, str] = {}
    for gene_id, symbol in zip(df[id_col].astype(str), df[symbol_col].astype(str)):
        gene_id = gene_id.strip()
        symbol = symbol.strip()
        if gene_id and symbol and symbol.lower() != "nan":
            mapping[gene_id] = symbol
            mapping.setdefault(gene_id.split(".", 1)[0], symbol)
    return mapping


def make_unique_gene_keys(gene_ids: list[str], symbols: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    keys: list[str] = []
    for gene_id, symbol in zip(gene_ids, symbols):
        base = symbol if symbol else gene_id
        count = seen.get(base, 0)
        if count == 0:
            key = base
        else:
            key = f"{base}__{gene_id}"
        while key in seen:
            count += 1
            key = f"{base}__{gene_id}__{count}"
        seen[base] = seen.get(base, 0) + 1
        seen[key] = 1
        keys.append(key)
    return keys


def write_gene_metadata(
    h5ad_path: Path,
    out_path: Path,
    gene_map_path: str | None,
    gtf_path: str | None,
) -> dict[str, Any]:
    tsv_mapping = load_tsv_gene_symbols(gene_map_path)
    gtf_mapping = load_gtf_gene_symbols(gtf_path)
    symbol_map = {**gtf_mapping, **tsv_mapping}

    with h5py.File(h5ad_path, "r") as handle:
        var_index = [str(value) for value in read_h5ad_index(handle, "var")]
        if h5ad_column_exists(handle, "var", "gene_id"):
            gene_ids = [
                str(value) for value in read_h5ad_column_full(handle, "var", "gene_id")
            ]
        else:
            gene_ids = var_index
        if h5ad_column_exists(handle, "var", "gene_symbol"):
            h5ad_symbols = [
                str(value)
                for value in read_h5ad_column_full(handle, "var", "gene_symbol")
            ]
        else:
            h5ad_symbols = []

    symbols = []
    for idx, gene_id in enumerate(gene_ids):
        h5ad_symbol = h5ad_symbols[idx].strip() if h5ad_symbols else ""
        symbol = (
            h5ad_symbol
            or symbol_map.get(gene_id)
            or symbol_map.get(gene_id.split(".", 1)[0])
            or (var_index[idx] if not var_index[idx].startswith("ENSG") else "")
            or gene_id
        )
        symbols.append(symbol)
    gene_keys = make_unique_gene_keys(gene_ids, symbols)
    mapped = sum(
        symbol != gene_id and not str(symbol).startswith("ENSG")
        for gene_id, symbol in zip(gene_ids, symbols)
    )

    df = pd.DataFrame(
        {
            "gene_index": np.arange(len(gene_ids), dtype=np.int32),
            "gene_id": gene_ids,
            "gene_symbol": symbols,
            "gene_key": gene_keys,
            "h5ad_var_name": var_index,
        }
    )
    df.to_parquet(out_path, index=False)
    return {
        "n_genes": int(len(gene_ids)),
        "n_gene_symbols_mapped": int(mapped),
        "gene_map_path": gene_map_path or None,
        "gtf_path": gtf_path or None,
    }


def assign_balanced_batches(
    perturbation_counts: pd.Series,
    batch_size: int,
) -> list[list[str]]:
    perturbations = perturbation_counts.sort_values(ascending=False).index.tolist()
    if not perturbations:
        return []

    n_batches = int(math.ceil(len(perturbations) / batch_size))
    batches: list[list[str]] = [[] for _ in range(n_batches)]
    batch_cells = [0 for _ in range(n_batches)]

    for perturbation in perturbations:
        candidates = [
            idx for idx, batch in enumerate(batches) if len(batch) < batch_size
        ]
        target_batch = min(candidates, key=lambda idx: batch_cells[idx])
        batches[target_batch].append(perturbation)
        batch_cells[target_batch] += int(perturbation_counts.loc[perturbation])

    return [sorted(batch) for batch in batches if batch]


def prepare_batches(args: argparse.Namespace) -> dict[str, Any]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    batches_dir = outdir / "batches"
    batches_dir.mkdir(exist_ok=True)

    h5ad_path = Path(args.h5ad)
    with h5py.File(h5ad_path, "r") as handle:
        n_cells = int(handle["X"].attrs["shape"][0])
        n_genes = int(handle["X"].attrs["shape"][1])
        target_labels = read_h5ad_column_full(handle, "obs", args.target_col)
        gene_counts = read_h5ad_column_full(handle, "obs", args.gene_count_col).astype(
            np.int64
        )
        if not h5ad_column_exists(handle, "obs", args.control_probe_count_col):
            raise KeyError(
                f"Missing required control annotation obs['{args.control_probe_count_col}']. "
                "Re-run data_sources/perturb-seq/comparison/comparison.py after the "
                "non-targeting control annotation update."
            )
        control_probe_counts = read_h5ad_column_full(
            handle, "obs", args.control_probe_count_col
        ).astype(np.int64)
        call_types = (
            read_h5ad_column_full(handle, "obs", args.call_type_col)
            if h5ad_column_exists(handle, "obs", args.call_type_col)
            else np.full(n_cells, "", dtype=object)
        )

    if (
        len(target_labels) != n_cells
        or len(gene_counts) != n_cells
        or len(control_probe_counts) != n_cells
    ):
        raise ValueError("Observation metadata lengths do not match X shape.")

    control_mask = (gene_counts == 0) & (control_probe_counts > 0)
    single_mask = (
        (gene_counts == 1) & (control_probe_counts == 0) & (target_labels != "")
    )
    multi_mask = gene_counts > 1
    mixed_gene_control_mask = (gene_counts > 0) & (control_probe_counts > 0)
    no_guide_mask = (gene_counts == 0) & (control_probe_counts == 0)
    classified_mask = (
        control_mask
        | single_mask
        | multi_mask
        | mixed_gene_control_mask
        | no_guide_mask
    )

    control_indices = np.flatnonzero(control_mask).astype(np.int64)
    single_indices = np.flatnonzero(single_mask).astype(np.int64)
    if control_indices.size == 0:
        raise ValueError(
            "No explicit non-targeting-only control cells found. Expected cells with "
            f"{args.gene_count_col} == 0 and {args.control_probe_count_col} > 0."
        )
    np.save(outdir / "control_indices.npy", control_indices)

    single_df = pd.DataFrame(
        {
            "cell_index": single_indices,
            "perturbed_target_symbol": target_labels[single_indices],
        }
    )
    perturbation_counts = single_df["perturbed_target_symbol"].value_counts()
    keep_counts = perturbation_counts[
        perturbation_counts >= args.min_cells_per_perturbation
    ]
    n_after_min_cell_filter = int(len(keep_counts))
    if args.limit_perturbations > 0:
        keep_counts = keep_counts.sort_values(ascending=False).head(
            args.limit_perturbations
        )
    dropped_small = perturbation_counts[
        perturbation_counts < args.min_cells_per_perturbation
    ].sort_values(ascending=False)

    perturbation_to_indices = {
        perturbation: group["cell_index"].astype(int).tolist()
        for perturbation, group in single_df.groupby(
            "perturbed_target_symbol", sort=True
        )
        if perturbation in keep_counts.index
    }

    batches = assign_balanced_batches(keep_counts, args.batch_size)
    batch_summaries: list[dict[str, Any]] = []
    for batch_idx, perturbations in enumerate(batches):
        batch_id = f"batch_{batch_idx:04d}"
        perturbation_indices = {
            perturbation: perturbation_to_indices[perturbation]
            for perturbation in perturbations
        }
        n_perturbation_cells = sum(len(v) for v in perturbation_indices.values())
        payload = {
            "batch_id": batch_id,
            "dataset_id": args.dataset_id,
            "perturbations": perturbations,
            "perturbation_indices": perturbation_indices,
            "n_control_cells": int(len(control_indices)),
            "n_perturbation_cells": int(n_perturbation_cells),
        }
        with open(batches_dir / f"{batch_id}.json", "w") as handle:
            json.dump(payload, handle)
        batch_summaries.append(
            {
                "batch_id": batch_id,
                "n_perturbations": int(len(perturbations)),
                "n_perturbation_cells": int(n_perturbation_cells),
            }
        )

    gene_summary = write_gene_metadata(
        h5ad_path=h5ad_path,
        out_path=outdir / "gene_metadata.parquet",
        gene_map_path=args.gene_map,
        gtf_path=args.gtf,
    )

    manifest = {
        "dataset_id": args.dataset_id,
        "h5ad_path": str(h5ad_path),
        "target_col": args.target_col,
        "gene_count_col": args.gene_count_col,
        "control_probe_count_col": args.control_probe_count_col,
        "call_type_col": args.call_type_col,
        "batch_size": int(args.batch_size),
        "min_cells_per_perturbation": int(args.min_cells_per_perturbation),
        "n_cells_total": int(n_cells),
        "n_genes": int(n_genes),
        "n_control_cells": int(len(control_indices)),
        "n_single_gene_perturbation_cells": int(len(single_indices)),
        "n_cells_removed_multiple_perturbations": int(multi_mask.sum()),
        "n_cells_removed_gene_control_mixed": int(mixed_gene_control_mask.sum()),
        "n_cells_excluded_no_called_guides": int(no_guide_mask.sum()),
        "n_cells_ignored_unassigned_noncontrol": int((~classified_mask).sum()),
        "perturbation_call_type_counts": {
            str(call_type): int(count)
            for call_type, count in pd.Series(call_types)
            .value_counts()
            .sort_index()
            .items()
        },
        "n_perturbations_before_min_cell_filter": int(len(perturbation_counts)),
        "n_perturbations_after_min_cell_filter": n_after_min_cell_filter,
        "limit_perturbations": int(args.limit_perturbations),
        "n_perturbations": int(len(keep_counts)),
        "n_perturbations_dropped_min_cell_filter": int(len(dropped_small)),
        "dropped_perturbations_min_cell_filter": {
            str(k): int(v) for k, v in dropped_small.items()
        },
        "gene_metadata": gene_summary,
        "batches": batch_summaries,
    }
    with open(outdir / "manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    print(f"Total cells: {manifest['n_cells_total']:,}")
    print(
        "Control cells (non-targeting guides only, 0 called genes): "
        f"{manifest['n_control_cells']:,}"
    )
    print(
        "Single-gene/no-control perturbation cells: "
        f"{manifest['n_single_gene_perturbation_cells']:,}"
    )
    print(
        "Cells removed because >1 gene was perturbed: "
        f"{manifest['n_cells_removed_multiple_perturbations']:,}"
    )
    print(
        "Cells removed because gene-targeting and control guides were mixed: "
        f"{manifest['n_cells_removed_gene_control_mixed']:,}"
    )
    print(
        "Cells excluded because no guide was called: "
        f"{manifest['n_cells_excluded_no_called_guides']:,}"
    )
    print(
        "Perturbations retained: "
        f"{manifest['n_perturbations']:,} "
        f"({manifest['n_perturbations_dropped_min_cell_filter']:,} dropped below "
        f"{args.min_cells_per_perturbation} cells)"
    )
    if args.limit_perturbations > 0:
        print(f"Perturbation limit applied: largest {args.limit_perturbations:,}")
    print(f"Batches: {len(batches):,}")
    print(
        "Gene symbols mapped: "
        f"{gene_summary['n_gene_symbols_mapped']:,}/{gene_summary['n_genes']:,}"
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare H5AD metadata and balanced perturbation batches."
    )
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--target-col", default="perturbed_target_symbol")
    parser.add_argument("--gene-count-col", default="called_knockout_gene_count")
    parser.add_argument(
        "--control-probe-count-col", default="called_control_probe_count"
    )
    parser.add_argument("--call-type-col", default="perturbation_call_type")
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--min-cells-per-perturbation", type=int, default=10)
    parser.add_argument(
        "--limit-perturbations",
        type=int,
        default=0,
        help="Keep only the largest N perturbations after filtering. Use 0 for no limit.",
    )
    parser.add_argument("--gene-map", default="")
    parser.add_argument("--gtf", default="")
    args = parser.parse_args()

    if args.batch_size < 1:
        parser.error("--batch-size must be >= 1")
    if args.min_cells_per_perturbation < 2:
        parser.error("--min-cells-per-perturbation must be >= 2")
    if args.limit_perturbations < 0:
        parser.error("--limit-perturbations must be >= 0")
    args.gene_map = args.gene_map or None
    args.gtf = args.gtf or None
    return args


def main() -> None:
    prepare_batches(parse_args())


if __name__ == "__main__":
    main()
