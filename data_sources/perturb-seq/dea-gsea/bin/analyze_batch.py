#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from io_schemas import DEA_SCHEMA, GSEA_SCHEMA, write_parquet


CONTROL_LABEL = "__control__"
GROUP_COL = "dea_group"


def decode_value(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def decode_array(values: np.ndarray) -> list[str]:
    return [decode_value(value) for value in values]


def read_h5ad_index_selected(
    handle: h5py.File,
    group: str,
    indices: np.ndarray,
) -> list[str]:
    group_obj = handle[group]
    index_name = decode_value(group_obj.attrs["_index"])
    obj = handle[f"{group}/{index_name}"]
    values = obj[indices]
    if values.dtype.kind in {"O", "S", "U"}:
        return decode_array(values)
    return [str(value) for value in values]


def contiguous_runs(indices: np.ndarray) -> list[tuple[int, int]]:
    if indices.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(indices) != 1) + 1
    starts = np.r_[0, breaks]
    ends = np.r_[breaks, indices.size]
    return [
        (int(indices[start]), int(indices[end - 1]) + 1)
        for start, end in zip(starts, ends)
    ]


def read_csr_rows(
    handle: h5py.File,
    matrix_key: str,
    row_indices: np.ndarray,
    dtype: np.dtype,
) -> sp.csr_matrix:
    if matrix_key not in handle:
        raise KeyError(f"Matrix key not found in H5AD: {matrix_key}")
    group = handle[matrix_key]
    if group.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError(
            f"Only CSR matrices are supported. {matrix_key} has "
            f"encoding {group.attrs.get('encoding-type')}"
        )

    n_cols = int(group.attrs["shape"][1])
    indptr_ds = group["indptr"]
    indices_ds = group["indices"]
    data_ds = group["data"]

    data_chunks: list[np.ndarray] = []
    index_chunks: list[np.ndarray] = []
    new_indptr = [0]
    nnz = 0

    for start, end in contiguous_runs(row_indices):
        ptr = indptr_ds[start : end + 1].astype(np.int64, copy=False)
        data_start = int(ptr[0])
        data_end = int(ptr[-1])
        counts = np.diff(ptr)
        if data_end > data_start:
            data_chunks.append(data_ds[data_start:data_end].astype(dtype, copy=False))
            index_chunks.append(
                indices_ds[data_start:data_end].astype(np.int32, copy=False)
            )
        cumulative = nnz + np.cumsum(counts)
        new_indptr.extend(cumulative.astype(np.int64).tolist())
        nnz = int(cumulative[-1]) if cumulative.size else nnz

    if data_chunks:
        data = np.concatenate(data_chunks)
        col_indices = np.concatenate(index_chunks)
    else:
        data = np.array([], dtype=dtype)
        col_indices = np.array([], dtype=np.int32)

    return sp.csr_matrix(
        (data, col_indices, np.asarray(new_indptr, dtype=np.int64)),
        shape=(int(row_indices.size), n_cols),
    )


def load_gmt(path: str | Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    with open(path, "r") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = []
            seen = set()
            for gene in parts[2:]:
                gene = gene.strip()
                if gene and gene not in seen:
                    genes.append(gene)
                    seen.add(gene)
            if genes:
                gene_sets[term] = genes
    return gene_sets


def build_batch_adata(
    h5ad_path: Path,
    matrix_key: str,
    control_indices: np.ndarray,
    batch: dict[str, Any],
    gene_metadata: pd.DataFrame,
    float_dtype: str,
) -> ad.AnnData:
    perturbation_arrays = [
        np.asarray(batch["perturbation_indices"][perturbation], dtype=np.int64)
        for perturbation in batch["perturbations"]
    ]
    selected = np.unique(np.concatenate([control_indices, *perturbation_arrays]))
    selected.sort()

    dtype = np.dtype(float_dtype)
    with h5py.File(h5ad_path, "r") as handle:
        matrix = read_csr_rows(handle, matrix_key, selected, dtype=dtype)
        barcodes = read_h5ad_index_selected(handle, "obs", selected)

    labels = np.full(selected.size, CONTROL_LABEL, dtype=object)
    for perturbation, cell_indices in zip(batch["perturbations"], perturbation_arrays):
        positions = np.searchsorted(selected, cell_indices)
        labels[positions] = perturbation

    categories = [CONTROL_LABEL] + list(batch["perturbations"])
    obs = pd.DataFrame(
        {
            GROUP_COL: pd.Categorical(labels, categories=categories),
            "cell_index": selected.astype(np.int64),
        },
        index=barcodes,
    )
    var = gene_metadata.set_index("gene_key", drop=False).copy()
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    return adata


def finite_series(values: Any, fill: float) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    arr[~np.isfinite(arr)] = fill
    return arr


def rank_genes_groups_to_dea(
    adata: ad.AnnData,
    dataset_id: str,
    perturbations: list[str],
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    rg = adata.uns["rank_genes_groups"]
    available_groups = set(rg["names"].dtype.names)
    gene_lookup = adata.var["gene_symbol"].astype(str).to_dict()

    dea_frames: list[pd.DataFrame] = []
    ranking_by_perturbation: dict[str, pd.DataFrame] = {}

    for perturbation in perturbations:
        if perturbation not in available_groups:
            continue

        gene_keys = np.asarray(rg["names"][perturbation]).astype(str)
        genes = pd.Series(gene_keys).map(gene_lookup).fillna(pd.Series(gene_keys))
        scores = finite_series(rg["scores"][perturbation], 0.0)
        if "pvals_adj" in rg:
            padj = finite_series(rg["pvals_adj"][perturbation], 1.0)
        elif "pvals" in rg:
            padj = finite_series(rg["pvals"][perturbation], 1.0)
        else:
            padj = np.ones_like(scores, dtype=np.float64)

        if "logfoldchanges" in rg:
            logfc = finite_series(rg["logfoldchanges"][perturbation], 0.0)
        else:
            logfc = np.zeros_like(scores, dtype=np.float64)

        df = pd.DataFrame(
            {
                "dataset_id": dataset_id,
                "perturbed_target_symbol": perturbation,
                "gene": genes.astype(str).to_numpy(),
                "padj": padj,
                "log2foldchange": logfc,
                "score_name": "Wilcoxon Score",
                "score_value": scores,
                "cell_type": None,
            }
        )
        dea_frames.append(df)

        ranking = df[["gene", "score_value"]].copy()
        ranking = ranking[ranking["gene"].astype(str).str.len() > 0]
        ranking["_abs_score"] = ranking["score_value"].abs()
        ranking = (
            ranking.sort_values(["_abs_score", "score_value"], ascending=[False, False])
            .drop_duplicates(subset=["gene"], keep="first")
            .sort_values("score_value", ascending=False)
            .drop(columns="_abs_score")
        )
        ranking_by_perturbation[perturbation] = ranking

    if dea_frames:
        dea = pd.concat(dea_frames, ignore_index=True)
    else:
        dea = pd.DataFrame(columns=[field.name for field in DEA_SCHEMA])
    return dea, ranking_by_perturbation


def split_gene_field(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, float) and math.isnan(value):
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if str(item)]
    text = str(value).strip()
    if not text:
        return []
    delimiter = ";" if ";" in text else ","
    return [item.strip() for item in text.split(delimiter) if item.strip()]


def result_value(row: pd.Series, candidates: list[str], default: Any = None) -> Any:
    for candidate in candidates:
        if candidate in row.index:
            return row[candidate]
    return default


def run_gsea_for_perturbation(
    perturbation: str,
    ranking: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    args: argparse.Namespace,
) -> pd.DataFrame:
    import gseapy as gp

    if ranking.shape[0] < args.gsea_min_size:
        return pd.DataFrame(columns=[field.name for field in GSEA_SCHEMA])

    seed_offset = int(hashlib.sha1(perturbation.encode("utf-8")).hexdigest()[:8], 16)
    seed = int(args.gsea_seed + seed_offset % 1_000_000)
    rank_input = ranking.rename(columns={"score_value": "score"})[["gene", "score"]]

    pre_res = gp.prerank(
        rnk=rank_input,
        gene_sets=gene_sets,
        outdir=None,
        min_size=args.gsea_min_size,
        max_size=args.gsea_max_size,
        permutation_num=args.gsea_permutations,
        weight=1.0,
        ascending=False,
        threads=args.threads,
        seed=seed,
        no_plot=True,
        verbose=False,
    )
    res = pre_res.res2d.copy()
    if res.empty:
        return pd.DataFrame(columns=[field.name for field in GSEA_SCHEMA])

    rank_genes = set(rank_input["gene"].astype(str))
    rows: list[dict[str, Any]] = []
    for _, row in res.iterrows():
        term = str(result_value(row, ["Term", "term"], row.name))
        raw = pre_res.results.get(term, {})
        pval = float(result_value(row, ["NOM p-val", "pval", "P-value"], 1.0))
        if not math.isfinite(pval):
            pval = 1.0

        leading_edge = split_gene_field(
            raw.get("lead_genes")
            or raw.get("leading_edge")
            or result_value(row, ["Lead_genes", "lead_genes"], "")
        )
        matched = split_gene_field(
            raw.get("matched_genes")
            or raw.get("matched genes")
            or result_value(row, ["Matched genes", "matched_genes"], "")
        )
        if not matched:
            matched = sorted(set(gene_sets.get(term, [])) & rank_genes)

        rows.append(
            {
                "dataset_id": args.dataset_id,
                "term": term,
                "perturbed_target_symbol": perturbation,
                "es": float(result_value(row, ["ES", "es"], 0.0) or 0.0),
                "nes": float(result_value(row, ["NES", "nes"], 0.0) or 0.0),
                "pval": pval,
                "sidak": float(1.0 - math.pow(1.0 - pval, len(res))),
                "fdr": float(
                    result_value(row, ["FDR q-val", "fdr", "Adjusted P-value"], 1.0)
                    or 1.0
                ),
                "geneset_size": int(len(matched)),
                "leading_edge": leading_edge,
                "cell_type": None,
            }
        )

    out = pd.DataFrame(rows)
    for col in ["es", "nes", "pval", "sidak", "fdr"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out[["es", "nes"]] = out[["es", "nes"]].fillna(0.0)
    out[["pval", "sidak", "fdr"]] = out[["pval", "sidak", "fdr"]].fillna(1.0)
    return out


def run_gsea(
    ranking_by_perturbation: dict[str, pd.DataFrame],
    args: argparse.Namespace,
) -> pd.DataFrame:
    gene_sets = load_gmt(args.gmt)
    frames: list[pd.DataFrame] = []
    for perturbation, ranking in ranking_by_perturbation.items():
        try:
            result = run_gsea_for_perturbation(perturbation, ranking, gene_sets, args)
        except Exception as exc:
            print(f"[WARN] GSEA failed for {perturbation}: {exc}")
            result = pd.DataFrame(columns=[field.name for field in GSEA_SCHEMA])
        if not result.empty:
            frames.append(result)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=[field.name for field in GSEA_SCHEMA])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run DEA and GSEA for one batch.")
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--batch-json", required=True)
    parser.add_argument("--control-indices", required=True)
    parser.add_argument("--gene-metadata", required=True)
    parser.add_argument("--gmt", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--matrix-key", default="X")
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument(
        "--float-dtype", default="float32", choices=["float32", "float64"]
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--tie-correct", action="store_true")
    parser.add_argument("--gsea-permutations", type=int, default=1000)
    parser.add_argument("--gsea-min-size", type=int, default=15)
    parser.add_argument("--gsea-max-size", type=int, default=500)
    parser.add_argument("--gsea-seed", type=int, default=1)
    args = parser.parse_args()

    if not args.gmt:
        parser.error("--gmt is required")
    if args.gsea_permutations < 1:
        parser.error("--gsea-permutations must be >= 1")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    return args


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.batch_json) as handle:
        batch = json.load(handle)
    batch_id = batch["batch_id"]
    perturbations = list(batch["perturbations"])
    control_indices = np.load(args.control_indices)
    gene_metadata = pd.read_parquet(args.gene_metadata)

    print(
        f"{batch_id}: {len(perturbations)} perturbations, "
        f"{len(control_indices):,} controls, "
        f"{batch['n_perturbation_cells']:,} perturbation cells"
    )

    adata = build_batch_adata(
        h5ad_path=Path(args.h5ad),
        matrix_key=args.matrix_key,
        control_indices=control_indices,
        batch=batch,
        gene_metadata=gene_metadata,
        float_dtype=args.float_dtype,
    )
    print(f"{batch_id}: loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(
        adata,
        groupby=GROUP_COL,
        groups=perturbations,
        reference=CONTROL_LABEL,
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        tie_correct=args.tie_correct,
        use_raw=False,
        n_genes=adata.n_vars,
    )

    dea, ranking_by_perturbation = rank_genes_groups_to_dea(
        adata, args.dataset_id, perturbations
    )
    gsea = run_gsea(ranking_by_perturbation, args)
    max_ingested_at = pd.Timestamp.now(tz="UTC").floor("s")
    if not dea.empty:
        dea["max_ingested_at"] = max_ingested_at
    if not gsea.empty:
        gsea["max_ingested_at"] = max_ingested_at

    dea_path = outdir / f"{batch_id}.dea.parquet"
    gsea_path = outdir / f"{batch_id}.gsea.parquet"
    metrics_path = outdir / f"{batch_id}.metrics.json"
    write_parquet(dea, dea_path, DEA_SCHEMA)
    write_parquet(gsea, gsea_path, GSEA_SCHEMA)

    metrics = {
        "batch_id": batch_id,
        "n_controls": int(len(control_indices)),
        "n_perturbations": int(len(perturbations)),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_dea_rows": int(len(dea)),
        "n_gsea_rows": int(len(gsea)),
        "max_ingested_at": max_ingested_at.isoformat(),
    }
    with open(metrics_path, "w") as handle:
        json.dump(metrics, handle, indent=2)

    print(f"{batch_id}: wrote {len(dea):,} DEA rows and " f"{len(gsea):,} GSEA rows")


if __name__ == "__main__":
    main()
