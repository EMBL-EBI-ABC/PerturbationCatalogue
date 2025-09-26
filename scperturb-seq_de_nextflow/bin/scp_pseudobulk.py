#!/usr/bin/env python3
import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import sys

p = argparse.ArgumentParser()
p.add_argument("--input", required=True, help="Path to pre-processed parquet file")
p.add_argument(
    "--pert-list",
    help="Comma-separated perturbations to include (non-control). If set, only these + control are kept.",
)
p.add_argument("--counts-layer", default="X")
p.add_argument(
    "--min-reps",
    type=int,
    default=3,
    help="Minimum replicate groups per perturbation to keep",
)
p.add_argument(
    "--replicate-mode",
    choices=["sample", "auto"],
    default="auto",
    help="Use biological samples, pseudo-replicate bags, or auto-select per perturbation",
)
p.add_argument(
    "--auto_min_cells_per_bag",
    type=int,
    default=20,
    help="Target bag size in auto mode",
)
p.add_argument(
    "--auto_min_bags", type=int, default=5, help="Lower bound for bags in auto mode"
)
p.add_argument(
    "--auto_max_bags", type=int, default=10, help="Upper bound for bags in auto mode"
)
p.add_argument("--random-seed", type=int, default=1, help="RNG seed for bagging")

# gene prefilter (independent filtering)
p.add_argument(
    "--min_gene_sum",
    type=int,
    default=10,
    help="Keep genes with total pseudobulk sum >= this",
)
p.add_argument(
    "--min_count_per_group",
    type=int,
    default=5,
    help="Keep genes with >= this in at least min_groups_expr groups",
)
p.add_argument(
    "--min_groups_expr",
    type=int,
    default=2,
    help="Number of groups required to meet min_count_per_group",
)
p.add_argument("--row-chunk", type=int, default=10_000, help="Rows per read chunk")
p.add_argument("--out-counts", required=True)
p.add_argument("--out-meta", required=True)
p.add_argument("--out-summary", required=True)
args = p.parse_args()

rng = np.random.default_rng(args.random_seed)
# ---------------------- Load (backed) ----------------------
try:
    A = ad.read_h5ad(args.input, backed="r")
except Exception as e:
    sys.exit(f"[PB-BACKED-FAST] Failed to read {args.input}: {e}")


obs = A.obs.copy()
if "perturbation" not in obs.columns:
    for alt in ("perturbed_target_symbol",):
        if alt in obs.columns:
            obs["perturbation"] = obs[alt].astype(str)
            break
if "perturbation" not in obs.columns:
    raise SystemExit(
        "Missing obs['perturbation']; tried perturbed_target_symbol as fallbacks."
    )

obs["perturbation"] = obs["perturbation"].str.replace(
    r"^control.*", "control", regex=True
)

if "sample_id" not in obs.columns:
    obs["sample_id"] = "S1"

A.obs = obs

# FIX: adata.var index structure.
A.var.drop(columns=["gene_symbol"], inplace=True)
A.var.index.name = "gene_symbol"
A.var.reset_index(inplace=True)
A.var.set_index("gene_symbol", inplace=True)
A.var["gene_symbol"] = A.var.index

# Put raw counts into X
if (
    args.counts_layer
    and args.counts_layer != "X"
    and args.counts_layer in (A.layers.keys() if A.layers is not None else [])
):
    A.X = A.layers[args.counts_layer]


obs = A.obs[["sample_id", "perturbation", "cell_type_label"]].copy()
obs["sample_id"] = obs["sample_id"].astype(str)
obs["perturbation"] = obs["perturbation"].astype(str)
obs["cell_type_label"] = obs["cell_type_label"].astype(str)

# ---------------------- Subset by pert-list ----------------------
if args.pert_list:
    keep_perts = set(s.strip() for s in args.pert_list.split(",") if s.strip())
    keep_mask = obs["perturbation"].isin(keep_perts.union({"control"}))
    sel_idx = np.where(keep_mask.values)[0]
else:
    sel_idx = np.arange(A.n_obs)

if sel_idx.size == 0:
    sys.exit(
        "[PB-BACKED-FAST] No rows selected after filtering; check --pert-list/--pert-column."
    )

# Restrict obs to selected rows
obs = obs.iloc[sel_idx].copy()


# Group assignment (sample / bag / auto), mirroring scp_pseudobulk.py
def auto_bags(nc: int) -> int:
    B = int(round(nc / max(1, args.auto_min_cells_per_bag)))
    B = int(np.clip(B, args.auto_min_bags, args.auto_max_bags))
    if nc >= args.auto_min_cells_per_bag * args.min_reps:
        B = max(B, args.min_reps)
    return max(1, B)


keys = pd.Series(index=obs.index, dtype="object")
modes = {}
summary = []
for pert, df in obs.groupby("perturbation", sort=False):
    idx = df.index.to_numpy()
    n_c = len(idx)
    n_s = df["sample_id"].nunique()
    if args.replicate_mode == "sample" or (
        args.replicate_mode == "auto" and n_s >= args.min_reps
    ):
        keys.loc[idx] = df["sample_id"].astype(str) + "::" + pert
        modes[pert] = "sample"
        kept = n_s
        total = n_s
    else:
        B = auto_bags(n_c)
        bag_ids = rng.integers(0, max(1, B), size=n_c)
        kept = 0
        for b in range(B):
            sel = idx[bag_ids == b]
            if sel.size >= args.auto_min_cells_per_bag:
                keys.loc[sel] = f"bag{b}::{pert}"
                kept += 1
        modes[pert] = f"bag[{kept}/{B}]"
        total = B
    summary.append(
        {
            "perturbation": pert,
            "cells": n_c,
            "unique_samples": int(n_s),
            "mode": modes[pert],
            "groups_kept": int(kept),
            "groups_planned": int(total),
        }
    )

keys = keys.dropna()
if keys.empty:
    sys.exit("[PB-BACKED-FAST] No groups formed.")

# Build group list and group id per cell (selected rows only)
groups = pd.Index(pd.unique(keys.values), name="group")
gid_map = pd.Series({g: i for i, g in enumerate(groups)}, dtype="int64")
cell_gid = keys.map(gid_map).to_numpy()  # len == number of selected cells

# Map selected obs labels to positional row indices in full matrix
obs_index_full = pd.Index(A.obs.index)
sel_pos = obs_index_full.get_indexer(keys.index)
if (sel_pos < 0).any():
    sys.exit("[PB-BACKED-FAST] Index alignment failed between obs and matrix rows.")

# ---------------------- Stream aggregation by row-chunks ----------------------
n_groups = len(groups)
n_genes = int(A.n_vars)
PB = np.zeros((n_groups, n_genes), dtype=np.int64)

# Sort positions so we read contiguous HDF5 blocks
order = np.argsort(sel_pos)
sel_pos_sorted = sel_pos[order]
cell_gid_sorted = cell_gid[order]

start = 0
N = sel_pos_sorted.size
while start < N:
    end = min(start + args.row_chunk, N)
    rows = sel_pos_sorted[start:end]

    # Read block (backed slice)
    block = A.X[rows, :]
    # Dense-ify block (works for both dense and scipy.sparse)
    block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)

    # Accumulate into group sums
    gids = cell_gid_sorted[start:end]
    np.add.at(PB, gids, block)

    start = end

# ---------------------- Meta & QC ----------------------
meta = pd.DataFrame(
    [g.split("::", 1) for g in groups], columns=["sample_id", "perturbation"]
)
meta.index = groups
meta["mode"] = [modes[p] for p in meta["perturbation"]]
meta["is_control"] = (meta["perturbation"] == "control").astype(int)

# Robust mode(cell_type_label) per group (alignment-safe)
align_idx = keys.index  # cells that contributed to groups
tmp = pd.DataFrame(
    {
        "group": keys.values,
        "cell_type_label": obs.loc[align_idx, "cell_type_label"].astype(str).values,
    }
)


def _mode_or_unspecified(s: pd.Series) -> str:
    m = s.mode()
    return m.iat[0] if not m.empty else "unspecified"


ct_map = tmp.groupby("group", sort=False)["cell_type_label"].apply(_mode_or_unspecified)
meta["cell_type_label"] = meta.index.map(ct_map).fillna("unspecified")

# Enforce min_reps per perturbation (count distinct sample_id across groups)
rep_counts = meta.groupby("perturbation")["sample_id"].nunique()
valid_perts = rep_counts[rep_counts >= args.min_reps].index
keep_groups = meta["perturbation"].isin(valid_perts).to_numpy()

PB = PB[keep_groups, :]
meta = meta.loc[keep_groups]

if meta.empty:
    sys.exit(
        "[PB-BACKED-FAST] No contrasts after min_reps filtering. Lower --min-reps or increase bag size."
    )

# ---------------------- Independent gene filtering (optional) ----------------------
genes = A.var_names.to_list()
counts_df = pd.DataFrame(PB.astype(np.int32), index=meta.index, columns=genes)

if args.min_gene_sum > 0 or args.min_groups_expr > 0:
    sums = counts_df.sum(axis=0)
    groups_expr = (counts_df >= args.min_count_per_group).sum(axis=0)
    keep_genes = (sums >= args.min_gene_sum) & (groups_expr >= args.min_groups_expr)
    kept_n = int(keep_genes.sum())
    if kept_n > 0:
        counts_df = counts_df.loc[:, keep_genes]
        print(
            f"[PB-BACKED-FAST] Gene filter kept {kept_n}/{len(genes)} genes.",
            file=sys.stderr,
        )
    else:
        print(
            "[PB-BACKED-FAST] Gene filter removed all genes; writing unfiltered matrix.",
            file=sys.stderr,
        )

# ---------------------- Write outputs ----------------------
out_counts = Path(args.out_counts)
out_meta = Path(args.out_meta)
out_summary = Path(args.out_summary)
out_counts.parent.mkdir(parents=True, exist_ok=True)
out_meta.parent.mkdir(parents=True, exist_ok=True)
out_summary.parent.mkdir(parents=True, exist_ok=True)

# Parquet writes
counts_df.to_parquet(out_counts, index=True)
meta.to_parquet(out_meta, index=True)

summary_df = pd.DataFrame(summary)
summary_df = summary_df.merge(
    rep_counts.rename("groups_before_minreps").reset_index(),
    on="perturbation",
    how="left",
)
summary_df.to_csv(out_summary, sep="\t", index=False)

print(f"[PB-BACKED-FAST] counts={counts_df.shape} -> {out_counts}")
print(f"[PB-BACKED-FAST] meta={meta.shape} -> {out_meta}")
print(f"[PB-BACKED-FAST] summary -> {out_summary}")
