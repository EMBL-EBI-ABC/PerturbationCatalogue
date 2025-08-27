#!/usr/bin/env python3
"""
Pseudobulk aggregation from AnnData Zarr with sample/bag/auto replicate logic.

- Loads .zarr (AnnData) lazily.
- Forms replicate groups per perturbation:
    * sample  : sample_id::perturbation when >= min_reps true samples exist
    * auto    : prefer sample; else pick K ~= cells/auto_min_cells_per_bag (clipped)
      and ensure at least min_reps bags if enough cells exist.
- Aggregates (groups x genes) out-of-core with Dask, writing dense Parquet.
- Writes a QC summary and carries group-level cell_type_label (mode) into metadata.
- Optional gene prefiltering (independent filtering) to boost DE power.

Outputs:
  --out-counts : Parquet (rows=groups, cols=genes; int32)
  --out-meta   : Parquet (index=groups; sample_id, perturbation, is_control, mode, cell_type_label)
  --out-summary: TSV with per-perturbation cell counts, #samples, chosen mode, groups formed
"""

import argparse
import numpy as np
import pandas as pd
import anndata as ad
import dask.array as da
import sys

# ---------------------- CLI ----------------------
p = argparse.ArgumentParser()
p.add_argument("--zarr", required=True, help="Path to AnnData Zarr store")

# replicate logic
p.add_argument(
    "--min-reps",
    type=int,
    default=3,
    help="Minimum replicate groups per perturbation to keep",
)
p.add_argument(
    "--replicate-mode",
    choices=["sample", "bag", "auto"],
    default="auto",
    help="Use biological samples, pseudo-replicate bags, or auto-select per perturbation",
)
p.add_argument(
    "--auto_min_cells_per_bag",
    type=int,
    default=200,
    help="Target bag size in auto mode",
)
p.add_argument(
    "--auto_min_bags", type=int, default=3, help="Lower bound for bags in auto mode"
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

# outputs
p.add_argument("--out-counts", required=True)
p.add_argument("--out-meta", required=True)
p.add_argument("--out-summary", required=True)
args = p.parse_args()

rng = np.random.default_rng(args.random_seed)

# ---------------------- Load AnnData Zarr lazily ----------------------
try:
    A = ad.read_zarr(args.zarr)  # obs/var in memory; X via zarr on disk
except Exception as e:
    sys.exit(f"[ZARR-PB] Failed to read Zarr store '{args.zarr}': {e}")

# AnnData Zarr stores X as (n_obs, n_vars) => (cells, genes)
X = da.from_array(A.X)  # dask array (cells x genes)

obs = A.obs.copy()
for col in ("sample_id", "perturbation"):
    if col not in obs.columns:
        sys.exit(f"[ZARR-PB] Missing obs['{col}'] in Zarr store.")
obs["sample_id"] = obs["sample_id"].astype(str)
obs["perturbation"] = obs["perturbation"].astype(str)

# optional covariate
if "cell_type_label" not in obs.columns:
    obs["cell_type_label"] = "unspecified"
else:
    obs["cell_type_label"] = obs["cell_type_label"].astype(str)

n_cells, n_genes = map(int, X.shape)


# ---------------------- Replicate logic (sample / bag / auto) ----------------------
def auto_choose_bags(nc: int) -> int:
    # B ~ cells / target_size, clipped; ensure >= min_reps if enough cells exist
    B = int(round(nc / max(1, args.auto_min_cells_per_bag)))
    B = int(np.clip(B, args.auto_min_bags, args.auto_max_bags))
    if nc >= args.auto_min_cells_per_bag * args.min_reps:
        B = max(B, args.min_reps)
    return max(1, B)


# keys: group label per cell (e.g., "S1::TP53" or "bag2::TP53")
keys = pd.Series(index=obs.index, dtype="object")
modes = {}  # perturbation -> "sample" or "bag[K_kept/...]"
summary_rows = []  # per-perturbation QC summary

for pert, df in obs.groupby("perturbation", sort=False):
    idx = df.index.to_numpy()
    n_c = len(idx)
    n_samp = df["sample_id"].nunique()

    if args.replicate_mode == "sample" or (
        args.replicate_mode == "auto" and n_samp >= args.min_reps
    ):
        # True sample replication
        keys.loc[idx] = df["sample_id"].astype(str) + "::" + pert
        modes[pert] = "sample"
        kept = n_samp
        total = n_samp
    else:
        # Bagging
        B = auto_choose_bags(n_c)
        bag_ids = rng.integers(low=0, high=max(1, B), size=n_c)
        kept = 0
        for b in range(B):
            mask = bag_ids == b
            if mask.sum() >= args.auto_min_cells_per_bag:
                keys.loc[idx[mask]] = f"bag{b}::{pert}"
                kept += 1
        modes[pert] = f"bag[{kept}/{B}]"
        total = B

    summary_rows.append(
        {
            "perturbation": pert,
            "cells": n_c,
            "unique_samples": int(n_samp),
            "mode": modes[pert],
            "groups_kept": int(kept),
            "groups_planned": int(total),
        }
    )

keys = keys.dropna()
if keys.empty:
    sys.exit("[ZARR-PB] No groups formed; relax thresholds or check metadata.")

# Align obs to selected cells (those that have a group key)
obs = obs.loc[keys.index]

# Build group labels and per-cell group ids (length = n_selected_cells)
group_labels = pd.Index(pd.unique(keys.values), name="group")
group_index = pd.Series({g: i for i, g in enumerate(group_labels)}, dtype="int64")
cell_group_ids = keys.map(group_index).to_numpy()  # len == n_selected_cells

# Subset the Dask array X to the selected rows ONLY, so mask lengths match
# Map selected obs labels to positional indices in the original AnnData
sel_pos = pd.Index(A.obs.index).get_indexer(keys.index)
if (sel_pos < 0).any():
    raise IndexError(
        "[ZARR-PB] Some selected cell indices were not found in A.obs.index"
    )
X_sel_all = X[sel_pos, :]  # (n_selected_cells x n_genes)

# ---------------------- Fast blockwise aggregation ----------------------
n_genes = int(A.n_vars)
n_groups = int(len(group_labels))

# Align selected rows
sel_pos = pd.Index(A.obs.index).get_indexer(keys.index)
if (sel_pos < 0).any():
    raise IndexError(
        "[ZARR-PB] Some selected cell indices were not found in A.obs.index"
    )

# Subset once and rechunk rows
X_sel_all = X[sel_pos, :]  # (n_selected_cells x n_genes)
X_sel_all = X_sel_all.rechunk((100000, -1))

# Per-row group ids (aligned to selected rows), chunked like X_sel_all rows
group_id = cell_group_ids.astype(np.int64)
gid_da = da.from_array(group_id, chunks=X_sel_all.chunks[0])


def sum_by_group(block, gid_block, n_groups, *, _np=np):
    """Sum rows in `block` by group ids in `gid_block`. Returns shape (1, n_groups, n_genes)."""
    # Dense-ify if sparse
    if hasattr(block, "toarray"):
        block = block.toarray()
    else:
        block = _np.asarray(block)
    gid_block = _np.asarray(gid_block, dtype=_np.int64).reshape(-1)
    out = _np.zeros((n_groups, block.shape[1]), dtype=_np.int64)
    _np.add.at(out, gid_block, block)
    return out[_np.newaxis, ...]  # add leading chunk axis


# Map per row-chunk -> (1, n_groups, genes_chunk)
partials = da.map_blocks(
    sum_by_group,
    X_sel_all,
    gid_da,
    n_groups,
    dtype=np.int64,
    new_axis=0,  # the first axis is the chunk axis of length 1 per block
    chunks=(1, n_groups, X_sel_all.chunks[1][0]),
)

# Reduce over the leading chunk axis -> (n_groups, n_genes)
PB = partials.sum(axis=0).compute().astype(np.int32)

# Always 2-D (defensive, in case n_groups == 1)
PB = np.atleast_2d(PB)

# Sanity check: rows of PB must equal number of groups in meta
if PB.shape[0] != len(group_labels):
    raise RuntimeError(
        f"[ZARR-PB] Aggregation shape mismatch: PB rows={PB.shape[0]} vs groups={len(group_labels)}"
    )

# ---------------------- Build meta and filter by min_reps ----------------------
meta = pd.DataFrame(
    [g.split("::", 1) for g in group_labels], columns=["sample_id", "perturbation"]
)
meta.index = group_labels
meta["mode"] = [modes[p] for p in meta["perturbation"]]
meta["is_control"] = (meta["perturbation"] == "control").astype(int)

# cell_type_label (mode per group)
ct_mode = []
for g in meta.index:
    cells_mask = keys == g
    cts = obs.loc[cells_mask, "cell_type_label"]
    ct_mode.append(cts.mode().iat[0] if not cts.empty else "unspecified")
meta["cell_type_label"] = ct_mode

# Enforce min_reps per perturbation at the *group* level
rep_counts = meta.groupby("perturbation")["sample_id"].nunique()
valid_perts = rep_counts[rep_counts >= args.min_reps].index
keep_groups = meta["perturbation"].isin(valid_perts).to_numpy()

# keep_groups must be boolean of length == PB.shape[0]
if keep_groups.dtype != bool:
    keep_groups = keep_groups.astype(bool)

if PB.shape[0] != keep_groups.shape[0]:
    raise RuntimeError(
        f"[ZARR-PB] keep_groups length {keep_groups.shape[0]} "
        f"does not match PB rows {PB.shape[0]}"
    )

meta = meta.loc[keep_groups]
PB = PB[keep_groups, :]

if meta.empty:
    sys.exit(
        "[ZARR-PB] After enforcing min_reps, no contrasts remain. Lower --min-reps or increase bag size."
    )

# ---------------------- Optional gene prefiltering ----------------------
# Independent filtering: keep genes with sufficient support across groups.
genes = A.var_names.to_list()
counts_df = pd.DataFrame(PB, index=meta.index, columns=genes).astype("int32")

if args.min_gene_sum > 0 or args.min_groups_expr > 0:
    sums = counts_df.sum(axis=0)
    groups_expr = (counts_df >= args.min_count_per_group).sum(axis=0)
    keep_genes = (sums >= args.min_gene_sum) & (groups_expr >= args.min_groups_expr)
    kept_n = int(keep_genes.sum())
    if kept_n == 0:
        print(
            "[ZARR-PB] Gene filter removed all genes; writing unfiltered matrix.",
            file=sys.stderr,
        )
    else:
        counts_df = counts_df.loc[:, keep_genes]
        print(
            f"[ZARR-PB] Gene filter kept {kept_n}/{len(genes)} genes.", file=sys.stderr
        )

# ---------------------- Write outputs ----------------------
counts_df.to_parquet(args.out_counts, index=True)
meta.to_parquet(args.out_meta, index=True)

summary_df = pd.DataFrame(summary_rows)
summary_df = summary_df.merge(
    rep_counts.rename("groups_after_filter").reset_index(),
    on="perturbation",
    how="left",
)
summary_df.to_csv(args.out_summary, sep="\t", index=False)

print(f"[ZARR-PB] counts={counts_df.shape} -> {args.out_counts}")
print(f"[ZARR-PB] meta={meta.shape} -> {args.out_meta}")
print(f"[ZARR-PB] summary -> {args.out_summary}")
