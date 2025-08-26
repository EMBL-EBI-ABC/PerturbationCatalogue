#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp

p = argparse.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--min-reps", type=int, default=3)
p.add_argument("--replicate-mode", choices=["sample", "auto"], default="auto")
p.add_argument("--auto_min_cells_per_bag", type=int, default=20)
p.add_argument("--auto_min_bags", type=int, default=5)
p.add_argument("--auto_max_bags", type=int, default=10)
p.add_argument("--random-seed", type=int, default=1)
p.add_argument("--out-counts", required=True)
p.add_argument("--out-meta", required=True)
p.add_argument("--out-summary", required=True)
a = p.parse_args()

rng = np.random.default_rng(a.random_seed)
A = ad.read_h5ad(a.input)

obs = A.obs[["sample_id", "perturbation", "cell_type_label"]].copy()
obs["sample_id"] = obs["sample_id"].astype(str)
obs["perturbation"] = obs["perturbation"].astype(str)
obs["cell_type_label"] = obs["cell_type_label"].astype(str)

X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)  # (cells x genes)
cells_n, genes_n = X.shape


# --- group assignment ---
def bags_for(number_of_cells: int) -> int:
    bags_number = int(round(number_of_cells / max(1, a.auto_min_cells_per_bag)))
    bags_number = int(np.clip(bags_number, a.auto_min_bags, a.auto_max_bags))
    # ensure at least min_reps if there are enough cells to populate them
    if number_of_cells >= a.auto_min_cells_per_bag * a.min_reps:
        bags_number = max(bags_number, a.min_reps)
    return max(1, bags_number)


keys = pd.Series(index=obs.index, dtype="object")
modes = {}  # perturbation -> "sample" or "bag"
per_summ = []

for pert, df in obs.groupby("perturbation", sort=False):
    idx = df.index.to_numpy()
    n_cells = len(idx)
    n_samps = df["sample_id"].nunique()

    if a.replicate_mode == "sample" or (
        a.replicate_mode == "auto" and n_samps >= a.min_reps
    ):
        # sample mode
        keys.loc[idx] = df["sample_id"].astype(str) + "::" + pert
        modes[pert] = "sample"
    else:
        # bag mode
        bags = bags_for(n_cells)
        bag_ids = rng.integers(low=0, high=max(1, bags), size=n_cells)
        kept = 0
        for b in range(bags):
            sel = idx[bag_ids == b]
            if len(sel) >= a.auto_min_cells_per_bag:  # threshold for auto & bag
                keys.loc[sel] = f"bag{b}::{pert}"
                kept += 1
        modes[pert] = f"bag[{kept}/{bags}]"

    per_summ.append(
        {"perturbation": pert, "cells": n_cells, "unique_samples": int(n_samps)}
    )

keys = keys.dropna()
if keys.empty:
    raise SystemExit("[PSEUDOBULK] No groups formed; relax thresholds?")

# Align obs and build positional index
obs = obs.loc[keys.index]
pos = obs.index.get_indexer(keys.index)  # integer positions
X = X[pos, :]

# Group list & meta
groups = pd.Index(pd.unique(keys.values), name="group")
g_index = pd.Series({g: i for i, g in enumerate(groups)})
col_g = keys.map(g_index).to_numpy()

# Aggregate (groups x cells) @ (cells x genes) -> (groups x genes)
grp = sp.csr_matrix(
    (np.ones_like(col_g), (col_g, np.arange(len(col_g)))),
    shape=(len(groups), len(col_g)),
)
PB = (grp @ X).toarray().astype("int32")

meta = pd.DataFrame(
    [g.split("::", 1) for g in groups], columns=["sample_id", "perturbation"]
)
meta.index = groups
meta["mode"] = [modes[p] for p in meta["perturbation"]]
meta["is_control"] = (meta["perturbation"] == "control").astype(int)

# map group -> mode cell_type_label
# keys.index are cell indices per group; we have obs aligned to keys.index
group_to_cells = pd.Series(
    keys.index, index=keys.values
)  # group_label -> cell_index (many-to-one)
group_cells = {g: [] for g in meta.index}
for g, cell_idx in group_to_cells.items():
    group_cells[g].append(cell_idx)
# mode cell type per group
ct_mode = []
for g in meta.index:
    cells = group_cells.get(g, [])
    if len(cells) == 0:
        ct_mode.append("unspecified")
    else:
        cts = obs.loc[cells, "cell_type_label"].astype(str)
        ct_mode.append(cts.mode().iat[0] if not cts.empty else "unspecified")
meta["cell_type_label"] = ct_mode

# Filter to perturbations with >= min_reps replicates
rep_counts = meta.groupby("perturbation")["sample_id"].nunique()
valid_perts = rep_counts[rep_counts >= a.min_reps].index
keep = meta["perturbation"].isin(valid_perts).to_numpy()
meta = meta.loc[keep]
PB = PB[keep, :]

if meta.empty:
    raise SystemExit(
        "[PSEUDOBULK] After enforcing min_reps, no contrasts remain. Lower --min-reps or thresholds."
    )


# Outputs
genes = A.var_names.to_list()
counts = pd.DataFrame(PB, index=meta.index, columns=genes)
counts.to_parquet(a.out_counts, index=True)
meta.to_parquet(a.out_meta, index=True)

summary = pd.DataFrame(per_summ).merge(
    rep_counts.rename("groups").reset_index(), how="left", on="perturbation"
)
summary["mode"] = summary["perturbation"].map({k: v for k, v in modes.items()})
summary.to_csv(a.out_summary, sep="\t", index=False)
print(
    f"[PSEUDOBULK] counts={counts.shape} written; meta={meta.shape}; summary -> {a.out_summary}"
)
