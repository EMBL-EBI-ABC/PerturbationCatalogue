#!/usr/bin/env python3
import argparse, anndata as ad

p = argparse.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--counts-layer", default="counts")  # use "X" if raw counts already in X
p.add_argument("--out", required=True)
a = p.parse_args()


A = ad.read_h5ad(a.input)

# Ensure required obs fields
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
    a.counts_layer
    and a.counts_layer != "X"
    and a.counts_layer in (A.layers.keys() if A.layers is not None else [])
):
    A.X = A.layers[a.counts_layer]

# Do NOT normalise or filter here; DE needs raw counts
A.write_h5ad(a.out)
print(f"[PREPROCESS] Wrote {a.out} :: cells={A.n_obs} genes={A.n_vars}")
