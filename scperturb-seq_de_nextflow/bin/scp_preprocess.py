#!/usr/bin/env python3
import argparse
import anndata as ad
import pandas as pd
import math

p = argparse.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--others-per-file", type=int, default=99)
p.add_argument("--out", required=True)
args = p.parse_args()

A = ad.read_h5ad(args.input, backed="r")

obs = A.obs
obs["perturbed_target_symbol"] = obs["perturbed_target_symbol"].str.replace(
    r"^control.*", "control", regex=True
)
labels = obs["perturbed_target_symbol"].astype(str)

counts = labels.value_counts()
non_ctrl = counts.drop(index="control", errors="ignore")
perts = non_ctrl.index.tolist()
per_file = max(1, args.others_per_file)
n_batches = math.ceil(len(perts) / per_file)

rows = []
for i in range(n_batches):
    batch_id = f"batch{i+1:03d}"
    start, end = i * per_file, min((i + 1) * per_file, len(perts))
    chunk = perts[start:end]
    rows.append({"batch_id": batch_id, "control": "control", "perts": ",".join(chunk)})

pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
print(f"[BATCHES] Wrote {args.out} with {n_batches} batches")
