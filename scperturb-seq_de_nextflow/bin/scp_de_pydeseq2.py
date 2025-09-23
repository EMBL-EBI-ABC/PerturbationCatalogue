#!/usr/bin/env python3
import argparse, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

p = argparse.ArgumentParser()
p.add_argument("--counts", required=True)  # rows = groups, cols = genes
p.add_argument("--meta", required=True)  # index = same groups
p.add_argument("--covariates", default="")
p.add_argument("--out-prefix", default="de_")
p.add_argument("--summary", default="summary.tsv")
a = p.parse_args()

counts = pd.read_parquet(a.counts)  # (G x genes)
meta = pd.read_parquet(a.meta)  # (G x meta)
# align
counts = counts.loc[meta.index]

meta = meta.assign(perturbation=meta["perturbation"].astype("category"))
if "control" not in list(meta["perturbation"].cat.categories):
    raise SystemExit(
        f"[DE] control 'control' not found in meta. Categories: {list(meta['perturbation'].cat.categories)}"
    )

# order categories with control first
meta["perturbation"] = meta["perturbation"].cat.reorder_categories(
    ["control"] + [x for x in meta["perturbation"].cat.categories if x != "control"],
    ordered=True,
)

design_factors = ["perturbation"]
if a.covariates.strip():
    design_factors += [c.strip() for c in a.covariates.split(",") if c.strip()]

print(design_factors)

dds = DeseqDataSet(
    counts=counts.astype(int),
    metadata=meta,
    design_factors=design_factors,
    ref_level=["perturbation", "control"],
)
dds.deseq2()

perts = [p for p in meta["perturbation"].cat.categories if p != "control"]

rows = []
for pert in perts:
    contrast = ["perturbation", pert, "control"]
    try:
        stat = DeseqStats(dds, contrast=contrast, alpha=0.05, cooks_filter=True)
        stat.summary()
    except Exception as e:
        print(f"[DE] Skipping {pert}: {e}")
        continue

    res = stat.results_df.sort_values("padj", na_position="last").copy()
    outp = f"{a.out_prefix}{pert}.csv"
    res.to_csv(outp)
    rows.append(
        {
            "contrast": pert,
            "n_genes": int(res.shape[0]),
            "n_sig_padj_0.05": int((res["padj"] < 0.05).sum()),
        }
    )

if len(perts) == 0:
    res = pd.DataFrame(
        data={
            "baseMean": [],
            "log2FoldChange": [],
            "lfcSE": [],
            "stat": [],
            "pvalue": [],
            "padj": [],
        },
        index=[],
    )
    outp = f"{a.out_prefix}None.csv"
    res.to_csv(outp)

pd.DataFrame(rows).to_csv(a.summary, sep="\t", index=False)
print(f"[DE] Wrote {len(rows)} contrasts; summary -> {a.summary}")
