# Perturb-seq DEA/GSEA analysis

These Python analysis scripts are used by the
[unified Nextflow pipeline](../pipeline/README.md), after QC and guide calling.
They run in the pipeline's single `../pipeline/perturb_seq.sif` image; this
directory has no separate container.

## Method

For each perturbation, the pipeline compares perturbation cells against controls
using Scanpy's Wilcoxon rank-sum implementation on library-size normalized,
log1p-transformed counts.

GSEA is run with GSEApy prerank, using the Wilcoxon score as the ranking metric.
The `sidak` column is computed as `1 - (1 - pval) ^ n_terms` within each
perturbation's GSEA result table.

The pipeline batches perturbations across Slurm jobs. Each batch reads only the
control rows plus that batch's perturbation rows from the H5AD CSR matrix.

## Download Gene Sets

```bash
mkdir -p ${HPS_PATH}/cache/msigdb
wget -O ${HPS_PATH}/cache/msigdb/h.all.v2025.1.Hs.symbols.gmt \
  https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt
```

## Run

Use the unified pipeline entry point and its `--gmt`, `--batch_size`,
`--min_cells_per_perturbation` and GSEA options. The workflow preserves preparation
metadata, batch products and merged results under `${OUTDIR}/dea_gsea/`.

Increase `--batch_size` to reduce repeated control-cell reads, or decrease it to
lower per-job memory. The default Slurm profile uses 4 CPUs and 64 GB per
`ANALYZE_BATCH` job. Leave `--limit_perturbations 0` for production; positive
values are only intended for smoke tests.

## Outputs

Final files are published under `${OUTDIR}/dea_gsea/`:

- `${DATASET_ID}.dea.parquet`
- `${DATASET_ID}.gsea.parquet`
- `${DATASET_ID}.summary.json`

Per-batch Parquet files and metrics are published under:

- `${OUTDIR}/dea_gsea/batch_results/`

Preparation metadata is published under:

- `${OUTDIR}/dea_gsea/prep/analysis_inputs/manifest.json`

The preparation step logs cell filtering and batching metrics.

Output schemas are defined in `bin/io_schemas.py`.
