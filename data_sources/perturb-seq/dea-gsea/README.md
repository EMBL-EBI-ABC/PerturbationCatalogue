# Perturb-seq DEA + GSEA Pipeline

This Nextflow pipeline runs differential expression and preranked GSEA from the
filtered H5AD produced by the comparison/guide-calling step.

## Method

For each perturbation, the pipeline compares perturbation cells against controls
using Scanpy's Wilcoxon rank-sum implementation on library-size normalized,
log1p-transformed counts.

GSEA is run with GSEApy prerank, using the Wilcoxon score as the ranking metric.
The `sidak` column is computed as `1 - (1 - pval) ^ n_terms` within each
perturbation's GSEA result table.

The pipeline batches perturbations across Slurm jobs. Each batch reads only the
control rows plus that batch's perturbation rows from the H5AD CSR matrix.

## Build Image

Build the image locally or in an interactive cluster session with Singularity:

```bash
cd ${HPS_PATH}/PerturbationCatalogue/data_sources/perturb-seq/dea-gsea
singularity build --force dea_gsea.sif Singularity.def
```

## Download Gene Sets

```bash
mkdir -p ${HPS_PATH}/cache/msigdb
wget -O ${HPS_PATH}/cache/msigdb/h.all.v2025.1.Hs.symbols.gmt \
  https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt
```

## Run

Set dataset name, for example `DATASET_ID=nadig_2025_jurkat`, then run:

```bash
cd ${HPS_PATH}/PerturbationCatalogue/data_sources/perturb-seq/dea-gsea
module load nextflow/25.04.6
mkdir -p logs
time srun --mem=16G --time=7-00:00:00 --unbuffered \
  nextflow -log logs/${DATASET_ID}.dea_gsea.nextflow.log \
    run main.nf \
    -profile slurm,singularity \
    -name ${DATASET_ID}_dea_gsea \
    -work-dir work/${DATASET_ID} \
    --dataset_id ${DATASET_ID} \
    --h5ad ${HPS_PATH}/perturb_seq_fastq/results/${DATASET_ID}/experiment_final.filtered.h5ad \
    --gmt ${HPS_PATH}/cache/msigdb/h.all.v2025.1.Hs.symbols.gmt \
    --outdir ${HPS_PATH}/perturb_seq_fastq/results/${DATASET_ID}/dea_gsea \
    --batch_size 50 \
    --limit_perturbations 0 \
    --min_cells_per_perturbation 10 \
    --gsea_permutations 1000
```

Increase `--batch_size` to reduce repeated control-cell reads, or decrease it to
lower per-job memory. The default Slurm profile uses 4 CPUs and 64 GB per
`ANALYZE_BATCH` job. Leave `--limit_perturbations 0` for production; positive
values are only intended for smoke tests.

## Outputs

Final files are published directly under `--outdir`:

- `${DATASET_ID}.dea.parquet`
- `${DATASET_ID}.gsea.parquet`
- `${DATASET_ID}.summary.json`

Per-batch Parquet files and metrics are published under:

- `${OUTDIR}/batch_results/`

Preparation metadata is published under:

- `${OUTDIR}/prep/analysis_inputs/manifest.json`

The preparation step logs cell filtering and batching metrics.

## Output Schemas

DEA columns:

```text
dataset_id, perturbed_target_symbol, perturbed_target_ensg, effect_gene_symbol,
effect_gene_ensg, padj, log2foldchange, score_name, score_value, cell_type,
max_ingested_at
```

GSEA columns:

```text
dataset_id, term, perturbed_target_symbol, perturbed_target_ensg, es, nes,
pval, sidak, fdr, geneset_size, leading_edge, cell_type, max_ingested_at
```

`leading_edge` is written as a Parquet list of strings for the BigQuery
`REPEATED STRING` field.
