# Perturb-seq DEA + GSEA Pipeline

This Nextflow pipeline runs differential expression and preranked GSEA from the
filtered H5AD produced by the comparison/guide-calling step.

It is intended for `nadig_2025_jurkat` and expects:

- `obs["called_knockout_gene_count"]`
- `obs["called_control_probe_count"]`
- `obs["perturbation_call_type"]`
- `obs["perturbed_target_symbol"]`
- symbol-based expression feature names in `var_names`, with `var["gene_symbol"]`
- raw counts in `X` by default

Re-run `data_sources/perturb-seq/comparison/comparison.py` before this pipeline
if the filtered H5AD was produced before the explicit non-targeting control
annotations were added.

Cells with no called guide are not controls. Controls are cells with one or more
called non-targeting probes and zero called gene-targeting probes. Perturbation
groups are cells with exactly one called gene and zero called non-targeting
probes. Cells with multiple gene calls, or mixed gene-targeting and
non-targeting calls, are excluded before DEA/GSEA.

## Method

For each perturbation, the pipeline compares all valid single-gene perturbation
cells against explicit non-targeting-only control cells using Scanpy's Wilcoxon rank-sum
implementation on library-size normalized, log1p-transformed counts.

GSEA is run with GSEApy prerank, using the Wilcoxon score as the ranking metric.
The `sidak` column is computed as `1 - (1 - pval) ^ n_terms` within each
perturbation's GSEA result table.

The pipeline batches perturbations across Slurm jobs. Each batch reads only the
control rows plus that batch's perturbation rows from the H5AD CSR matrix.

## Build Image

Build the image locally or in an interactive cluster session with Singularity:

```bash
cd $HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/dea-gsea
singularity build dea_gsea.sif Singularity.def
```

## Inputs

For `nadig_2025_jurkat` on the cluster:

```bash
export DATASET_ID=nadig_2025_jurkat
export H5AD=/hps/nobackup/mfreeberg/perturb_seq_fastq/results/nadig_2025_jurkat/experiment_final.filtered.h5ad
export OUTDIR=/hps/nobackup/mfreeberg/perturb_seq_fastq/results/nadig_2025_jurkat/dea_gsea
export GMT=$HPS_PATH/cache/msigdb/h.all.v2025.1.Hs.symbols.gmt
```

The filtered H5AD should be produced by the updated comparison/QC pipeline,
which writes gene symbols before this pipeline runs.

## Run

```bash
cd $HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/dea-gsea
module load nextflow/25.04.6
mkdir -p logs

time srun --mem=16G --time=7-00:00:00 --unbuffered \
  nextflow -log logs/${DATASET_ID}.dea_gsea.nextflow.log \
    run main.nf \
    -profile slurm,singularity \
    -name ${DATASET_ID}_dea_gsea \
    -work-dir work/${DATASET_ID} \
    --dataset_id $DATASET_ID \
    --h5ad $H5AD \
    --gmt $GMT \
    --outdir $OUTDIR \
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

The preparation step logs:

- total cells
- non-targeting-only control cells
- single-gene/no-control perturbation cells
- cells removed because more than one gene was perturbed
- cells removed because gene-targeting and non-targeting guides were mixed
- cells excluded because no guide was called
- retained perturbation count

## Output Schemas

DEA columns:

```text
dataset_id, perturbed_target_symbol, gene, padj, log2FoldChange,
score_name, score_value, cell_type, ingested_at
```

GSEA columns:

```text
dataset_id, term, perturbed_target_symbol, es, nes, pval, sidak, fdr,
geneset_size, leading_edge, cell_type, ingested_at
```

`leading_edge` is written as a Parquet list of strings.

## Smoke Test

For a fast Nextflow wiring test without GSEA, use a tiny batch size and skip
GSEA:

```bash
nextflow run main.nf \
  -profile standard \
  --dataset_id $DATASET_ID \
  --h5ad ./experiment_final.filtered.h5ad \
  --outdir /tmp/${DATASET_ID}_dea_gsea_smoke \
  --batch_size 1 \
  --limit_perturbations 2 \
  --skip_gsea true
```

Do not use `--skip_gsea true` for the production run.
