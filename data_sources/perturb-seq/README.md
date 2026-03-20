# Perturb-Seq

## Background

Perturb-Seq data is processed from a CSV file containing differential expression results.

Before processing, set the following environment secrets:

- WAREHOUSE_BUCKET
- ELASTIC_ENDPOINT

## Ingest metadata

```bash
# Process.
python3 process.py \
  --input-filename /tmp/perturb-seq.csv \
  --output-filename /tmp/perturb-seq.jsonl
gsutil -q -m rm -r "gs://${WAREHOUSE_BUCKET}/perturb-seq"
gsutil -q cp /tmp/perturb-seq.jsonl "gs://${WAREHOUSE_BUCKET}/perturb-seq/metadata.jsonl"

# Ingest into Elastic.
python3 ../elastic_load.py \
  --elastic-endpoint "${ELASTIC_ENDPOINT}" \
  --elastic-index "perturb-seq" \
  --jsonl-data "gs://${WAREHOUSE_BUCKET}/perturb-seq/metadata.jsonl" \
  --id-field "record_id" \
  --field-properties '{
    "study_id": {
      "type": "keyword"
    },
    "perturbation": {
      "type": "keyword"
    },
    "gene": {
      "type": "keyword"
    }
  }'
```

## Stats for the four files

These are the four studies curated and currently used:

```
Processing complete:
  Records written: 29476
  Records skipped: 2803366
  Filter criteria:
    padj <= 0.05
    log2fc < -1.0 or > 1.0
```

# Perturb-seq Raw Data Processing Pipeline

This Nextflow pipeline processes raw Perturb-seq FASTQ files (downloaded from ENA/SRA) into `h5ad` format count matrices using the robust and fast [kallisto-bustools (`kb-python`)](https://www.kallistobus.tools/) suite. The pipeline supports both standard Single-Cell RNA-seq (cDNA) analysis and Feature Barcode (CRISPR/guide RNA) quantification using the KITE workflow.

## Requirements
- **Nextflow**
- **Conda** or **Mamba** (the pipeline will create an environment with `kb-python` automatically, provided `conda` is accessible). Alternatively, ensure `kb-python` is installed and available in your `$PATH`.

## The "Whitelist of Probes" (Features List)
When processing CRISPR guides or other feature barcodes, you must use the **KITE** workflow. This workflow requires a "whitelist of probes", which is a simple tab-separated values (TSV) file mapping the guide name to its sequence.

**Format for `features.tsv`:**
```tsv
sgRNA_A    ATCGATCGATCGATCG
sgRNA_B    GCTAGCTAGCTAGCTA
```
*Note: The file should NOT contain a header. Column 1 is the feature ID (e.g., guide name), and Column 2 is the sequence.*

## Running the Pipeline on SLURM
The pipeline is pre-configured to run on SLURM clusters. Ensure you specify the `-profile slurm` flag.

### Example 1: Standard cDNA Processing
```bash
nextflow run main.nf \
    -profile slurm \
    --fastq_dir /path/to/perturb_seq_fastq/SAMN40972597 \
    --outdir /path/to/results/SAMN40972597_cDNA \
    --workflow standard \
    --chemistry 10x_v3 \
    --reads_pattern "*_{1,2,3}.fastq.gz" \
    --transcriptome_fa /path/to/transcriptome.fa \
    --gtf /path/to/annotation.gtf
```

### Example 2: KITE Workflow (CRISPR Guides)
```bash
nextflow run main.nf \
    -profile slurm \
    --fastq_dir /path/to/perturb_seq_fastq/SAMN40972597 \
    --outdir /path/to/results/SAMN40972597_guides \
    --workflow kite \
    --chemistry 10x_v3 \
    --reads_pattern "*_{1,2,3}.fastq.gz" \
    --features_tsv /path/to/features.tsv
```

### Important Parameters:
- `--fastq_dir`: Path to the directory containing downloaded `fastq.gz` files.
- `--reads_pattern`: Glob pattern to pair reads. SRA often splits reads into three files (`_1`, `_2`, `_3`). If your SRA download has three files but your chemistry expects only two (e.g. `10x_v3`), you MUST use `reads_pattern` to select only the two relevant read files containing the barcode and the biological sequence (e.g., `*_{2,3}.fastq.gz`). Passing 3 files to a 2-file chemistry will cause an error. If your SRA download has two files, leave it as default or `*_{1,2}.fastq.gz`.
- `--chemistry`: Single-cell chemistry version. E.g., `10x_v2`, `10x_v3`, `10x_v3_multi`.
- `--workflow`: Choose `standard` for cDNA or `kite` for guides/features.

## Outputs
The pipeline output directories are structured as follows:
- `results/reference/`: Contains the generated Kallisto index and `t2g` (transcript-to-gene) mapping files.
- `results/counts/<SRR_ID>/`: Contains the quantification outputs. The raw counts and the processed `adata.h5ad` format matrix can be found in the `counts_unfiltered/` subdirectory.
