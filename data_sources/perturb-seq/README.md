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

## End-to-End Example: Processing the Nadig 2025 Jurkat Dataset

The Nadig 2025 Jurkat dataset (`SAMN40972597`) uses a multiplexed CRISPRi library with two guides per cell. To process this using the KITE workflow, we must extract the individual guide sequences from the authors' supplementary Excel file and create a `features.tsv` whitelist.

### 1. Generate the Guide Whitelist (`features.tsv`)

The script `generate_features_nadig.py` is included in this folder to extract the sequences for you.

```bash
# Ensure you have pandas and openpyxl installed
pip install pandas openpyxl

# Generate the whitelist using the curated supplementary table (assuming you run this from the pipeline dir)
python3 generate_features_nadig.py ../../data_exploration/Perturbseq/supplementary/nadig_2025_guide_info.xlsx features.tsv
```

*Note: The generated `features.tsv` will be a headerless file with two columns: `<sgID>` and `<sequence>`.*

### 2. Run the Pipeline on the SLURM Cluster

Assuming the raw FASTQ files are downloaded to `$HPS_PATH/perturb_seq_fastq/SAMN40972597`, run the guide quantification (KITE workflow):

```bash
# Process CRISPR Guide RNAs
nextflow run main.nf \
    -profile slurm \
    --fastq_dir $HPS_PATH/perturb_seq_fastq/SAMN40972597 \
    --outdir $HPS_PATH/results/SAMN40972597_guides \
    --workflow kite \
    --chemistry 10x_v3 \
    --reads_pattern "*_{1,2}.fastq.gz" \
    --features_tsv features.tsv
```

If you also need to re-quantify the Gene Expression (cDNA) from the same or corresponding reads:

```bash
# Process Gene Expression (cDNA)
nextflow run main.nf \
    -profile slurm \
    --fastq_dir $HPS_PATH/perturb_seq_fastq/SAMN40972597 \
    --outdir $HPS_PATH/results/SAMN40972597_cDNA \
    --workflow standard \
    --chemistry 10x_v3 \
    --reads_pattern "*_{1,2}.fastq.gz" \
    --transcriptome_fa /path/to/human_transcriptome.fa \
    --gtf /path/to/human_annotation.gtf
```


# Perturb-seq Raw Data Processing Pipeline

This Nextflow pipeline processes raw Perturb-seq FASTQ files (downloaded from ENA/SRA) into `h5ad` format count matrices using the [kallisto-bustools (`kb-python`)](https://www.kallistobus.tools/) suite. The pipeline supports both standard Single-Cell RNA-seq (cDNA) analysis and Feature Barcode (CRISPR/guide RNA) quantification using the KITE workflow.

## Requirements
- **Nextflow**: Install via `wget -qO- https://get.nextflow.io | bash` or your cluster's module system.
- **kb-python**: The core kallisto-bustools wrapper. 
  - **Option 1 (Virtual Environment):** Create a Python virtual environment and run `pip install kb-python`. Ensure the `kb` command is in your `$PATH`.
  - **Option 2 (Singularity/Apptainer):** The pipeline includes a Singularity profile. If your cluster has Singularity, you can run the pipeline with `-profile slurm,singularity` and it will automatically pull and use the `kallistobustools/kb_python:latest` image.
- **Python Data Stack (for extraction script only)**: `pip install pandas openpyxl`

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

## End-to-End Example: Processing the Nadig 2025 Jurkat Dataset

The Nadig 2025 Jurkat dataset (`SAMN40972597`) uses a multiplexed CRISPRi library with two guides per cell. To process this using the KITE workflow, we must extract the individual guide sequences from the authors' supplementary Excel file and create a `features.tsv` whitelist.

### 1. Generate the Guide Whitelist (`features.tsv`)

The script `generate_features_nadig.py` is included in this folder to extract the sequences for you.

```bash
# Ensure you have pandas and openpyxl installed
pip install pandas openpyxl

# Generate the whitelist using the curated supplementary table (assuming you run this from the pipeline dir)
python3 generate_features_nadig.py ../../data_exploration/Perturbseq/supplementary/nadig_2025_guide_info.xlsx features.tsv
```

*Note: The generated `features.tsv` will be a headerless file with two columns: `<sgID>` and `<sequence>`.*

### 2. Run the Pipeline on the SLURM Cluster

Assuming the raw FASTQ files are downloaded to `$HPS_PATH/perturb_seq_fastq/SAMN40972597`, run the guide quantification (KITE workflow):

```bash
# Process CRISPR Guide RNAs
nextflow run main.nf \
    -profile slurm \
    --fastq_dir $HPS_PATH/perturb_seq_fastq/SAMN40972597 \
    --outdir $HPS_PATH/results/SAMN40972597_guides \
    --workflow kite \
    --chemistry 10x_v3 \
    --reads_pattern "*_{1,2}.fastq.gz" \
    --features_tsv features.tsv
```

If you also need to re-quantify the Gene Expression (cDNA) from the same or corresponding reads:

```bash
# Process Gene Expression (cDNA)
nextflow run main.nf \
    -profile slurm \
    --fastq_dir $HPS_PATH/perturb_seq_fastq/SAMN40972597 \
    --outdir $HPS_PATH/results/SAMN40972597_cDNA \
    --workflow standard \
    --chemistry 10x_v3 \
    --reads_pattern "*_{1,2}.fastq.gz" \
    --transcriptome_fa /path/to/human_transcriptome.fa \
    --gtf /path/to/human_annotation.gtf
```

# Perturb-seq Raw Data Processing Pipeline

This Nextflow pipeline processes raw Perturb-seq FASTQ files (downloaded from ENA/SRA) into a **single, unified `h5ad` count matrix**. 

In Perturb-seq, the FASTQ files often contain a mixture of standard Single-Cell RNA-seq (cDNA) reads and Feature Barcode (CRISPR/guide RNA) reads. This pipeline leverages the [kallisto-bustools (`kb-python`)](https://www.kallistobus.tools/) suite to process **both modalities simultaneously** across all sequencing runs (SRRs), and automatically merges the results. 

The final output is a single `experiment_final.h5ad` where:
- The standard cell-by-gene expression matrix is in the main `.X` layer.
- The cell-by-guide CRISPR counts are securely mapped to the exact same cells and stored in `.obsm['guides']`.
- The guide names/IDs are stored in `.uns['guide_names']`.

## Requirements
- **Nextflow**: Install via `wget -qO- https://get.nextflow.io | bash` or your cluster's module system.
- **kb-python**: The core kallisto-bustools wrapper. 
  - **Option 1 (Virtual Environment):** Create a Python virtual environment and run `pip install kb-python`. Ensure the `kb` command is in your `$PATH`.
  - **Option 2 (Singularity/Apptainer):** The pipeline includes a Singularity profile. If your cluster has Singularity, you can run the pipeline with `-profile slurm,singularity` and it will automatically pull and use the `kallistobustools/kb_python:latest` image.
- **Python Data Stack (for extraction script only)**: `pip install pandas openpyxl`

## The "Whitelist of Probes" (Features List)
When processing CRISPR guides, the pipeline's internal KITE workflow requires a "whitelist of probes". This is a simple tab-separated values (TSV) file mapping the guide name to its sequence.

**Format for `features.tsv`:**
```tsv
sgRNA_A    ATCGATCGATCGATCG
sgRNA_B    GCTAGCTAGCTAGCTA
```
*Note: The file should NOT contain a header. Column 1 is the feature ID (e.g., guide name), and Column 2 is the sequence.*

## End-to-End Example: Processing the Nadig 2025 Jurkat Dataset

The Nadig 2025 Jurkat dataset (`SAMN40972597`) uses a multiplexed CRISPRi library with two guides per cell. To process this, we must extract the individual guide sequences from the authors' supplementary Excel file and create a `features.tsv` whitelist.

### 1. Generate the Guide Whitelist (`features.tsv`)

The script `generate_features_nadig.py` is included in this folder to extract the sequences.

```bash
# Ensure you have pandas and openpyxl installed
pip install pandas openpyxl

# Generate the whitelist using the curated supplementary table (assuming you run this from the pipeline dir)
python3 generate_features_nadig.py ../../../data_exploration/Perturbseq/supplementary/nadig_2025_guide_info.xlsx features.tsv
```

*Note: The generated `features.tsv` will be a headerless file with two columns: `<sgID>` and `<sequence>`.*

### 2. Run the Unified Pipeline on the SLURM Cluster

Assuming the raw FASTQ files are downloaded to `$HPS_PATH/perturb_seq_fastq/SAMN40972597`, run the unified dual-modality pipeline:

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --fastq_dir $HPS_PATH/perturb_seq_fastq/SAMN40972597 \
    --outdir $HPS_PATH/results/SAMN40972597 \
    --chemistry 10x_v3 \
    --reads_pattern "*_{1,2}.fastq.gz" \
    --transcriptome_fa /path/to/human_transcriptome.fa \
    --gtf /path/to/human_annotation.gtf \
    --features_tsv features.tsv
```

### Important Parameters:
- `--fastq_dir`: Path to the directory containing downloaded `fastq.gz` files.
- `--reads_pattern`: Glob pattern to pair reads. SRA often splits reads into three files (`_1`, `_2`, `_3`). If your SRA download has three files but your chemistry expects only two (e.g. `10x_v3`), you MUST use `reads_pattern` to select only the two relevant read files containing the barcode and the biological sequence (e.g., `*_{2,3}.fastq.gz`). 
- `--chemistry`: Single-cell chemistry version. E.g., `10x_v2`, `10x_v3`, `10x_v3_multi`.
- `--transcriptome_fa` / `--gtf`: Reference genome files for the standard expression matrix.
- `--features_tsv`: Whitelist mapping guides for the KITE matrix.

## Outputs
- `results/reference/standard/`: cDNA Kallisto index.
- `results/reference/kite/`: CRISPR Guide Kallisto index.
- `results/experiment_final.h5ad`: The fully combined, merged matrix containing all cells across all FASTQs, with both gene expression and guide assignments.