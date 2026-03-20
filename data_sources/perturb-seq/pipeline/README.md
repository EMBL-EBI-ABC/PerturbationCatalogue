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