# Perturb-seq Raw Data Processing Pipeline

This Nextflow pipeline processes raw Perturb-seq FASTQ files (downloaded from ENA/SRA) into a **single, unified `h5ad` count matrix**. 

In Perturb-seq, the FASTQ files often contain a mixture of standard Single-Cell RNA-seq (cDNA) reads and Feature Barcode (CRISPR/guide RNA) reads. This pipeline leverages the [kallisto-bustools (`kb-python`)](https://www.kallistobus.tools/) suite to process **both modalities simultaneously** across all sequencing runs (SRRs), and automatically merges the results. 

The final output is a single `experiment_final.h5ad` where:
- The standard cell-by-gene expression matrix is in the main `.X` layer.
- The cell-by-guide CRISPR counts are securely mapped to the exact same cells and stored in `.obsm['guides']`.
- The guide names/IDs are stored in `.uns['guide_names']`.

## Requirements
- **Nextflow**: On the cluster, load the module via `module load nextflow/25.04.6`.
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

The Nadig 2025 Jurkat dataset (`SAMN40972597`) uses a multiplexed CRISPRi library with two guides per cell. To process this, we must download the reference human transcriptome, extract the guide sequences from the authors' supplementary data, and then run the unified pipeline.

The following steps provide exact, copy-pasteable commands with no placeholders.

### 1. Download the Reference Transcriptome and GTF
We will use the standard Ensembl GRCh38 (Release 111) for mapping cDNA reads. These files will be stored in `$HPS_PATH/cache/reference`.

```bash
mkdir -p $HPS_PATH/cache/reference
cd $HPS_PATH/cache/reference

# Download the FASTA and GTF files
wget -q http://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -q http://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
```

### 2. Generate the Guide Whitelist (`features.tsv`)
The script `generate_features_nadig.py` extracts the individual guide sequences from the authors' supplementary Excel file.

```bash
cd $HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/pipeline

# Install required python packages
pip install pandas openpyxl

# Generate the whitelist (saving to the dataset directory)
python3 generate_features_nadig.py ../../../data_exploration/Perturbseq/supplementary/nadig_2025_guide_info.xlsx $HPS_PATH/perturb_seq_fastq/SAMN40972597/features.tsv
```

### 3. Run the Unified Pipeline on the SLURM Cluster
Assuming the raw FASTQ files are downloaded to `$HPS_PATH/perturb_seq_fastq/SAMN40972597`, you can now run the pipeline. 

*(Note: The pipeline automatically inspects the 10x FASTQ triplet files per SRR and dynamically detects which is the barcode read and which is the biological read based on their internal sequence lengths. You no longer need to specify read patterns.)*

```bash
# Load Nextflow module
module load nextflow/25.04.6

# Run the pipeline
nextflow run main.nf \
    -profile slurm,singularity \
    --fastq_dir $HPS_PATH/perturb_seq_fastq/SAMN40972597 \
    --outdir $HPS_PATH/perturb_seq_fastq/results \
    --chemistry 10x_v3 \
    --transcriptome_fa $HPS_PATH/cache/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --gtf $HPS_PATH/cache/reference/Homo_sapiens.GRCh38.111.gtf.gz \
    --features_tsv $HPS_PATH/perturb_seq_fastq/SAMN40972597/features.tsv
```

## Outputs
- `results/reference/standard/`: cDNA Kallisto index.
- `results/reference/kite/`: CRISPR Guide Kallisto index.
- `results/experiment_final.h5ad`: The fully combined, merged matrix containing all cells across all FASTQs, with both gene expression and guide assignments.