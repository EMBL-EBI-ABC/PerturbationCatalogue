# Perturb-seq Raw Data Processing Pipeline

This Nextflow pipeline processes raw Perturb-seq FASTQ files (downloaded from ENA/SRA) into a **single, unified `h5ad` count matrix**. 

The pipeline is **sample-aware**: it uses ENA metadata to group multiple sequencing runs (SRRs) and lanes (L001-L004) into their original physical libraries. This prevents barcode collisions across separate 10x reactions and ensures that CRISPR guide RNA counts are correctly linked to the matching cell's mRNA profile.

## Requirements
- **Nextflow**: Version 25.04.6 or newer.
- **kb-python**: Installed via `pip install kb-python` or provided via Singularity.
- **Python Data Stack**: `pip install pandas openpyxl anndata`

## End-to-End Example: Processing the Nadig 2025 Jurkat Dataset (SAMN40972597)

### 1. Download the Reference Transcriptome and GTF
We use Ensembl GRCh38 (Release 111).

```bash
mkdir -p $HPS_PATH/cache/reference
cd $HPS_PATH/cache/reference

wget -q http://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -q http://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
```

### 2. Generate the Guide Whitelist (`features.tsv`)
Extract guide sequences from the authors' supplementary data.

```bash
cd $HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/pipeline

python3 generate_features_nadig.py \
  ../../../data_exploration/Perturbseq/supplementary/nadig_2025_guide_info.xlsx \
  $HPS_PATH/perturb_seq_fastq/SAMN40972597/features.tsv
```

### 3. Fetch ENA Metadata
This is required for the pipeline to correctly group FASTQs by physical sample.

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SAMN40972597&result=read_run&fields=run_accession,library_name,fastq_ftp&format=tsv" \
  > $HPS_PATH/perturb_seq_fastq/SAMN40972597/ena_metadata.tsv
```

### 4. Build the Custom Singularity Image
Before running the pipeline, build the Singularity image from the provided definition file locally, then upload it to the cluster to `$HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/pipeline/kb_python.sif`.

```bash
# On your local machine with sudo access
sudo singularity build kb_python.sif Singularity.def

# Upload to the cluster
scp kb_python.sif user@cluster:$HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/pipeline/
```

### 5. Run the Pipeline
The pipeline will automatically identify samples by their primary group (e.g., `8`) and process mRNA and sgRNA modalities in parallel before merging and concatenating with unique barcode suffixes (e.g., `BARCODE-8`).

```bash
module load nextflow/25.04.6

# Run the pipeline head process via srun
time srun --mem=16G --time=7-00:00:00 --unbuffered \
  nextflow run main.nf \
    -profile slurm,singularity \
    --fastq_dir $HPS_PATH/perturb_seq_fastq/SAMN40972597 \
    --metadata_tsv $HPS_PATH/perturb_seq_fastq/SAMN40972597/ena_metadata.tsv \
    --outdir $HPS_PATH/perturb_seq_fastq/results \
    --chemistry 10xv3 \
    --transcriptome_fa $HPS_PATH/cache/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --gtf $HPS_PATH/cache/reference/Homo_sapiens.GRCh38.111.gtf.gz \
    --features_tsv $HPS_PATH/perturb_seq_fastq/SAMN40972597/features.tsv
```

## Outputs
- `results/merged_samples/`: Individual H5AD files for each physical 10x well.
- `results/experiment_final.h5ad`: The final unified matrix (Gzip compressed).

## Understanding the "Sample ID" logic
The ENA libraries use the notation `jurkat_<modality>_<sample_group>_<sub_sample>_L<lane>`.
Example: `jurkat_mRNA_8_4_L004` vs `jurkat_sgRNA_8_1_L001`.

The pipeline identifies **`8`** as the unique Sample ID (the main well or condition pool). It aggregates all sub-samples (`8_1`, `8_4`) and lanes for that specific group and ensures mRNA and sgRNA are merged correctly for that physical pool of cells. During the final concatenation, barcodes are suffixed with `-8` to prevent collisions with other pools (e.g., `1`).
