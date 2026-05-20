# Perturb-seq Raw Data Processing Pipeline

This Nextflow pipeline processes raw Perturb-seq FASTQ files downloaded from SRA into a **single, unified `h5ad` count matrix**.

The pipeline is indended to be run on the Slurm cluster.

## Set up (to be done once)

The commands below need to be run in the interactive session (`sinteractive`).

### 1. Download the Reference Transcriptome and GTF
```bash
mkdir -p $HPS_PATH/cache/reference
cd $HPS_PATH/cache/reference
wget -q http://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -q http://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
```

### 2. Clone repository on the cluster
```bash
cd $HPS_PATH
git clone https://github.com/EMBL-EBI-ABC/PerturbationCatalogue
# Switch to a relevant branch as necessary
```

### 3. Build and upload Singularity image
Before running the pipeline, build the Singularity image from the provided definition file locally, then upload it to the cluster to `$HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/pipeline/kb_python.sif`.

```bash
# On your local machine with sudo access
sudo singularity build kb_python.sif Singularity.def
# Then copy via scp to your local directory on the cluster

# On the cluster
mv ~/kb_python.sif $HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/pipeline/
```

## Run (for every individual dataset)

### 1. Prepare probe whitelist and sample sheet

For each dataset, add a dataset-specific script under `datasets/<dataset_or_group>/` that produces:

- `<dataset>_features.tsv`: headerless TSV with a 20 bp guide sequence and probe name.
- `<dataset>_samples.tsv`: TSV with `sample_id`, `mRNA_srrs`, and `sgRNA_srrs`. SRR lists are semicolon-separated. It is important to group all runs for a given sample into one row. The grouping is usually clear from the `library_name` field, which is a bit different from dataset to dataset, but includes a clearly identifiable sample name. It is this sample name which the dataset-specific script must extract. 

Then, run the script to generate the inputs (in this example for nadig_2025 dataset group):

```bash
cd $HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/pipeline
python3 datasets/nadig_2025/generate_inputs.py
```

This will produce the $FEATURES_PATH and $SAMPLE_SHEET_PATH files in the same directory as the script.

### 2. Set up pipeline parameters

#### nadig_2025_jurkat
```bash
# Dataset
export DATASET_ID=nadig_2025_jurkat
export FEATURES_PATH=datasets/nadig_2025/features.tsv
export SAMPLE_SHEET_PATH=datasets/nadig_2025/jurkat_samples.tsv
export FASTQ_DIR_PATH=$HPS_PATH/perturb_seq_fastq/SAMN40972597
export OUTPUT_DIR=$HPS_PATH/perturb_seq_fastq/results/$DATASET_ID
# Chemistry
export CHEMISTRY=10xv3
# Reference
export TRANSCRIPTOME_FA=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
export GTF=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz
```

#### nadig_2025_hepg2
```bash
# Dataset
export DATASET_ID=nadig_2025_hepg2
export FEATURES_PATH=datasets/nadig_2025/features.tsv
export SAMPLE_SHEET_PATH=datasets/nadig_2025/hepg2_samples.tsv
export FASTQ_DIR_PATH=$HPS_PATH/perturb_seq_fastq/SAMN40972598
export OUTPUT_DIR=$HPS_PATH/perturb_seq_fastq/results/$DATASET_ID
# Chemistry
export CHEMISTRY=10xv3
# Reference
export TRANSCRIPTOME_FA=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
export GTF=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz
```

#### replogle_2022_k562_essential_normalized
```bash
# Dataset
export DATASET_ID=replogle_2022_k562_essential_normalized
export FEATURES_PATH=datasets/replogle_2022/replogle_2022_k562_essential_normalized_features.tsv
export SAMPLE_SHEET_PATH=datasets/replogle_2022/replogle_2022_k562_essential_normalized_samples.tsv
export FASTQ_DIR_PATH=$HPS_PATH/perturb_seq_fastq/SAMN28561243
export OUTPUT_DIR=$HPS_PATH/perturb_seq_fastq/results/$DATASET_ID
# Chemistry
export CHEMISTRY=10xv3
# Reference
export TRANSCRIPTOME_FA=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
export GTF=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz
```

#### replogle_2022_rpe1_essential_normalized
```bash
# Dataset
export DATASET_ID=replogle_2022_rpe1_essential_normalized
export FEATURES_PATH=datasets/replogle_2022/replogle_2022_rpe1_essential_normalized_features.tsv
export SAMPLE_SHEET_PATH=datasets/replogle_2022/replogle_2022_rpe1_essential_normalized_samples.tsv
export FASTQ_DIR_PATH=$HPS_PATH/perturb_seq_fastq/SAMN28561244
export OUTPUT_DIR=$HPS_PATH/perturb_seq_fastq/results/$DATASET_ID
# Chemistry
export CHEMISTRY=10xv3
# Reference
export TRANSCRIPTOME_FA=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
export GTF=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz
```

#### replogle_2022_k562_gw_normalized
```bash
# Dataset
export DATASET_ID=replogle_2022_k562_gw_normalized
export FEATURES_PATH=datasets/replogle_2022/replogle_2022_k562_gw_normalized_features.tsv
export SAMPLE_SHEET_PATH=datasets/replogle_2022/replogle_2022_k562_gw_normalized_samples.tsv
export FASTQ_DIR_PATH=$HPS_PATH/perturb_seq_fastq/SAMN28561242
export OUTPUT_DIR=$HPS_PATH/perturb_seq_fastq/results/$DATASET_ID
# Chemistry
export CHEMISTRY=10xv3
# Reference
export TRANSCRIPTOME_FA=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
export GTF=$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz
```

### 3. Run the pipeline

```bash
cd $HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/pipeline
module load nextflow/25.04.6
time srun --mem=16G --time=7-00:00:00 --unbuffered \
  nextflow run main.nf \
    -profile slurm,singularity \
    --fastq_dir $FASTQ_DIR_PATH \
    --sample_sheet $SAMPLE_SHEET_PATH \
    --features_tsv $FEATURES_PATH \
    --transcriptome_fa $TRANSCRIPTOME_FA \
    --gtf $GTF \
    --chemistry $CHEMISTRY \
    --outdir $OUTPUT_DIR
```

## Details

### Processing and merging logic
KITE guide counts are merged from `counts_unfiltered`, then aligned to the filtered mRNA cell barcodes. This avoids independently filtering the guide barcode universe before mRNA/guide alignment. For 10xv3 chemistry, sgRNA-library cell barcodes are first corrected by complementing bases 8-9 before alignment to the mRNA barcodes; the per-sample diagnostics report both raw and corrected barcode overlap.

### Outputs
- `results/merged_samples/`: Individual H5AD files for each physical 10x well.
- `results/merged_samples/*_guide_diagnostics.json`: Per-sample KITE barcode overlap and guide UMI diagnostics.
- `results/experiment_final.h5ad`: The final unified matrix (Gzip compressed).
