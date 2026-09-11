# Perturb-seq SRA download and counting pipeline

This Nextflow pipeline downloads the SRA runs in a sample sheet, streams their
FASTQ records into kb-python, and produces a unified H5AD count matrix. FASTQ
files are neither stored nor gzip-compressed. Comparison, QC, probe assignment
and DEA/GSEA are separate workflows.

## Dependencies and references

Use Nextflow 25.04.6 and the container built from `Singularity.def`, which includes
kb-python 0.30.2 and SRA Toolkit 3.4.1. The `slurm,singularity` profiles use
`kb_python.sif` in this directory. Build the image where Singularity builds are
supported, and place it directly in the pipeline directory using your site's
approved transfer procedure.

Provide a genome FASTA, its GTF, and the dataset's guide-feature TSV. The reference
example uses Ensembl release 115:

```bash
mkdir -p "$HPS_PATH/cache/reference"
cd "$HPS_PATH/cache/reference"
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
```

Downloads, image installation and pipeline execution must follow your cluster's
allocation and storage policies. Keep HOME, caches and temporary directories in
allocated project storage.

## Sample sheet

Dataset preparation scripts under `datasets/` produce a headerless guide-feature
TSV (20 bp sequence and probe name carrying the target ENSG) and a sample sheet
with columns `sample_id`, `mRNA_srrs`, `sgRNA_srrs`. Each SRR list is
semicolon-separated. Put all runs for one physical sample/well in one row;
combine neither different wells nor different modalities.

For example, `python3 datasets/nadig_2025/generate_inputs.py` produces
`datasets/nadig_2025/features.tsv`, `jurkat_samples.tsv` and `hepg2_samples.tsv`.
Replogle inputs are under `datasets/replogle_2022/`, named
`<dataset_id>_features.tsv` and `<dataset_id>_samples.tsv`.

## Run Jurkat

Run the Nextflow controller inside a SLURM allocation, with its own writable HOME
and TMPDIR under project storage. From this pipeline directory:

```bash
module load nextflow/25.04.6
DATASET_ID=nadig_2025_jurkat
mkdir -p logs
nextflow -log "logs/${DATASET_ID}.nextflow.log" run main.nf \
  -profile slurm,singularity \
  -work-dir "work/$DATASET_ID" \
  -with-trace "logs/${DATASET_ID}.trace.tsv" \
  --sample_sheet datasets/nadig_2025/jurkat_samples.tsv \
  --features_tsv datasets/nadig_2025/features.tsv \
  --transcriptome_fa "$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
  --gtf "$HPS_PATH/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz" \
  --chemistry 10xv3 \
  --outdir "$HPS_PATH/perturb_seq_fastq/results/$DATASET_ID"
```

`--sra_bin /path/to/sratoolkit/bin` optionally uses an external SRA Toolkit
installation; make that path visible inside the container. Otherwise the tools
are found on PATH. `--limit N` selects the first N sample groups for small tests.

## Streaming and processing

Each sample/modality has one download/count task. It prefetches at most two SRA
archives: the current run and the next run. `fasterq-dump --split-spot --stdout
--include-technical` sends all reads through a pipe. The parser groups records
by spot and selects one 20–40 bp barcode/UMI read and one biological read longer
than 40 bp, matching this pipeline's supported read layout. Other reads, such as
short sample indexes, are consumed without being forwarded. Malformed records,
ambiguous or changing layouts, missing reads and duplicate spot/read IDs fail the
task. Supporting other layouts requires an explicit parser/chemistry update.

Pairs enter a single `kb count --inleaved` process through standard input. Sample
barcode correction, UMI deduplication and counting happen after every run has
been streamed. Consumed archives are deleted; SRA extraction scratch and BUS
intermediates still require disk space. Each task verifies that kallisto's
processed-pair count matches the producer's total and writes per-run counts and
timestamps to `stream_metrics.json`.

The SLURM profile allows eight concurrent count tasks per modality, each with 16
CPUs and 32 GB RAM. Two CPUs are assigned to extraction, one to read routing and
the remainder to kb. Download/count failures terminate the workflow; tasks are
not automatically retried. Successful Nextflow tasks can be reused with
`-resume`; an incomplete task needs its runs downloaded again.

mRNA counts use bustools cell filtering. KITE counts use `counts_unfiltered`, then
align to the mRNA cell barcodes. The existing guide-barcode transformation
complements bases 8–9 before alignment. All merged samples are concatenated on
disk; the final H5AD is compressed with HDF5 gzip compression.

## Outputs and checks

- `counts_standard/<sample>/` and `counts_kite/<sample>/`: count matrices and streaming metrics.
- `merged_samples/`: per-sample H5AD files and guide-overlap diagnostics.
- `experiment_final.h5ad`: unified count matrix.
- Nextflow trace: task timing, resource use and completion status.

Run `python3 -B test_stream_count.py` for parser checks. Stream-to-counter
integration should also be compared with ordinary paired-input counting using
the same references and runs before changing the streaming implementation.
