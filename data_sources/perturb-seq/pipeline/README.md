# Perturb-seq SRA download and counting pipeline

This Nextflow pipeline downloads the SRA runs in a sample sheet, streams their
FASTQ records into kb-python, and produces a unified H5AD count matrix. Raw
uncompressed FASTQs are buffered on disk and deleted after feeding the counter. Comparison, QC, probe assignment
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

Each sample/modality has one task with three overlapping stages: sequential SRA
archive retrieval using 32 parallel curl HTTPS byte ranges per file,
`fasterq-dump` extraction to an uncompressed FASTQ file containing all
technical and biological reads, and feeding one persistent `kb count --inleaved`
process. An archive is deleted immediately after successful extraction; the
FASTQ file is deleted after successful routing into the counter pipe.

Capacity is reserved before downloading/extracting. Each task has at most one
waiting archive in addition to the archive being extracted, and at most one
waiting FASTQ set in addition to the set being fed. In-progress downloads and
fasterq extraction scratch also consume storage. A completed FASTQ buffer prevents
further extraction; an occupied archive buffer prevents another download.

The NCBI locator supplies the full-quality SRA URL, size and MD5. Every range's
Content-Range and length is checked; the assembled archive must match the
published size and MD5 before extraction. Curl retries transient errors with
bounded timeouts. Range parts are removed as they are assembled, so assembly
adds at most one range's size to the archive footprint. No prefetch download or
SRA Lite substitution is used. Curl 7.68 or later is required for parallel transfers.

The C++ reader is compiled with `g++` in the task environment and processes files
using buffered I/O. It verifies FASTQ structure, sequence/quality lengths,
consecutive spot IDs, unique increasing read IDs and a fixed layout throughout
each accession. It selects exactly one 20–40 bp barcode/UMI read and one biological
read longer than 40 bp; short index reads are validated and consumed. Other
layouts require an explicit reader/chemistry update. Its spot/read totals must
match fasterq's extraction summary, and kallisto's processed-pair total must match
all runs. UMI deduplication and sample-level counting combine all runs.

The SLURM profile has no pipeline concurrency cap (`executor.queueSize=0`); SLURM
schedules tasks against account resources. Each count task requests 16 CPUs and
32 GB RAM, with four extraction threads, ten kb threads and capacity for the
reader and downloader. Smaller local allocations reduce those thread counts.
Download/count failures terminate the workflow. Completed Nextflow tasks can be
reused with `-resume`; an incomplete task needs its runs downloaded again.

`stream_events.jsonl` records stage starts/completions and samples task/buffer
disk bytes every five seconds (the larger of logical and allocated size). `stream_status.json` holds an atomic
current snapshot, including completed runs and per-stage timing. Final
`stream_metrics.json` includes per-run archive/FASTQ sizes, processing times,
spot counts and sampled disk peaks. Sampling can miss short peaks; extraction
scratch and final matrix/BUS intermediates must be included in storage estimates.

mRNA counts use bustools cell filtering. KITE counts use `counts_unfiltered`, then
align to the mRNA cell barcodes. The existing guide-barcode transformation
complements bases 8–9 before alignment. All merged samples are concatenated on
disk; the final H5AD is compressed with HDF5 gzip compression.

## Outputs and checks

- `counts_standard/<sample>/` and `counts_kite/<sample>/`: count matrices and streaming metrics.
- `merged_samples/`: per-sample H5AD files and guide-overlap diagnostics.
- `experiment_final.h5ad`: unified count matrix.
- Nextflow trace: task timing, resource use and completion status.

Run `python3 -B test_stream_count.py` for native reader, bounded-stage overlap
and failure-cleanup checks. Stream-to-counter
integration should also be compared with ordinary paired-input counting using
the same references and runs before changing the streaming implementation.
