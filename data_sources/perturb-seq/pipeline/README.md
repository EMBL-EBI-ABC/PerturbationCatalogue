# Unified Perturb-seq pipeline

This Nextflow pipeline downloads the SRA runs in a sample sheet, streams their
FASTQ records into kb-python, and produces raw and QC-filtered H5ADs, comparison
reports, probe calls, DEA and GSEA Parquet data products in one workflow.
Uncompressed FASTQs are buffered on disk and deleted after feeding the counter.

## Build and install the image

Build the image locally:

```bash
cd data_sources/perturb-seq/pipeline
singularity build --fakeroot perturb_seq.sif Singularity.def
```

The workflow runs with `perturb_seq.sif` built from `Singularity.def`.

Make the finished image available at
`$HPS_PATH/data_sources/perturb-seq/pipeline/perturb_seq.sif` in the execution
environment.

Provide a genome FASTA, its GTF, the dataset's guide-feature TSV, the curated
(author) H5AD for comparison, and a gene-set GMT for GSEA. The reference
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
TSV (20 bp sequence and probe name carrying the target) and a sample sheet with
columns `sample_id`, `mRNA_srrs`, and `sgRNA_srrs`. Each run list is
semicolon-separated. Put all runs for one physical sample/well in one row;
combine neither different wells nor different modalities. A dataset may add
`guide_feature_offset` and `guide_feature_length` to trim a captured guide from
a longer read; both default to zero.

For example, `python3 datasets/nadig_2025/generate_inputs.py` produces
`datasets/nadig_2025/features.tsv`, `jurkat_samples.tsv` and `hepg2_samples.tsv`.
Replogle inputs are under `datasets/replogle_2022/`, named
`<dataset_id>_features.tsv` and `<dataset_id>_samples.tsv`.
The K562 genome-wide sample sheet intentionally excludes the unavailable
sgRNA run `SRR19331204` from `KD8_17`; the exclusion is encoded in its input
generator as well as the committed sample sheet.

Adamson 2016 inputs are under `datasets/adamson_2016/`. Gasperini 2019 inputs
are under `datasets/gasperini_2019/`; run
`python3 datasets/gasperini_2019/generate_inputs.py` to regenerate its guide
reference and ENA-derived sample sheet.

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
  --dataset_id "$DATASET_ID" \
  --curated_h5ad "$HPS_PATH/perturb_seq_fastq/source_h5ad/$DATASET_ID.h5ad" \
  --gmt "$HPS_PATH/cache/msigdb/h.all.v2025.1.Hs.symbols.gmt" \
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

`stream_metrics.json` records per-run archive/FASTQ sizes, checksums, processing
times and spot counts for output verification. Command logs remain in
`stream_logs/`. Storage estimates must include extraction scratch and final
matrix/BUS intermediates as well as the bounded input buffers.

mRNA counts use bustools cell filtering. KITE counts use `counts_unfiltered`, then
align to the mRNA cell barcodes. The existing guide-barcode transformation
complements bases 8–9 before alignment. All merged samples are concatenated on
disk; the final H5AD is compressed with HDF5 gzip compression.

## QC, comparison, probe calling and DEA/GSEA

After raw H5AD compression, QC and Gaussian–Poisson probe calling produce
`experiment_final.filtered.h5ad` and curated-data comparison reports. Cell/gene
filtering and control assignment are described in [comparison](comparison/README.md).
The same GTF is used for counting and gene-symbol resolution.

DEA compares each eligible single-gene perturbation against the shared
non-targeting controls using normalized/log-transformed counts and Scanpy
Wilcoxon scores. GSEApy prerank uses those scores and the supplied GMT.
See [DEA/GSEA](dea-gsea/README.md) for methods and output schemas.

`--batch_size 50` and `--min_cells_per_perturbation 10` control analysis batching
and eligibility. `--limit_perturbations 0` analyzes every eligible perturbation.
GSEA defaults: 1,000 permutations, minimum/maximum gene-set sizes 15/500 and
seed 1. Existing analysis controls `--target_sum`, `--matrix_key`,
`--gene_map` and `--tie_correct` remain available.

Nextflow stages every downstream input from its producer, rather than reading
asynchronously published files. The filtered H5AD is shared by all analysis
batches. One `-resume` covers the entire workflow.

## Outputs and checks

- `counts_standard/<sample>/` and `counts_kite/<sample>/`: count matrices and streaming metrics.
- `merged_samples/`: per-sample H5AD files and guide-overlap diagnostics.
- `experiment_final.h5ad`: raw unified count matrix.
- `experiment_final.filtered.h5ad`: QC-filtered counts, guide calls and control annotations.
- `comparison_results/<dataset_id>/`: all comparison reports and plots; override the parent with `--comparison_outdir`.
- `dea_gsea/prep/analysis_inputs/`: manifest, gene metadata, controls and batch definitions.
- `dea_gsea/batch_results/`: per-batch DEA/GSEA Parquet files and metrics.
- `dea_gsea/<dataset_id>.{dea.parquet,gsea.parquet,summary.json}`: merged analysis products.
- Nextflow trace: task timing, resource use and completion status.

## Benchmarks

Clean full-run measurements are recorded for [Nadig Jurkat](benchmarks/nadig_2025_jurkat.md)
and [Replogle K562 GW](benchmarks/replogle_2022_k562_gw_normalized.md), including
wall time, resource footprint and final DEA/GSEA product sizes.
