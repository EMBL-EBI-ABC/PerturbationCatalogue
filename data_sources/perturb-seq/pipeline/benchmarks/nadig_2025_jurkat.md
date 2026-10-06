# Nadig Jurkat unified-pipeline benchmark

Clean full run on a Slurm cluster on 19 September 2026 using main-repository commit
`e1b2c84725b615c5865fa0efb55d700648ef459a`, Nextflow 25.04.6 and the
pipeline-wide image `perturb_seq.sif` (SHA-256
`b72e25f35aeeef3da9841103004dd54fadb349ff44cf12149600d28996e73e07`). The
dataset contains 56 samples and 1,792 accessions.

## Results

| Measurement | Result |
| --- | ---: |
| Pipeline wall time, input ready to final products | **1h 28m 08.989s** (5,288.989 s) |
| Fresh source H5AD download | 451.9 s, 8 ranged workers, 9,366,490,264 bytes |
| Successful / failed / retried tasks | **222 / 0 / 0** |
| Peak running / pending tasks | 57 / 48 |
| Peak allocated CPUs | 906 |
| Maximum sampled footprint | **1.459 TB / 1.327 TiB** |
| Largest task MaxRSS | **15.42 GiB**, `MERGE_RESULTS` |
| Final published output footprint | 72.499 GB |
| DEA/GSEA analysis batches | 47 |

The controller ran from **2026-09-19 01:12:58.704937 UTC** to
**02:41:07.694273 UTC**. It covers the clean unified workflow from the verified
curated H5AD through counting, modality merges, sample concatenation, H5AD
compression, QC/comparison/probe calling, input preparation, all DEA/GSEA
batches and final publication. The source H5AD transfer is recorded separately
because it completed before the controller started.

The profiler sampled every five seconds from 01:13:11.321854 to 02:41:11.105611
UTC. The peak includes Nextflow work, published outputs, comparison reports and
controller files; it does not include the immutable source H5AD. Final tracked
sizes were 1,374,834,122,974 bytes of work, 72,498,706,906 bytes of published
output, 2,580,480 bytes of comparison reports and 32,169,319 bytes of
controller files.

## Products

| Product | Bytes | Content |
| --- | ---: | --- |
| Raw final H5AD | 5,818,433,232 | published experiment H5AD |
| Filtered final H5AD | 4,737,039,769 | 588,498 cells × 12,926 genes |
| Merged DEA Parquet | 347,041,728 | 30,078,802 rows; 2,327 perturbations |
| Merged GSEA Parquet | 9,778,860 | 114,023 rows; 49 terms |

The QC/comparison report was produced in the run's configured comparison-results
output directory.

Machine-readable results: [jurkat-benchmark.json](jurkat-benchmark.json).

No comparison with the five production datasets has been run.
