# Replogle K562 GW unified-pipeline benchmark

Clean full run on a Slurm cluster on 19 September 2026 using main-repository commit `e1b2c84`,
Nextflow 25.04.6 and the pipeline-wide image `perturb_seq.sif` (SHA-256
`b72e25f35aeeef3da9841103004dd54fadb349ff44cf12149600d28996e73e07`). This
replaces the earlier resumed benchmark.

## Run

| Item | Value |
| --- | --- |
| Controller | SLURM `51765760`, run `39a3aedc089e4fd4a82c6c7bc8992dbc` |
| Profiler | SLURM `51765769`, run `0de1f5e3f6bd4a338547cc825c54f99c` |
| Start | `2026-09-19T02:59:43.640905+00:00` |
| Finish | `2026-09-19T05:34:22.474926+00:00` |
| Wall time | **9,278.834 s (2h 34m 38.834 s)** |
| Input | 96 sample groups, 4,591 included SRA accessions |
| Excluded input | `SRR19331204` (persistent NCBI S3 authorization failure) |
| Successful / failed attempts | 388 / 1 expected retry |
| Peak sampled footprint | **8,623,764,356,398 bytes (8.624 TB / 7.843 TiB)** |
| Largest task MaxRSS | **63,570,276K (60.6 GiB)**, successful `MERGE_RESULTS` retry |
| Final published output | 395,744,032,507 bytes |

The first `MERGE_RESULTS` attempt exited 137 and Nextflow automatically retried
it. The retry succeeded; this is the expected memory-growth behaviour for this
dataset, not an unrecovered pipeline failure. The profiler sampled every five
seconds and exited successfully.

## Products

The final filtered H5AD contains 1,314,345 cells and 12,420 genes. Its analysis
manifest records 29,492 control cells, 614,397 single-gene perturbation cells,
and 9,548 perturbations after the minimum-cell filter.

| Product | Bytes |
| --- | ---: |
| `experiment_final.h5ad` | 33,962,377,846 |
| `experiment_final.filtered.h5ad` | 10,227,397,126 |
| Merged DEA Parquet | 1,274,411,680 |
| Merged GSEA Parquet | 35,465,328 |
| Summary JSON | 78,447 |

The merged products contain 118,586,160 DEA rows for 9,548 perturbations and
458,304 GSEA rows for 48 terms. QC/comparison reports were produced in the
run's configured comparison-results output directory.

No new comparison against the five production datasets was part of this clean
benchmark run.
