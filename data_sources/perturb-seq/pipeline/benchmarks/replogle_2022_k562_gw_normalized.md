# `replogle_2022_k562_gw_normalized` benchmark

Large-dataset end-to-end validation completed on a Slurm cluster, 18 September 2026,
through backed comparison/QC/probe calling, DEA and GSEA. This was a resumed
completion: index and KITE stages were cached, while the standard count wave,
aggregate H5AD and all downstream products completed in the final run.

| Measurement | Result |
|---|---:|
| Main-repository commit | `57061b9548315e1996857336a40117b220bacc85` |
| Samples / SRA accessions | 96 / 4,591 |
| Excluded accession | `SRR19331204` (persistent NCBI S3 authorization failure) |
| End-to-end wall time | **8,539.851 s (2h 22m 19.851s)** |
| Successful / cached / failed / retried tasks | **388 / 98 / 1 / 1** |
| Nextflow peak tasks / CPUs / memory | 191 / 1,536 / 11.9 TB |
| Maximum sampled footprint | **12.073 TB (10.980 TiB)** |
| Largest observed task RSS | 63,522,880K (about 60.6 GiB) |
| Final observed output footprint | 395.744 GB |
| Raw H5AD | 1,314,345 cells; 33.962 GB |
| Filtered H5AD | 1,314,345 cells × 12,420 genes; 10.227 GB |
| Merged DEA | 118,586,160 rows; 9,548 perturbations; 1.274 GB |
| Merged GSEA | 458,304 rows; 48 terms; 35.465 MB |
| DEA/GSEA batches | 191 |

The one failed task was the expected first `MERGE_RESULTS` attempt (exit 137,
22,626,764K maximum RSS); its automatic memory-growth retry succeeded
(63,522,880K maximum RSS). It is recorded as expected dataset-specific
behaviour, not an unrecovered pipeline failure.

The run used the single image `perturb_seq.sif` (SHA-256
`8165d68e8f2f6b26ef5421879ccefb88de3206a3a32745c1105edf32edda935b`). Disk was
sampled every five seconds over work files, published outputs and controller
files. Finished products were published under the configured pipeline
`--outdir`.
