# `nadig_2025_jurkat` benchmark

Clean end-to-end run on Codon, 17 September 2026, from SRA retrieval through
merged DEA/GSEA. The run used main-repository commit `f20badf`, Nextflow
25.04.6, Ensembl release 115 and the single image
`perturb_seq.sif` (SHA-256
`8165d68e8f2f6b26ef5421879ccefb88de3206a3a32745c1105edf32edda935b`).

| Measurement | Result |
|---|---:|
| Samples / SRA accessions | 56 / 1,792 |
| End-to-end wall time | **5,124.989 s (1h 25m 24.989s)** |
| Successful / failed / retried tasks | **222 / 0 / 0** |
| Nextflow peak tasks / CPUs / memory | 57 / 904 / 2.9 TB |
| Maximum sampled footprint | **1.459 TB (1.327 TiB)** |
| Final published output footprint | 72.499 GB |
| Raw H5AD | 739,766 cells × 78,899 genes; 5,341 guides; 5.818 GB |
| Filtered H5AD | 588,498 cells × 12,926 genes; 4.737 GB |
| Merged DEA | 30,078,802 rows; 2,327 perturbations; 347.048 MB |
| Merged GSEA | 114,023 rows; 49 terms; 9.777 MB |
| DEA/GSEA batches | 47 |

Wall time includes fresh reference/index construction, all downloads, extraction,
counting, sample merges, H5AD compression, backed QC/probe calling, input
preparation, all DEA/GSEA batches and final publication. Disk was sampled every
five seconds over work files, published outputs and controller files.

Finished output root:
`/hps/nobackup/mfreeberg/perturb_seq_fastq/results/nadig_2025_jurkat/`.
