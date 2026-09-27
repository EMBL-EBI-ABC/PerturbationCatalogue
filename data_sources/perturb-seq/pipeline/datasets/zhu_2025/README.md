# Zhu 2025

The public raw archive contains per-lane Cell Ranger filtered feature matrices
rather than a complete set of sequencing reads. Eight of the twelve
donor/condition outputs have missing public guide-read lanes, and the libraries
use 10x Flex with Ultima sequencing, which the FASTQ counter does not support.

The unified pipeline therefore imports the complete official GEO Cell Ranger
matrices from `GSE314342_RAW.tar`, preserving raw gene and guide UMI counts.
`ProbeNTC-*` features are excluded, matching the authors' guide-assignment
code. The remaining 26,504 guides use the authors' curated target mapping in
`guide_targets.tsv`. Each lane is converted separately and the lane H5ADs are
concatenated on disk before the common compression, QC/comparison,
probe-calling, DEA and GSEA stages.

Use `--cellranger_h5_tar`, `--cellranger_sample` (for example `D1_Rest`) and
`--guide_targets datasets/zhu_2025/guide_targets.tsv` instead of a sequencing
sample sheet and count references. The author-provided `*.assigned_guide.h5ad`
for the same donor/condition remains the independent curated comparison input.
