# Zhu 2025 read-level inputs

These twelve donor/condition datasets use multiplexed 10x Flex libraries with
Ultima sequencing. Gene-expression inputs are the original paired FASTQs;
guide inputs require full-quality SRA extraction including technical reads.
Sources, physical lane grouping and donor/state probe-barcode assignments are
recorded in [the raw-data notes](../../../raw-data-notes.md).

Four outputs in the R1 pool have complete archived guide-lane coverage. Eight
outputs in the R2 pools have missing guide lanes; process their available reads
and retain this limitation in their provenance. The earlier Cell Ranger
matrix-import route has been removed. Read-level processing requires a
probe-aware counter and separate gene-expression/guide barcode demultiplexing.

`guide_targets.tsv` supplies the curated guide target metadata. Author
`*.assigned_guide.h5ad` files are independent comparison references.
