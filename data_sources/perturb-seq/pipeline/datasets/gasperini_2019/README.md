# Gasperini 2019 inputs

The 32 `at_scale_screen` GEX/guide pairs are submitted 10x Genomics BAMs and
are generated from the ENA run report. Gasperini uses the v2 layout: the
captured 20-base guide starts at offset 23 in the long guide read.

`generate_inputs.py` retrieves the published guide dictionary and writes the
feature table. It keeps gene-targeting TSS guides and negative/technical
controls, which are the representations supported by the current gene-level
DEA/GSEA schema; enhancer-only guide targets are not assigned an Ensembl gene
by this pipeline.
