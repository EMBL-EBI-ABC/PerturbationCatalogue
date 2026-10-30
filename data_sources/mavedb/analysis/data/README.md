# Homology identifier mapping

`mavedb_uniprot_to_ensembl.tsv` is the 538-row mapping snapshot from
`ot_perturbation_catalog_mave_experiments/meta/mavedb_uniprot_to_ensembl.tsv`,
copied unchanged to reproduce that project's reference selection. `From` is a
UniProt accession and `To` is an Ensembl gene identifier, sometimes versioned.

The reference generator removes version suffixes and duplicate pairs in its
output copy. This snapshot does not cover every protein in newer MaveDB dumps.
The generator reports metadata proteins absent from it; pass `--id-map` with an
updated mapping to extend coverage. No new mappings are inferred from gene names.
