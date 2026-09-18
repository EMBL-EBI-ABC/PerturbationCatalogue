# Perturb-seq validation

`compare_production.py` compares a bounded sample of the unified Replogle
products with the live Perturbation Catalogue production API. It reads the
cluster Parquet products with PyArrow, samples deterministically, queries DEA
through `/v1/perturb-seq/{dataset_id}/search` and GSEA through
`/v1/perturb-seq-gsea`, then reports key coverage, correlations, absolute
errors, sign agreement and GSEA leading-edge overlap.

The latest result is recorded in
[replogle-production-comparison.md](replogle-production-comparison.md).
