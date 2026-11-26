# Bounded GMM input

These canonical scores and metadata are the bundled input for the standalone
`--profile smoke` and Nextflow GMM smoke runs.

The source score set is `urn:mavedb:00000001-a-1` (`UBE2I imputed & refined`),
a human protein-coding MAVE. Its original publication attribution is retained
in `metadata.tsv`. MaveDB source data is CC0, as documented in the MaveDB README.

The 500 rows are a fixed reservoir sample selected with `random.Random(123)`
from 3,180 eligible original rows, sorted by `source_row`. Original `dataset_id`
and `source_row` values are preserved. The full source metadata reports 3,180
variants; the fixture deliberately contains only 500. This is already processed
GMM input, and must not be passed through raw-data processing again.

Use `--processed_dir` with `--gmm_profile smoke` to bound fitting to one study,
500 rows, two sequential chains, and 200 warmup/200 retained draws per chain.
Use the fixture to check execution and output schemas. Assess scientific
performance on representative score sets with production sampler settings.
