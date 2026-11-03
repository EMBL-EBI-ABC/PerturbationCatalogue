# MaveDB GMM validation record

## Acceptance status

The Python implementation has a working, bounded local smoke workflow for
both model families. Deterministic post-processing and the draw/LOO adapters
have shared-input R references. This establishes the local implementation
milestone; it does not establish production scientific equivalence or readiness
to replace the production R run.

## Validated environment

The measured runs used Python 3.12.8, NumPy 2.2.6, pandas 2.2.3, SciPy 1.14.1,
xarray 2025.3.0, ArviZ 0.22.0, CmdStanPy 1.2.5, stanio 0.5.1 and CmdStan 2.37.0.
The host was an arm64 Mac with 16 GiB of physical memory and Apple clang++
21.0.0 (`arm64-apple-darwin25.6.0`). The shared R references were generated
with R 4.6.1, CmdStanR 0.9.0, posterior 1.7.0, loo 2.10.1, sn 2.1.3 and
readr 2.2.0. See the [fixture provenance](../tests/fixtures/gmm/README.md) and
[Stan source hashes](stan/README.md).

CmdStan is installed explicitly under `$HOME/.cmdstan/cmdstan-2.37.0`. The fitting
CLI does not install packages, download CmdStan or compile unrelated candidate
models.

## Deterministic and shared-input parity

The checked-in [R-generated fixtures](../tests/fixtures/gmm/reference/) cover
the five candidate parameterizations (`flexible`, `ko`, `del`, `neut`, `pos`),
transforms, output densities and classifications, posterior summaries, chain-
aware relative efficiency, pointwise PSIS-LOO and modes comparison.

- Transform, density and ordinary summary values pass comparisons at
  `rtol=1e-8`, `atol=1e-10`. The checks retain the R output-density constants,
  unclipped post-fit logistic behavior, strict classification order, repeated
  observations, and no density Jacobian.
- On identical likelihood arrays, the 500-observation PSIS reference differs
  by at most `1.3e-14` for ELPD and `9.7e-14` for Pareto-k. The fixed-tail
  fixture also covers repeated values at the R PSIS tail cutoff.
- On the same binary draws, scalar R-hat differs from R by about `3e-15`.
  Decimals parsed from CmdStan CSVs can round a draw by one ULP; because R-hat
  is rank-based, tied ranks can then shift by `0.005`. The shared binary-array
  checks pass; raw CSV decimal text is not claimed to produce bit-identical
  rank diagnostics.
- Modes comparison is tested in source order (`ko`, `del`, `neut`), with
  `pos` retained in outputs but excluded from winner selection and
  `p_loo < 12` applied after comparison without rebasing differences.
- Both model families now have fixed-draw R goldens for all six post-processing
  tables in `../tests/fixtures/gmm/tables/`. The injected reference runner uses
  the shared R draws, summaries and pointwise LOO values; the tests check exact
  headers, candidate/row order, identities and labels, plus numeric values at
  `rtol=1e-8`, `atol=1e-10`. This also exercises the legacy parameter metadata
  order and neutral posterior fraction omission through the complete
  `fit_study` output path.
- The clean R reader harness passed against real flexible and modes smoke
  outputs and an all-failed modes placeholder fixture. The contract checks also
  cover mixed success/failure outputs and the literal identifier `NA`.
- The optional R-backed pytest checks in `../tests/test_gmm_readers.py` cover the
  full-table goldens, a filtered-best modes case with a retained negative ELPD
  difference, mixed and all-failed studies, and a flexible `p_loo > 12` case
  that must not use the modes-only cutoff. They skip when Rscript/readr is not
  available.

The R table-reader checks can be repeated after generating the smoke outputs:

```bash
Rscript --vanilla data_sources/mavedb/tests/read_gmm_outputs.R \
  data_sources/mavedb/output/gmm-validation/smoke/maves/flexible flexible 0
Rscript --vanilla data_sources/mavedb/tests/read_gmm_outputs.R \
  data_sources/mavedb/output/gmm-validation/smoke/maves/modes modes 0
Rscript --vanilla data_sources/mavedb/tests/read_gmm_outputs.R \
  data_sources/mavedb/tests/fixtures/gmm/reader_all_failed_modes modes 0
```

The recorded outcomes were one successful flexible study, one successful modes
study, and one failed study with valid sentinel tables.

## Real-sampler smoke runs

The smoke fixture is one fixed 500-row reservoir sample from the 3,180 eligible
rows in `urn:mavedb:00000001-a-1`. Both runs used one study, 500 selected and
fitted rows, two sequential chains, 200 warmup draws and 200 retained draws per
chain, seed `123`, and the 2,048 MiB memory budget. Models were already
compiled; timings below exclude initial toolchain/model compilation.

| Model | Elapsed | Observed aggregate peak RSS | Observed temporary fit scratch | Output rows: score / density / params / posterior / summary / LOO |
| --- | ---: | ---: | ---: | ---: |
| `flexible` | 10.16 s | 291,733,504 B (278.2 MiB) | 6,317,540 B (6.03 MiB) | 500 / 4,920 / 1 / 400 / 15 / 1 |
| `modes` | 21.56 s | 296,960,000 B (283.2 MiB) | 6,328,226 B (6.04 MiB) | 2,000 / 19,680 / 4 / 1,600 / 48 / 3 |

Aggregate RSS was sampled across the Python process and its CmdStan children at
50 ms intervals. Temporary scratch is the peak logical size of the per-fit
temporary directory; it excludes compiled model/cache files, library caches
and final/staged TSV output. These are observed values, not OS-enforced limits.
The six result tables and `qc/batch0/` were written for each successful run.

The flexible smoke fit reported 75 divergent transitions and Pareto-k above
0.7 for 3 of 500 observations (maximum `1.724734`). These warnings remain in
QC. Short chains and a small fixture validate execution, schemas and resource
use; they are not scientific acceptance criteria or evidence of convergence.

For bounded LOO storage, a separate probe processed a 512,000,000-byte logical
likelihood array (2 chains × 200 draws × 160,000 observations) in 62 blocks.
The largest copied block was 8,387,200 bytes; the measured process peak RSS was
141,737,984 bytes (135.2 MiB), with 3.88 seconds elapsed. This verifies the
disk-backed block reader and page eviction at that synthetic scale. It is not a
full CmdStan fit or a representative production-study benchmark.

Initial environment/model build costs were measured separately: the CmdStan
`make -j1` build after download and extraction took 70.94 seconds with
1,603,174,400 bytes peak aggregate RSS; serial compilation of all five
candidate models took 25.74 seconds with 896,696,320 bytes peak aggregate RSS.
The CLI compiles only models required by the selected family and reuses
matching cached executables.

## Test evidence

The focused GMM suite is run with:

```bash
.venv/bin/python -m pytest -q data_sources/mavedb/tests/test_gmm*.py
```

It covers deterministic R fixtures, expected numerical failures, injected
candidate runners, repeated rows, mode rejection and late-candidate failure,
family schemas, canonical/legacy input contracts, source identities,
zero-eligible studies, preflight rejection, output staging/rollback, streamed
batch tables, and CmdStan CSV/LOO adapters. The focused GMM suite has passed
**62 tests**. Warning propagation and seed handling regression coverage did
not alter the fixed-draw outputs.

## Remaining production gates

- A representative full-size study has not been sampled end to end to measure
  peak RAM, scratch disk use and runtime. The 160,000-observation block-reader
  probe does not close this gate.
- Independent Python-versus-R sampling comparisons are not complete. Fixed
  seed integers do not imply the same random samples or draws; compare a shared
  explicit row subset and assess posterior differences against Monte Carlo
  uncertainty before changing production commands.
- No production workflow cutover or Nextflow/Slurm integration has been made.

Keep the R workflow as the production reference until these gates are met.
