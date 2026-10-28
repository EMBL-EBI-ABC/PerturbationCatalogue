# Python implementation plan for MaveDB GMM batch fitting

Refined on 2026-10-28 against `ot_perturbation_catalog_mave_experiments` commit `da4b24f2f2ad6a8d86d10831e34699199c748300` (clean working tree). This document plans the port; it does not implement or run the models.

**Local development requirement:** the first usable implementation must run the real GMM workflow on a tiny subset on a Mac with 16 GB RAM. Provide an explicit `smoke` profile with small input fixtures, bounded rows/draws, sequential execution and separate outputs. A full-size production benchmark must not block delivery of this local workflow. Section 7 defines the settings and acceptance checks.

## 1. Scope and implementation location

Replace the R orchestration and post-processing used by [fit_gmm_batch.R](/Users/zakirov/Documents/GitHub/ot_perturbation_catalog_mave_experiments/bin/fit_gmm_batch.R:1) and the fitting helpers in [gmms.R](/Users/zakirov/Documents/GitHub/ot_perturbation_catalog_mave_experiments/src/gmms.R:1). Keep the five existing Stan models and their shared `functions.stan` as the statistical definitions, including their priors, constraints, transformed parameters and generated quantities.

Implement the migrated code **in PerturbationCatalogue**, following the ownership decision in [the pipeline plan](mavedb-pipeline-implementation-plan.md). The draft's `src/mavenn/gmm.py` and `pyproject.toml` refer to the sibling analysis repository; this repository has neither that package nor that project file. Use this proposed layout instead:

```text
data_sources/mavedb/
  scripts/fit_gmm_batch.py       # CLI and path resolution
  analysis/
    gmm.py                      # public workflow, transforms and compatibility tables
    _gmm_stan.py                 # CmdStan execution and draw/LOO adaptation
    stan/
      functions.stan
      fixed_gmm_flexible.stan
      fixed_gmm_ko.stan
      fixed_gmm_del.stan
      fixed_gmm_neut.stan
      fixed_gmm_pos.stan
  requirements-gmm.txt
  tests/                        # fit existing test conventions during implementation
    fixtures/gmm/
      smoke/                    # tiny catalogue inputs for the local smoke profile
```

Copy the Stan sources with their include relationship intact and preserve attribution/licences. Resolve model files from the installed/source module location, independently of the working directory. Production execution must not read a sibling checkout. CmdStanPy can compile and run the retained Stan definitions; the production Python path needs no R runtime. [CmdStanPy overview](https://mc-stan.org/cmdstanpy/users-guide/overview.html)

Initial scope is `maves` with `flexible` and `modes`. Keep score-set normalization distinct from cross-study meta-analysis. Porting plotting functions, `make_gmm_singles`, `make_gmm_wide`, training-data filters, new model families, publication to databases, and Nextflow/Slurm integration are subsequent work. The CLI should be callable by a future Nextflow task without scheduler-specific code.

## 2. Small module interface

Expose a typed interface with these responsibilities:

```python
plan_batches(study_counts, variants_per_batch) -> BatchPlan
fit_study(study_id, rows, model, settings, runner) -> StudyFit
run_batch(config) -> BatchSummary
```

- `BatchConfig` owns validated input/output paths, input format, model family, execution profile, study/row selection, batch selection, sampler settings, memory budget and optional seed.
- `fit_study` operates on one score set, shares one fitting subset across candidate models, selects the winner and produces the six logical result tables plus analysis QC. Keep transforms and classification callable without CmdStan.
- The internal runner accepts raw scores, candidate name and sampler settings. It returns compact parameter draws, generated scale values, compatible summaries and pointwise LOO results. Keep CmdStan objects and temporary CSV handling inside the adapter.
- `run_batch` selects eligible studies, handles expected study failures and writes outputs incrementally. `BatchSummary` contains counts, failures and paths, rather than all study fits.

The six logical outputs need not all be materialized DataFrames: retain compact fit tables, and generate large per-variant results in chunks after choosing the winner. An injectable runner permits workflow tests using fixed draws or controlled failures without invoking Stan. Start with these two implementation modules; split further only when responsibilities warrant it.

## 3. Input and CLI contract

### Inputs

Support explicit `catalogue` and `legacy` input adapters. Do not guess the format or infer a score column from its position.

| Adapter | Inputs | Required identity and selection |
| --- | --- | --- |
| `catalogue` (default) | Current processor's `scores.tsv` and `metadata.tsv` | Join on full `dataset_id`; preserve `(dataset_id, source_row)`; read `meta_study` and `analysis_types`; validate the pipeline's single-human-target eligibility from source metadata |
| `legacy` | Original `data/maves/full.tsv` and `data/maves/meta.tsv`, supplied explicitly | Join on `score_urn`; metadata requires `meta_study`; score rows require `score_urn`, `hgvs_nt`, `hgvs_pro`, `score`; retain original row order and multiplicity |

The [current processor contract](../data_sources/mavedb/README.md) retains all original variant records in the canonical scores table and removes obsolete category-specific `full.tsv`/`meta.tsv` outputs. The Python fitter must therefore consume the new tables directly; it must not require regenerating the retired tables.

For catalogue input, select datasets whose `analysis_types` includes `maves`, exclude meta studies, and validate the pipeline's single-human-target requirement. Do not impose that new cohort restriction on legacy parity fixtures: cohort selection and fitting parity are separate checks. Reject ambiguous/multi-target catalogue metadata with a recorded eligibility reason rather than pooling targets.

The fitting and transformation subset consists of finite numeric `score` values with nonmissing protein HGVS in the MAVE dataset. This reproduces the original upstream MAVE subset; it includes WT, synonymous, nonsense and multiple-mutant records. Do not restrict to simple substitutions or aggregate repeated scores/HGVS. Preserve excluded canonical rows in the source tables and emit identity-keyed exclusion QC with reasons. Every eligible row receives a transformed result, including rows outside the fitting sample. The smoke profile deliberately restricts both fitting and transformation to its selected tiny subset, with the restriction recorded in QC. Missing nucleotide HGVS is allowed.

Validate required headers, metadata key uniqueness and catalogue observation-key uniqueness before fitting. Define explicit handling for unknown metadata joins and missing eligibility fields; structural ambiguity is an input error. Preserve literal identifiers such as `NA` in text columns by configuring TSV missing-value parsing explicitly. An empty fitting subset or zero/nonfinite sampled MAD is an expected study failure.

### CLI

Keep the existing option names and defaults:

| Option | Contract |
| --- | --- |
| `--profile` | `production` (default) or explicit `smoke`; resolve and print all effective settings before fitting |
| `--study_type` | Default `maves`; reject MPRA/splicing before loading or compiling |
| `--model` | `flexible` or `modes`; default `modes` |
| `--index` | Integer; `0` processes all eligible score sets; positive values select a planned batch; reject negatives |
| `--variants_per_batch` | Positive integer; default `50000`; batching target, not a fitting cap |
| `--out` | Default `intermediates`; production results go below `<out>/maves/<model>/`, smoke results below `<out>/smoke/maves/<model>/` |

Add explicit `--scores`, `--meta`, `--input-format`, `--seed`, and `--list-batches` options. Production catalogue input paths default to `data_sources/mavedb/output/processed/{scores,metadata}.tsv`, resolved from the repository/module location like the current processor. Smoke input defaults to the small fixture directory in section 7, avoiding a scan of the full dataset. Explicit input paths override either default. Legacy input requires explicit paths. `--list-batches` validates/counts inputs and prints assignments and batch count without compiling or sampling.

Expose `--study-id` (repeatable), `--max-studies`, `--max-study-rows`, `--max-fit-rows`, `--chains`, `--parallel-chains`, `--iter-warmup`, `--iter-sampling`, and `--memory-budget-mb`. Keep these in the typed settings object. `--max-study-rows` caps the actual selected data used by fitting and post-processing; `--max-fit-rows` only limits the fit sample. Keep production defaults separate from smoke settings, and record every explicit override. Reject invalid settings before creating result files.

Example intended production invocation from this repository root; use the smoke commands in section 7 during development on the 16 GB Mac:

```bash
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --input-format catalogue --study_type maves --model modes \
  --index 1 --variants_per_batch 50000 --seed 123 --out intermediates
```

## 4. Numerical behavior to preserve

These details are part of the compatibility contract, even where the source implementation may be internally inconsistent. **The Python parity port must preserve the behaviors below as written. Diagnosing or correcting them is out of scope; consider any change separately with scientific review and new reference outputs.**

### Potential quirks and bugs explicitly out of scope

- **Skew-normal locations differ:** Stan's likelihood uses deleterious/positive locations `0.75`/`1.25`; R's post-fit density and `mode_density` helpers use `0.8`/`1.2`. Keep the two definitions in their current roles.
- **Neutral noise weight differs:** the neutral Stan candidate uses `theta_noise=0.03`, but R's output-density and density-classification helpers default to `0.01`. Preserve both; do not pass the fitted weight through in this port.
- **Positive cannot win:** the `pos` model is fit and exported, but excluded from the modes LOO comparison. Keep it in outputs with `top_model=FALSE` and keep it out of `loo_tbl`.
- **LOO filter follows comparison:** compare first, then retain rows with `p_loo < 12`. Keep `elpd_diff` and `se_diff` relative to the original comparison winner, even if that winner is filtered out.
- **KO threshold adds the location twice:** `classify_point_sd` uses `0.1 + qnorm(0.99, 0.1, 0.05)`, giving a threshold of about `0.316` rather than the single-location 99th percentile of about `0.216`. Preserve the threshold and assignment order.
- **Transform clipping differs:** Stan clips the logistic argument to `[-20, 20]` in the likelihood; R's post-fit score transform does not. Preserve this distinction for extreme scores.
- **Transformed density omits a Jacobian:** the R helper maps the `w` density onto the `s` axis without changing the density by a Jacobian. Preserve the generated `s`-view values.
- **Density scaling is family-specific:** flexible density outputs divide by `100`; modes outputs divide by each candidate's maximum. Keep these display values as they are.

These are compatibility requirements for this port, not endorsements that each behavior is statistically preferable. Record them in regression fixtures so a later, explicitly scoped correction can identify and review its effects.

### Sampling and scale

- Production batch defaults are four chains, four parallel chains, 1,000 warmup iterations, 2,000 retained iterations per chain, `adapt_delta=0.9`, `max_treedepth=10`, and `refresh=4000`. The helpers' 1,000-sample defaults are overridden by the batch script. The explicit smoke profile uses the reduced settings in section 7.
- In production, fit at most **256,000 eligible rows**, sampled without replacement. All four modes use exactly the same sampled rows and ordering. Apply fitted parameters to all eligible study rows, not just sampled rows. Smoke applies this same workflow inside its bounded data subset.
- Pass raw `y` and `N` to Stan. Each model computes `median(y)` and `median(abs(y - median(y))) * 1.4826` on the sampled data. Do not standardize twice, recompute scale on the full study, or substitute a slightly different normal-consistency constant.
- Post-processing uses posterior **medians** of the selected scalar parameters, including generated `y_median`/`y_mad` renamed to `median`/`mad`.
- Preserve candidate initial values from the R helpers:

| Candidate | Initial values |
| --- | --- |
| `flexible` | `sigma_neut=0.1`, `beta=1.5`, `delta=0`, `frac_neut=0.8`, `frac_non_neut=[0.8, 0.19, 0.01]` |
| `ko` | `sigma_neut=0.1`, `beta=2`, `delta=-0.4`, `frac_neut=0.65` |
| `del` | `sigma_neut=0.1`, `beta=1.25`, `delta=-0.1`, `frac_neut=0.7` |
| `neut` | `frac_neut=0.9`, `beta=1.01`, `delta=0`; `sigma_neut` is fixed in Stan |
| `pos` | `sigma_neut=0.1`, `beta=1.25`, `delta=0`, `frac_neut=0.75` |

Keep the R production default of an unspecified seed when `--seed` is absent; smoke defaults to seed `123`. When supplied, seed both subsampling and Stan, using reproducible study/candidate seeds independent of batch membership and Python's randomized `hash()`. Matching R and NumPy seed integers does not reproduce their sampled rows: parity fixtures must supply the actual chosen row indices.

### Transforms, density and classifications

Port `g`, `g_inv`, `f`, `f_inv`, `dflexgmm`, `transformed_density`, `transformed_scores`, `classify_point_density`, and `classify_point_sd` from [the R helpers](/Users/zakirov/Documents/GitHub/ot_perturbation_catalog_mave_experiments/src/gmms.R:1).

- Standardize `s` with the sampled median/MAD, then compute `w = 2 * expit(beta * (g(s) + delta))` using a stable logistic implementation. R post-processing is unclipped; the Stan likelihood clips its logistic argument to `[-20, 20]`. Preserve that distinction, including extreme-score behavior.
- R density helpers use deleterious/positive skew-normal locations **0.8/1.2**; `functions.stan` uses **0.75/1.25**. Keep the former for output densities/labels and the latter for likelihood/LOO. Do not reconstruct likelihoods with the output-density helper.
- Output-density/classification helpers default to `theta_noise=0.01`, even for the neutral candidate whose Stan `theta_noise` is `0.03`. They pass fitted `sigma_neut` and the four fitted component weights. Preserve the helper behavior; do not renormalize component weights or substitute fitted noise silently.
- `mode_density` selects the maximum weighted component density, with first-tie order `ko`, `del`, `neut`, `pos`, `noise`. Repeated values must retain separate rows. Verify the R `sn::dsn` and SciPy skew-normal parameter conventions with fixed examples.
- The scoring caller uses `mode_sd` threshold **0.99**, overriding the helper's 0.95 default. Preserve strict inequalities and assignment order: start neutral, assign positive, deleterious, then KO. The source's KO boundary is `0.1 + qnorm(0.99, 0.1, 0.05)`, including both `0.1` terms.
- Density grids contain 500 evenly spaced `w` values from `0.001` through `2`. The `s` view transforms only points satisfying `0.03 < w < 1.97`, keeping their density ordinates. It does **not** apply a change-of-variable Jacobian. Flexible outputs use `scaled_density = density / 100`; modes outputs divide by the maximum density within each candidate.
- Retain the Stan likelihood and generated `log_lik` exactly as defined. Adding Jacobian corrections or changing the positive candidate's selection eligibility would change the analysis and requires a separate scientific decision.

## 5. Draw summaries and PSIS-LOO: prove parity early

Add a dedicated `requirements-gmm.txt` using NumPy, pandas, SciPy, CmdStanPy and a pinned, tested ArviZ statistics stack. Follow the repository's requirements-file convention; do not edit the sibling project's `pyproject.toml`. Record compatible Python/CmdStan/toolchain versions in the analysis setup instructions. Use this repository's `.venv` for local Python work. Provision CmdStan and the C++ toolchain explicitly; do not download/install them implicitly during batch fitting.

### Draw and summary adaptation

Preserve chain structure until diagnostics and LOO are complete. Normalize posterior identifiers to R's one-based `iteration`, `chain`, `draw` convention using a fixture that exposes flattening order. Export retained draws only. Rename `frac_non_neut[1:3]` to `frac_ko`, `frac_del`, `frac_pos` in the relevant outputs, and `y_median`/`y_mad` to `median`/`mad` in parameter summaries.

The R summary is produced by `posterior::summarise_draws`, rather than an arbitrary CmdStan summary. Reproduce `mean`, `median`, sample `sd`, scaled `mad`, `q5`, `q95`, rank-based `rhat`, `ess_bulk`, and `ess_tail`, including constant-parameter missing diagnostics. Compute only required scalar summaries; summarizing every `log_lik[n]` wastes substantial work. A column rename alone is insufficient. [CmdStanR summary contract](https://mc-stan.org/cmdstanr/reference/fit-method-summary.html)

The neutral R posterior export omits `frac_neut`, although the parameter exists and its fit summary includes it. Capture this omission in the legacy contract and keep that posterior cell missing until a deliberate schema correction is made. Optional fraction columns differ between model families; do not invent values for absent exports.

### LOO compatibility adapter

Treat LOO compatibility as an early implementation gate, before running large studies. Use the sampled observations' generated `log_lik`, aligned identically across candidates. Retain pointwise results and Pareto-k diagnostics. Do not enable model averaging, moment matching or a different LOO approximation by default.

CmdStanR's default `$loo()` computes **per-observation relative efficiency** for likelihood draws. Reproduce that chain-aware calculation; do not substitute `reff=1` or a scalar average of parameter ESS. A numerically shifted exponentiation may be used to avoid overflow without changing relative efficiency. [CmdStanR LOO defaults](https://mc-stan.org/cmdstanr/reference/fit-method-loo.html), [likelihood relative efficiency](https://mc-stan.org/loo/reference/relative_eff.html)

ArviZ's documented high-level `loo` accepts a scalar `reff` and does not directly supply every historical R output column. Use its PSIS primitives through a small compatibility adapter that handles each observation's efficiency, rather than assuming a high-level call matches R. Pin the chosen API and demonstrate its behavior on shared arrays. If necessary, process observations individually or in bounded blocks. [ArviZ statistics LOO API](https://python.arviz.org/projects/stats/en/v1.3.0/api/generated/arviz_stats.loo.html)

Construct the existing LOO columns from pointwise quantities, checking formulas against the pinned R reference:

```text
lpd_i       = logmeanexp(log_lik[:, i])
p_loo_i     = lpd_i - elpd_loo_i
elpd_loo    = sum(elpd_loo_i)
p_loo       = sum(p_loo_i)
looic       = -2 * elpd_loo
se_elpd_loo = sqrt(N * sample_variance(elpd_loo_i))
se_p_loo    = sqrt(N * sample_variance(p_loo_i))
se_looic    = 2 * se_elpd_loo
```

Handle degenerate/insufficient-data cases explicitly. These standard errors measure variation across observations, not sampler MCSE. [R loo output definitions](https://mc-stan.org/loo/reference/loo.html)

For `modes`, preserve the exact sequence in [process_study_modes](/Users/zakirov/Documents/GitHub/ot_perturbation_catalog_mave_experiments/src/gmms.R:506):

1. Fit `ko`, `del`, `neut`, and `pos`, including their helper-level LOO computation.
2. Compare only `ko`, `del`, `neut`, sorting by decreasing ELPD. Calculate `elpd_diff = elpd(model) - elpd(best)` and paired `se_diff` from pointwise differences. Use deterministic source-order tie handling.
3. Apply **`p_loo < 12` after comparison**, without rebasing the surviving difference columns. The selected model can therefore have a negative `elpd_diff` if the original comparison winner was filtered out.
4. Select the first remaining model. If none remain, return a failed study.
5. Retain all four successful candidates in score, density, parameter, posterior and summary outputs. Only the selected candidate has `top_model=TRUE`; `pos` is never eligible to win. The LOO table contains only retained comparison rows.

The flexible family has one candidate, always marked top on success; its `elpd_diff` and `se_diff` are missing and it has no modes-specific `p_loo` filter. Pairwise standard errors must use aligned observations, not differences between marginal standard errors. [R loo comparison](https://mc-stan.org/loo/reference/loo_compare.html)

## 6. Stable outputs and expected failures

Write these six historical files below `<out>/maves/<model>/` for production and `<out>/smoke/maves/<model>/` for smoke:

| File | Required content |
| --- | --- |
| `score{index}.tsv` | `score_urn`, original `hgvs_nt`/`hgvs_pro`, raw score as `s`, transformed `w`, `mode_density`, `mode_sd`, `model`, `top_model` |
| `density{index}.tsv` | `score_urn`, `model`, `top_model`, component `mode`, view `name` (`s`/`w`), `x`, `density`, `scaled_density` |
| `fitted_params{index}.tsv` | `score_urn`, `model`, `top_model`, posterior-median scalar parameters, including scale, weights and `lp__` |
| `posterior{index}.tsv` | `score_urn`, `model`, `top_model`, draw identifiers, `lp__`, exported fractions, weights, `beta`, `delta`, `sigma_neut` |
| `fit_summary{index}.tsv` | `score_urn`, `model`, `top_model`, `variable`, `mean`, `median`, `sd`, `mad`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail` |
| `loo{index}.tsv` | `score_urn`, `model`, `elpd_diff`, `se_diff`, `elpd_loo`, `se_elpd_loo`, `p_loo`, `se_p_loo`, `looic`, `se_looic` |

Freeze column order and nullable dtypes per family using reference fixtures. Use an explicit, stable schema union so headers do not depend on which studies failed. Preserve source ordering within each candidate, with top candidate first and the R candidate ordering for the others. Write `NA`, `TRUE`/`FALSE`, tab delimiters, a single header, no DataFrame index, adequate numeric precision and round-trippable text quoting.

For catalogue input, `score_urn` is an explicit compatibility alias of the full `dataset_id`. Additionally carry `dataset_id` and `source_row` in score outputs so repeated observations remain traceable. These additive columns are intentional and tested separately from exact legacy headers. Keep canonical source files unchanged. Put study diagnostics and identity-keyed exclusions below a `qc/` subdirectory so legacy filename patterns cannot consume them accidentally. Record input/eligible/fitted/excluded row counts, candidate failure reasons, sampler warnings and Pareto-k summaries; this is analysis QC, not a new release/versioning system. Smoke QC also records the profile, effective settings, selected source-row identities and subsetting status; reduced-data outputs must be identifiable even after copying them elsewhere.

Expected numerical/data failures produce one score placeholder per eligible source row and one sentinel row in each other logical table. Retain identifiers/raw `s`, set fitted values, classifications, `model` and `top_model` to missing, and populate the complete declared headers. Preserve the R all-or-nothing study behavior: if any required candidate's fitting/LOO fails, do not silently select among survivors. Extend safe failure handling to later candidates and post-processing, where the current R implementation can instead abort. Fix `empty_fit`'s omitted columns, including `theta_ko` in posterior placeholders.

Catch known runner and numerical failures at the study boundary; do not swallow programming errors. Missing dependencies, compile failures, malformed global configuration and output I/O errors fail the batch with a nonzero exit. A completed batch may contain failed-study placeholders; report their count explicitly. Convergence/Pareto warnings remain visible QC and do not introduce new selection filters during the parity port.

Stage files in an isolated run directory, validate all six before publishing them, and preserve existing outputs on computation failure. Document that separate filesystem replacements do not form a multi-file transaction; downstream jobs must wait for successful batch completion. Avoid combining old R files and new Python files in one run directory, especially because batch numbering is being corrected.

## 7. Batching, memory and execution

### Local smoke profile for a 16 GB Mac

Make this a supported CLI workflow using the real Stan sampler, LOO, post-processing and six output writers. Ship small catalogue `scores.tsv`/`metadata.tsv` fixtures under `data_sources/mavedb/tests/fixtures/gmm/smoke/`, with at least one study having nonzero MAD and a representative score distribution. Supply attribution for any extracted public data. Generate fixtures once with fixed selection and retain source identities; repeated smoke runs should read these small files directly.

| Setting | Smoke default |
| --- | --- |
| Studies | At most `1`; optional explicit `--study-id` |
| Selected rows per study | At most `500`, for both fitting and transformation |
| Fitting rows | At most `500` within the selected rows |
| Chains | `2`, with `parallel_chains=1` |
| Warmup / retained draws | `200` / `200` per chain; do not save warmup draws |
| Candidate execution | Sequential; `flexible` runs one model, `modes` runs all four unchanged |
| Seed | `123`, overridable explicitly |
| Input read chunks | At most `5000` rows, with only required columns materialized |
| Memory preflight budget | `2048` MiB for sampling/post-processing; reject estimates above this before fitting |
| Results | `<out>/smoke/maves/<model>/`, always separate from production |

Apply study selection and row caps **before** fitting or generating per-variant outputs. With custom full input files, filter by selected study IDs during chunked reads and use a seeded bounded reservoir to select rows without loading a whole study. Default study selection is deterministic; fail clearly if an explicitly requested study is absent or ineligible. Preserve original row identities/order after selection. All candidates must receive the same selected rows. An empty/zero-MAD subset fails explicitly rather than triggering a larger-data fallback. Small-subset exclusions must not be reported as missing or invalid source scores.

`--variants_per_batch` is not a memory control: an oversized study still occupies a batch by itself. Likewise, reducing only `--max-fit-rows` leaves full-study post-processing unchanged. The smoke profile must cap both the number of studies and the data passed to the entire workflow. Default to `--index 0` over the bounded selection; `0` must not bypass smoke limits.

At these defaults, one float64 likelihood array is `2 * 200 * 500 * 8 = 1,600,000` bytes, about **1.53 MiB**, before overhead. Sequential candidates avoid multiplying this by four. Use a conservative preflight estimate that includes parsing, array copies, PSIS workspace and concurrent chains. An estimate is not an OS-enforced limit: validate the implementation by measuring the aggregate peak memory of Python and its CmdStan children, with a sampling/post-processing acceptance target below **2 GiB**. Report compilation memory separately and compile one required model at a time with no parallel builds. Reject overrides that exceed the profile's budget; never fall back to production limits or silently reduce requested settings.

Intended commands **after implementation and explicit environment setup**:

```bash
# Fastest real-sampler check: the bundled subset and one model.
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile smoke --model flexible

# Full modes workflow on the same bounded subset; all four candidates remain.
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile smoke --model modes
```

Compile only the requested family's models and reuse cached executables. First-run compilation may dominate elapsed time; benchmark it separately from repeat runs. Aim for a repeat run measured in minutes, then document actual elapsed time and peak memory on the development Mac rather than promising an unmeasured runtime. Short chains and tiny datasets validate execution, schemas and resource use; their convergence and selected modes are not production scientific results. Keep diagnostics enabled and never relax the LOO selection rules simply to make a smoke run pass.

### Production batching and larger studies

Use deterministic decreasing eligible-row counts with full study ID as the tie breaker, then the source's sequential accumulator strategy. Define positive batches as **1-based** and reserve `0` exclusively for all-studies mode. Oversized studies receive their own batch. Reject a positive index outside the generated range before fitting; an entirely empty selection is a clear no-data error. Verify that every selected study appears once and that indexed-batch unions equal all-studies selection.

This deliberately fixes the R accumulator's initial batch-0 behavior when early studies fit under the target. Recompute job-array ranges from `--list-batches`; do not copy the README's historical `1-94` range or positional `1` example. Keep whole studies together and count eligible rows before applying the fitting cap.

Read metadata and count/select score rows in chunks. Stage or reread selected studies as necessary; do not load all inputs separately for each batch job or assume metadata's total-row counts equal analysis-eligible counts. Hold one study and one candidate's large likelihood state at a time, release it after reducing LOO, and append outputs rather than concatenating a full batch in memory. Default to sequential studies/candidates with four parallel chains in production to avoid accidental nested parallelism. Smoke overrides chain parallelism to one.

At the production cap, `4 * 2000 * 256000` float64 log-likelihood values occupy **16.384 GB (about 15.3 GiB) per candidate**, before CSV text, parser overhead, PSIS work arrays or copies. Four such arrays alone exceed 61 GiB. The original script's large R heap limits are evidence that selected-input reads alone do not solve the memory problem.

Prototype bounded likelihood extraction and PSIS processing early. A streaming conversion of CmdStan CSVs to a disk-backed array, followed by observation-block processing, is a suitable approach if the selected APIs cannot avoid full allocations. Verify actual peak memory: requesting selected variables may still cause the backend to parse/cache all draws. Keep scalar posterior draws and pointwise LOO reductions after releasing raw likelihoods. Preserve the Stan likelihood; memory optimization must not switch to the mismatched R output densities, reduce sampling iterations or thin draws silently.

Compile each required model once per environment/process, with a writable build location and the shared include available. Rebuild when `functions.stan` or compilation settings change; prevent concurrent jobs from racing to build the same executable. Keep sampling files in job-local scratch, retain useful logs on failures, and define scratch retention/cleanup. Measure scratch disk use as well as RAM and runtime on a representative large study. [CmdStanPy compilation and output API](https://mc-stan.org/cmdstanpy/api.html)

## 8. Implementation phases and exit criteria

### Phase 1 — Capture the R contract

Create compact fixtures for both model families: selected input rows, fixed parameter/draw arrays, generated scales, pointwise log likelihoods, six R output tables and failure examples. No existing GMM result directories were present in the inspected checkout, so generating reference outputs is a prerequisite. Preserve the source revision and R/CmdStanR/posterior/loo/CmdStan versions with the fixtures.

Include small/tied samples, extreme scores, repeated HGVS/scores, a sampled study, zero MAD, invalid/missing values and candidates rejected by `p_loo`. Capture mixed successful/failed and all-failed schemas. Freeze the catalogue adapter's identity and exclusion behavior separately. Prepare the bounded smoke input files in this phase so local development never requires the full dataset.

**Exit:** reviewed input/output schemas and reference numerical examples, with every intended departure listed in section 9.

### Phase 2 — Prove the dependency, LOO and memory approach

Provision a pinned Python/CmdStan development environment, compile an unchanged small model and build the LOO adapter. Compare R and Python on identical likelihood arrays, including correlated chains, high Pareto-k values and rejection of the original best model. Prototype selective scalar extraction and bounded likelihood storage.

**Exit:** every LOO output field and winner-selection step agrees within declared tolerances on shared small inputs; smoke-scale memory/disk use is measured and larger allocations are estimated without making them. Resolve API/algorithm differences here before scaling up. Actual full-size memory benchmarking is a production gate in phase 6, not a prerequisite for local smoke runs.

### Phase 3 — Port deterministic calculations

Implement transforms, both density definitions where needed, classification, summary adaptation and typed placeholders. Keep fitting-likelihood and presentation-density constants explicitly distinct.

**Exit:** deterministic fixtures match R, including exact labels, strict thresholds, density grids/scaling, missing diagnostics and stable failure schemas. Run focused Python tests and Black on changed Python files.

### Phase 4 — Integrate candidate fitting and study workflows

Implement model-specific initializations, sampling/subsampling, compact draw extraction, expected-error handling and both model families. Unit-test workflow decisions with an injected runner, then run small real CmdStan integrations.

**Exit:** one flexible and one modes study produce all six outputs on the smoke fixture; all four modes share the same sample, positive rows remain present, and failures in any candidate are handled consistently. This is the first local workflow milestone; avoid expensive production-size runs during development.

### Phase 5 — Add adapters, batching, CLI and staged writes

Implement explicit input formats, catalogue row tracing/QC, deterministic batching, `--list-batches`, validation, streamed output and cleanup. Deliver the `--profile smoke` CLI and bundled-input defaults early in this phase. Test study/row caps before allocation, sequential chains/candidates, both input adapters, deterministic subsets, output isolation, and memory-preflight rejection. Include an oversized input study to verify that smoke loads/transforms only the capped rows; test `--index 0` cannot bypass these limits. Test boundaries at/exceeding 50,000 rows with lightweight fixtures, the first small study, tie ordering, invalid/empty indices, all-failed batches and interrupted writes.

**Exit:** both smoke commands work on the 16 GB Mac with the specified caps, all six outputs, at least one successful fit per family, and measured sampling/post-processing peak memory below 2 GiB. Report first-run compilation and repeat-run time separately. For production planning, no eligible studies disappear or repeat; retained source identity survives chunking; successful outputs round-trip through Python and R TSV readers; failed runs preserve prior complete outputs.

### Phase 6 — Validate parity and downstream use

Separate two kinds of evidence:

- **Identical inputs/draws:** use the same sampled rows, CmdStan CSVs and likelihood arrays. Require exact identities, labels, schemas, candidate order and winner/filter decisions. Start numeric checks at `rtol=1e-8`, `atol=1e-10` for deterministic transforms/densities and ordinary moments; calibrate and document summary-diagnostic/PSIS tolerances against shared arrays before using them. Investigate discrepancies instead of broadly loosening thresholds.
- **Independent sampling:** use the same Stan sources/toolchain/settings and explicit row subset. Compare posterior summaries and transformed-score distributions against Monte Carlo uncertainty. Matching seed integers alone does not establish identical draws. Require stable winner agreement for clearly separated candidates; flag near-tied comparisons for examination.

Exercise downstream readers in a clean, uncached fixture directory. The original QC scripts require specific attention: [both parameter readers](/Users/zakirov/Documents/GitHub/ot_perturbation_catalog_mave_experiments/analysis/mave_gmm_flexible_qc.R:40) search for `params.*.tsv`, while the writer emits `fitted_params{index}.tsv`; [the modes reader](/Users/zakirov/Documents/GitHub/ot_perturbation_catalog_mave_experiments/analysis/mave_gmm_modes_qc.R:52) comments out the posterior/density loads but still uses the resulting objects. Fix these in the repository-owned consumer or compatibility harness before claiming end-to-end compatibility. Do not mask them with cached results or change the established output filename to match the broken pattern.

Run a representative large-study check for peak RAM, scratch space and runtime on a suitably provisioned host, not the 16 GB development Mac. This gates production migration, not delivery of the local smoke workflow. Update the MaveDB README with environment setup, smoke commands first, measured local resource use, explicit production input examples, failure/QC interpretation and scheduler-neutral invocation. Provide corrected LSF migration examples only if retaining that execution route; derive their array bounds from the Python batch plan. Keep the R reference available until acceptance.

**Exit:** both model families pass deterministic and sampling comparisons, downstream fixtures load successfully, and large-study resource requirements are measured. Only then switch production batch commands. Nextflow/Slurm wiring remains a follow-up to this port.

## 9. Deliberate departures and completion checklist

The initial port deliberately changes operational behavior in these documented ways: repository-owned code/dependencies, canonical-input support, an explicit bounded smoke profile, 1-based positive batching, reproducible optional seeds, stable nullable headers, complete study-failure handling, source-row tracing and bounded/staged output handling. Repairs to the downstream readers are compatibility work. Scientific quirks described in section 4 and the modes comparison sequence remain preserved; reduced rows/draws are explicit smoke settings and do not alter production defaults.

- [ ] All five Stan models and their include are retained without scientific changes.
- [ ] Bundled smoke inputs and both real-sampler smoke commands run on the 16 GB Mac, with one study, at most 500 rows, sequential chains/candidates, separate outputs and measured memory/time.
- [ ] Smoke limits apply before large allocations and full-study post-processing; production remains an explicit separate profile with its original defaults.
- [ ] Canonical and legacy inputs have explicit, validated selection/identity contracts.
- [ ] Production sampler settings, candidate initializations and shared subsamples match the R workflow.
- [ ] All six files have declared schemas for success, mixed failure and all-failed batches.
- [ ] Both summary diagnostics and every LOO field pass shared-input parity checks.
- [ ] Modes selection preserves the positive exclusion, `p_loo < 12` cutoff and compare-before-filter behavior.
- [ ] Batch planning covers every eligible score set exactly once, with positive indices starting at 1.
- [ ] Large-study memory and disk use are measured, with no hidden reduction in sampling fidelity.
- [ ] Clean R reader tests pass after the documented consumer fixes.
- [ ] Setup, CLI examples, deviations and acceptance evidence are documented before production migration.
