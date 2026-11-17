# MaveDB Nextflow workflow

This DSL2 workflow implements raw-data ingest, canonical processing, optional
GMM fitting, compatibility tables, basic QC and per-study reports. Verified
completed outputs are reused by default; regeneration requires `--rerun`.

```text
pinned Zenodo ZIP / local ZIP -> ingest -> raw dump + manifest
                                             |
local / downloaded taxonomy ---------> process -> scores + metadata + audit
                                             |                  |
                        --compatibility_tables: tables          |
                                             |                  |
                                  --run_qc: basic QC             |
                                                                |
                  --run_gmm: frozen plan -> batch fits -> collection -> merge
                                                                         |
                                                     --run_reports: reports
```

Taxonomy acquisition remains independent of ingest and processor code. The
workflow passes exact durable result paths between stages; published paths are
convenient views of a consistent result set. Biological validation benchmarks,
homology analysis, other reference acquisition, and cloud/database loading
remain separate commands or future pipeline stages.

## Module layout

`main.nf` connects the stage workflows. Each module owns its actual commands,
staged file inputs, native outputs, resource labels, and reuse-or-execute routing:

| Module | Responsibility |
|---|---|
| `modules/ingest.nf` | Download, validate and unpack the source archive |
| `modules/taxonomy.nf` | Fetch or import the taxonomy snapshot |
| `modules/process.nf` | Produce canonical scores, metadata and audit |
| `modules/tables.nf` | Generate optional compatibility tables |
| `modules/gmm.nf` | Plan, fit and collect GMM batches |
| `modules/merge.nf` | Package and merge completed fits |
| `modules/qc.nf` | Generate basic QC |
| `modules/reports.nf` | Generate per-study reports |

`lib/MavedbCache.groovy` provides the shared interface to the Python cache policy.
`pipeline_runtime.py` resolves reusable results, registers successful producer
outputs, and publishes a consistent result set. It does not execute scientific
commands. Producer scripts validate their own data contracts; registration
checks the module's required files and records the immutable output inventory.
Add or change a scientific command in its `.nf` module, keeping shared cache
policy independent of how the command runs.

## Local setup

Nextflow supports macOS and Linux and uses its local executor without a cluster
or container runtime. Install Bash 3.2+ and a supported Java version (currently
17–26; Java 21 LTS is a suitable choice), then install Nextflow as described in
the [official installation guide](https://docs.seqera.io/nextflow/install/).
This workflow requires Nextflow 24.10 or later and Python 3.10 or later; use the
repository's Python environment for development.

Install commands below are documentation, not an automatic setup step:

```bash
# After installing a supported Java version:
curl -fsSL https://get.nextflow.io | bash
mkdir -p "$HOME/.local/bin"
mv nextflow "$HOME/.local/bin/nextflow"
export PATH="$HOME/.local/bin:$PATH"

# Run from the repository root; includes the processor's Python dependencies.
.venv/bin/python -m pip install -r data_sources/mavedb/requirements.txt
```

The workflow defaults to the repository's `.venv/bin/python`, resolved from the
pipeline directory independently of the launch directory. For another machine,
set `--python /absolute/path/to/python`. Nextflow itself is installed separately
from that Python environment. Record the exact Nextflow and Java versions used
when beginning validation.

## Running later

From the repository root, ingest and process the fixed release currently used
by `scripts/download_zenodo.py`:

```bash
nextflow run data_sources/mavedb/pipeline -profile local
```

This runs canonical processing by default. Add `--ingest_only` to stop at the
raw bundle without taxonomy acquisition or processing. Repeat the command to
reuse completed results; `-resume` is optional and does not control the durable
result policy described below.

## Reusing results and explicit regeneration

Keep the same `--outdir` to continue one saved analysis. Each stage validates
its saved output hashes and the identities of the files it actually consumes.
An unchanged invocation reuses completed outputs without fitting again, even
after Nextflow's disposable work and task-cache directories have been removed.
Missing requested stages are generated. Adding compatibility tables, QC or
reports does not by itself invalidate completed fits.

```bash
# Reuse completed results and generate any missing requested GMM products.
nextflow run data_sources/mavedb/pipeline -profile local \
  --outdir /absolute/path/to/results --run_gmm

# Repackage and merge saved fits without fitting again.
nextflow run data_sources/mavedb/pipeline -profile local \
  --outdir /absolute/path/to/results --rerun merge

# Add source compatibility tables and basic QC, reusing canonical data and fits.
nextflow run data_sources/mavedb/pipeline -profile local \
  --outdir /absolute/path/to/results --run_qc

# Explicitly regenerate more than one stage.
nextflow run data_sources/mavedb/pipeline -profile local \
  --outdir /absolute/path/to/results --rerun gmm,merge,reports

# Regenerate all stages enabled for this invocation, including GMMs.
nextflow run data_sources/mavedb/pipeline -profile local \
  --outdir /absolute/path/to/results --run_gmm --rerun all
```

Repeat the original source, model, cohort and sampler arguments in these
commands when they differ from the defaults. An individual `--rerun` selector
enables its stage and required dependencies. Dependencies with compatible saved
results are reused. `--rerun all` regenerates only enabled stages; it does not
enable optional analyses. Unknown names and incompatible selections fail before
execution.

| Selector | Product | Enable without forcing regeneration |
| --- | --- | --- |
| `ingest` | Verified, unpacked source archive | Default, or `--ingest_only` |
| `taxonomy` | Taxonomy snapshot | Default except ingest-only/external processed input |
| `process` | Canonical scores, metadata and audit | Default |
| `tables` | Raw-score singles, doubles and wide compatibility tables | `--compatibility_tables` |
| `gmm` | Frozen plan, complete batch fits and checked collection | `--run_gmm` |
| `merge` | Packaged batches and merged GMM tables | `--run_gmm` |
| `qc` | Basic source-score and dataset QC figures | `--run_qc` (also enables `tables`) |
| `reports` | Per-study GMM reports | `--run_reports` (also enables GMM products) |

Fitting includes its existing model diagnostics. The separate `qc` stage is
basic source-data QC, rather than a rerun of those diagnostics.

A new `--rerun gmm` starts a fresh attempt for the whole selected fitting cohort;
individual batch or study rerun selectors are not supported. To continue an
interrupted attempt, repeat its arguments **without another `--rerun gmm`**.
Its frozen code and plan, and each completed batch, are stored durably. Only
unfinished batches run, independently of the original work directories.

Changed scientific settings, model/cohort selection or upstream input content
never silently trigger regeneration of a completed stage. The pipeline stops
and identifies the required explicit reruns. Changes to implementation code
alone retain the saved result's provenance and reuse compatible completed
outputs. A registered result with missing or corrupt files is an error requiring
an explicit rerun; it is not treated as a cache miss that authorizes fitting.

When regenerated upstream files change downstream inputs, the pipeline stops
before publishing an inconsistent combination and lists affected selectors,
such as `--rerun gmm,merge,reports`. Completed regeneration work remains
checkpointed for the next invocation. The preceding published result set stays
available until a consistent replacement is ready. Concurrent writers to one
output directory are rejected.

The durable registry and immutable result generations live under
`<outdir>/.pipeline/`. Published `ingest/current`, `processed`, `compatibility`,
`gmm`, `qc`, `reports`, and `references/taxonomy` paths are relative symlinks
through `.pipeline/current` to an immutable publication view. Switching that
single link makes a verified replacement result set visible atomically.
Preserve or transfer the **whole output directory**, including
`.pipeline` and symlinks. Deleting `.pipeline` removes the durable cache;
deleting Nextflow work/cache does not. Older generations are retained, so allow
space for explicit reruns. There is no automatic garbage collection.

Existing canonical processing outputs and packaged GMM results in the output
directory can be adopted after their manifests, receipts, hashes and coverage
checks pass. GMM adoption verifies retained batches and merged tables and does
not require the old work directories or invoke the sampler. Pre-migration
publication directories are preserved under `.pipeline/legacy/` when replaced
with publication links. Unverifiable existing results require an explicit
rerun; arbitrary file existence does not count as completion.

The default URL is the February 2026 dump in versioned Zenodo record `18511521`;
there is no automatic switch to a newer record version. Override it with
`--zenodo_url URL`. For a direct versioned Zenodo file URL, the downloader reads
the record's public `/api/records/<id>` metadata, matches the URL's filename, and
requires its MD5 and byte size. It checks those values for fresh downloads and
cache hits before passing the archive downstream. URLs must use the
`/records/<id>/files/<filename>` form so the version and archive can be
identified unambiguously.

`--archive_sha256` remains optional and is checked by the ingest step. It
supplies an additional SHA-256 check when you have an independent expected
value; the workflow always records the archive's computed SHA-256 in its
manifest. Zenodo's MD5 is used as its published transfer-integrity check and
cache identity; it is not a cryptographic authenticity guarantee.

For offline development, use a local ZIP:

```bash
nextflow run data_sources/mavedb/pipeline -profile local \
  --archive /absolute/path/to/small-mavedb-dump.zip \
  --taxonomy /absolute/path/to/taxcat.dmp \
  --outdir /absolute/path/to/development-results
```

For URL-based runs, the verified ZIP is cached under
`<outdir>/downloads/md5-<zenodo-md5>/<filename>` by default. The record metadata
request is made on each download task execution; the archive transfer is skipped
when the checksum-keyed cache entry still matches both the published size and
MD5. A mismatched entry is discarded and downloaded again. Set `--download_dir
PATH` to use a different persistent cache directory. The verified cache file is
copied into the Nextflow task workspace for validation and ingestion.

A two-row synthetic fixture is supplied in `assets/ingest-smoke/`. To prepare
and run it later, from the repository root:

```bash
repo="$PWD"
cd data_sources/mavedb/pipeline/assets/ingest-smoke
zip -r /tmp/mavedb-ingest-smoke.zip main.json csv
cd "$repo"
nextflow run data_sources/mavedb/pipeline -profile local \
  --archive /tmp/mavedb-ingest-smoke.zip \
  --ingest_only \
  --outdir data_sources/mavedb/output/ingest-smoke
```

This fixture exercises raw ingest only; it is not biological validation data or
a fixture for processing/GMM stages.

The local-archive branch does not schedule a download task. A small fixture must
contain valid `main.json` metadata and at least one `csv/*.scores.csv` with a
`score` column. Metadata can be a list of experiment sets or an object with an
`experimentSets` list; nested experiments use `scoreSets`, matching the existing
processor. The files may be at the ZIP root or inside one enclosing dump tree.
Ingest does not sample, filter, or otherwise modify the raw records. Using the
full release URL on a fresh output directory still downloads and extracts the
entire release, even locally.

To continue using durable outputs, repeat the command with the same `--outdir`.
Nextflow's separate `-resume` facility can reuse matching task executions while
its `.nextflow/` metadata and work directories remain available; explicit
`--rerun` still forces the selected computation. The default work location is
`data_sources/mavedb/tmp/nextflow/work`; override it with Nextflow's
`-work-dir /absolute/path` option. Neither this work directory nor the launch
directory identifies the durable saved analysis.

## Canonical processing

The `process` stage wraps the existing `scripts/process_mavedb.py`, preserving its
dataset selection and score-export rules. For every included dataset, all
original variant records remain in `scores.tsv`, including duplicate HGVS
observations, missing scores, and invalid scores retained in `score_raw`.
Nonqualifying datasets remain excluded from both canonical tables.

Taxonomy is an explicit task input. Supply `--taxonomy /path/to/taxcat.dmp` to
use a local snapshot. Otherwise the `taxonomy` stage downloads NCBI's live
taxcat archive in its own task, extracts `categories.dmp`, and publishes it as
`references/taxonomy/taxcat.dmp` with a checksum manifest. Processing always
passes this file explicitly to the processor, and never uses its implicit
repository cache or downloads inside the processing task. For repeatable runs,
reuse the saved snapshot by default; `--rerun taxonomy` can acquire updated
NCBI content and require explicit dependent reruns. Override the source with `--taxonomy_url URL`;
this is ignored when `--taxonomy` supplies an existing file. No other references
are downloaded yet.

The downloader can also be run independently of Nextflow and processing:

```bash
.venv/bin/python data_sources/mavedb/scripts/download_taxonomy.py \
  --output /absolute/path/to/new-taxonomy-bundle
```

Its manifest records both the downloaded archive identity and the extracted
taxonomy identity. Nextflow calls `process_mavedb.py` directly with the resulting
taxonomy file; there is no separate processing adapter.

Processing writes `<outdir>/processed/` containing:

| File | Contents |
| --- | --- |
| `scores.tsv` | Canonical variant records, keyed by `(dataset_id, source_row)` |
| `metadata.tsv` | One row per included dataset |
| `dataset_processing.tsv` | One row per input score file, inclusion/exclusion reason and available row counts |
| `processing_manifest.json` | Input identities, taxonomy source/hash, code hashes, runtime versions, output hashes and reconciled totals |

Audit reasons are `explicitly_excluded`, `missing_metadata`, `missing_targets`,
and `no_eligible_analysis_rows`. `source_rows` is blank for files excluded before
reading their variants; `exported_rows` is zero for every exclusion. Ingest's
metadata/file reconciliation is not implied by the audit: metadata records with
no corresponding CSV are outside the processor's file-based selection.

Before publishing, the processor verifies canonical score columns, unique
metadata IDs, score-to-metadata joins, sequential `source_row` values from 1 to
N per dataset, metadata row-count agreement, and audit coverage of every input
score file. A failed check prevents the processed directory from being emitted.
The raw ingest bundle is read-only. Compatibility tables are not generated by
the canonical stage. Add `--compatibility_tables` to generate
`<outdir>/compatibility/maves/{singles,doubles,wide}.tsv` independently. These
tables aggregate raw source scores and are distinct from GMM-derived summaries.
The stage reuses the processor's target selection, positional offsets and
aggregation helpers, records its own completion manifest, and leaves canonical
scores, metadata and their processing manifest unchanged.

The processor's `--ingest-manifest PATH` option links the raw provenance and
automatically enables the dataset audit, validation, and processing manifest.
It requires an explicit `--taxonomy` file; `--taxonomy-source URI` records its
original source. These checks run on staged outputs before publication. A new
output directory is published by rename; updating an existing directory uses
separate file replacements. Existing standalone processing arguments remain
available, and `--audit` also validates the canonical tables without requiring
an ingest manifest.

For an entirely offline processing smoke run, package the supplied four-dataset
synthetic fixture (from the repository root):

```bash
repo="$PWD"
cd data_sources/mavedb/pipeline/assets/processing-smoke
zip -r /tmp/mavedb-processing-smoke.zip main.json csv
cd "$repo"
nextflow run data_sources/mavedb/pipeline -profile local \
  --archive /tmp/mavedb-processing-smoke.zip \
  --taxonomy "$repo/data_sources/mavedb/pipeline/assets/processing-smoke/taxcat.dmp" \
  --outdir data_sources/mavedb/output/processing-smoke
```

Expect one included dataset with four score rows, plus three exclusions:
explicit exclusion, absent targets, and absent analysis-eligible rows. This
synthetic fixture is only for software checks, not scientific validation.

## GMM batch fitting

Fitting is opt-in. Add `--run_gmm` to the ingest/processing command to continue
into GMMs, or use an existing canonical input directory directly:

```bash
nextflow run data_sources/mavedb/pipeline -profile local \
  --processed_dir /absolute/path/to/processed --run_gmm \
  --gmm_study_ids urn:mavedb:00000001-a-1 --outdir /absolute/path/to/results
```

`--processed_dir` must initially contain `scores.tsv` and `metadata.tsv`; these
files are imported into a durable canonical snapshot. It bypasses raw ingest,
taxonomy and source processing and is useful when transferring processed data
to HPC. Repeating the same arguments reuses the verified snapshot even if the
original input directory is no longer available. If the supplied canonical
files change, explicitly refresh their snapshot with `--rerun process` (or
`--rerun all`); in this mode `process` imports the supplied tables rather than
running source processing. Include downstream selectors such as
`--rerun process,gmm,merge` when their consumed data change.

External processed input also supports GMM reports with `--run_reports`. Generating source
compatibility tables or basic QC requires the matching raw archive and taxonomy;
`--processed_dir` cannot be combined with `--compatibility_tables`, `--run_qc`
or `ingest`, `taxonomy`, or `tables` reruns. `--ingest_only` cannot be combined with fitting or other
analyses. Omitting optional analysis flags keeps ingest/processing behavior.

The GMM stages use the existing `fit_gmm_batch.py` entry point:

1. Planning snapshots the five GMM library modules, CLI and Stan
   sources, and writes a deterministic plan. It records input/code SHA-256
   identities, model/settings, selected studies and their whole-study batch
   assignments. Planning does not import CmdStan or fit models.
2. Fitting runs one positive batch index per task, with identical
   selection/settings. It verifies the frozen input/code/settings identities
   before fitting and checks the resulting assignment. Each task retains its
   native six tables, QC, and a portable result summary. Each successful task is
   registered in durable storage before downstream collection.
3. Collection verifies every planned batch arrived exactly once,
   every assigned study was processed, and all required tables/QC match their
   recorded hashes.
4. The independent `merge` stage puts each batch's outputs under
   `maves/<model>/batches/batch<N>/`, merges each TSV family with one header in
   batch-index order under `maves/<model>/merged/`, and publishes the complete
   result tree.

Expected study-level numerical failures (such as zero MAD or degenerate scores)
retain placeholder rows, failure summaries and QC in the published batch
outputs, matching the standalone fitter. Pass `--gmm_fail_on_study_error` to
make study-level failures exit nonzero and block collection. Successful partial
batches remain in durable storage for continuation; the final GMM
bundle is emitted only after all batches complete. Sampler/convergence/Pareto-k
warnings remain in QC and do
not become scientific acceptance gates or change model selection.

Defaults are the `production` profile, `modes` model, 50,000 eligible variants
per batch, and deterministic dataset-derived per-study seeds, e.g., `10101` for `urn:mavedb:00000001-a-1` (override with
`--gmm_seed`). Control selection with `--gmm_study_ids` (a comma list
or a list in a params file) and `--gmm_max_studies`. Configure bounded runs with
`--gmm_max_study_rows`, `--gmm_max_fit_rows`, `--gmm_chains`,
`--gmm_parallel_chains`, `--gmm_iter_warmup`, `--gmm_iter_sampling`, and
`--gmm_memory_budget_mb`. Unspecified sampler settings retain the existing
profile defaults. Whole studies stay together; large studies occupy a batch
alone, and row count is only an approximate cost measure.

Provision the pinned `requirements-gmm.txt`, CmdStan 2.37.0, and a compatible
C++ toolchain on the execution host before fitting (see the main MaveDB README).
The pipeline does not install or download toolchains. Compiled models use the
existing locked cache, configurable with `--gmm_cache_dir`. Batch processes
limit BLAS/OpenMP threads to one and request CPUs matching parallel chains.
Memory requests are the configured fitter budget plus 1,024 MiB headroom;
production defaults therefore request 17 GiB, which may exceed laptop capacity.
These are initial requests, not measured full-study requirements. Each task
still scans the canonical inputs; concurrency can increase shared-filesystem I/O.

For a bounded real-sampler check, use the included 500-row processed fixture:

```bash
nextflow run data_sources/mavedb/pipeline -profile local \
  --processed_dir data_sources/mavedb/pipeline/assets/gmm-smoke --run_gmm \
  --gmm_profile smoke --gmm_model flexible \
  --outdir data_sources/mavedb/output/gmm-nextflow-smoke
```

This fits one study using two sequential chains, 200 warmup/200 retained draws
per chain, a 500-row cap, and a 2,048 MiB fitter budget. It uses the supplied
inputs explicitly; the `smoke` profile does not substitute different data.

For Slurm use `-profile slurm -c /absolute/path/to/site.config`. Downloads stay
on the launch host; other stages use Slurm. Supply account/partition/time and
storage settings in the site config. Source inputs, work, output, Python,
CmdStan and compiled cache must be accessible on execution nodes. Use
`--gmm_max_parallel` to bound fitting concurrency (default 2), and
`--gmm_scratch true` to use Nextflow task scratch when available. Sampling and
likelihood temporaries live under that task's working directory. Allow ample
temporary disk space, preserve `.pipeline/` for continuation, and benchmark full-size studies
before increasing concurrency. Slurm execution has not been exercised here.

Outputs are published beneath `<outdir>/gmm/`: `batch_plan.json`,
`completion.json`, and per-batch receipts in `batches/`. This published directory
is a relative link to a durable generation. Native tables and QC
are under `maves/<model>/` (or `smoke/maves/<model>/` for smoke):

```text
maves/<model>/
├── batches/
│   ├── batch1/
│   │   ├── score.tsv
│   │   ├── density.tsv
│   │   ├── fitted_params.tsv
│   │   ├── posterior.tsv
│   │   ├── fit_summary.tsv
│   │   ├── loo.tsv
│   │   └── qc/
│   └── batch2/
│       └── ...
└── merged/
    ├── score.tsv
    ├── density.tsv
    ├── fitted_params.tsv
    ├── posterior.tsv
    ├── fit_summary.tsv
    ├── loo.tsv
    └── qc/
```

Each batch folder contains its own six native tables and QC files. `merged/`
contains the corresponding tables concatenated in batch-index order, with one
header per file. The batch `score.tsv` files retain `dataset_id`; the merged
`score.tsv` omits that redundant column. `completion.json` records each merged
table's row count, size, and SHA-256. The frozen planning checkpoint and code
snapshot are retained separately in the durable registry.
The top-level `gmm/batches/batch<N>.json` files are receipts; the
`maves/<model>/batches/batch<N>/` directories contain the batch data. The
packaging process is separate from model fitting. Use `--rerun merge` after a
packaging change to rebuild it from verified retained fits. Consumers should use
the emitted bundle and its manifest rather than combine files from different
result generations. Use a new `--outdir` to retain another model/cohort as a
separate saved analysis, or explicitly rerun the affected stages.

## Optional QC and reports

Add `--run_qc` to generate basic source-score and dataset figures. This enables
the separate compatibility-table stage because `analysis/mave_qc.py` requires
`singles.tsv`; doubles and wide tables are not inputs to QC. Figures are
published at `<outdir>/qc/`. GMM fitting is not a dependency of basic QC.

Add `--run_reports` to generate per-study reports from the packaged GMM model
and canonical metadata, published at `<outdir>/reports/`. It enables the GMM
products, reusing completed fits and merged tables where available. If those
products have never been generated, fitting is required to create them. Use
`--rerun qc` or `--rerun reports` to regenerate only the selected figures using
compatible saved inputs.

## Outputs and checks

The default output root is `data_sources/mavedb/output/nextflow`. The registry
under `.pipeline/` tracks the durable raw bundle, published through
`<outdir>/ingest/current/`:

```text
<outdir>/ingest/current/
├── manifest.json
└── mavedb_dump/
    ├── main.json
    └── csv/
        ├── *.scores.csv
        └── *.counts.csv   # when supplied by the archive
```

URL downloads are cached separately at
`<outdir>/downloads/md5-<zenodo-md5>/<filename>` by default; this directory can
be changed with `--download_dir`.

The manifest records the source URI, archive size and SHA-256, optional supplied
SHA-256, verified Zenodo filename/record/MD5/size (or `null` for a local ZIP),
ingestion time in UTC, relative dump paths, and metadata/score/count-file totals.
The download task separately checks the Zenodo MD5 and advertised byte size
before extraction.
Archive identities remain recorded in the raw manifest. One `--outdir`
represents a saved analysis; explicit regeneration creates new retained result
generations rather than overwriting the old stage files.
The downloaded ZIP remains in the persistent download cache; a task-local copy
is also kept in Nextflow's work directory. Local input archives remain at their
original location. Allow disk space for the ZIP, the extracted work copy, durable
results and retained generations. Do not edit durable or published bundles.

The ingest task rejects unsafe/duplicate ZIP paths, links and special files,
corrupt archives, ambiguous or invalid metadata, duplicate metadata score-set
IDs, missing score files, and score files without a `score` column. It extracts
into a temporary directory and only exposes a bundle after validation succeeds.
All source files in the detected dump tree are retained byte-for-byte. These
checks establish a readable raw dump, not scientific eligibility or complete
metadata-to-file reconciliation; those belong to subsequent processing/QC.

Processes request one CPU and 2 GB RAM (4 GB for processing), with a 12-hour time
request as initial configuration values, not measured requirements. The local
profile allows one
task at a time. Site-specific executor/resource settings can later be supplied
with `-c /path/to/site.config`; ingest has not yet been validated on HPC.

## Validation status

On 2026-11-17, the durable-reuse implementation passed 97 local pipeline
contract tests and four offline Nextflow 26.10.5 integration tests. The latter
exercise automatic reuse after deleting native work/cache, selective and all-stage
reruns with `-resume`, QC, packaged-GMM adoption, and interrupted batch continuation.
A bounded real flexible GMM smoke fit, explicit refit with `-resume`, report-only
generation, and subsequent no-task reuse also completed. The stage-module
refactor passed the four integration cases and a fresh bounded GMM-to-report
run. Slurm and full-size
production runs remain untested for this implementation. Old result generations
are retained; there is no automatic garbage collection.

The pinned 1.8 GB Zenodo archive completed an ingest run on 2026-11-06 before
this checksum integration. Its SHA-256 matched the published raw-bundle
manifest, and all 3,381 extracted files matched the ZIP sizes and CRCs. The new
Zenodo API lookup and MD5-verified network cache have not been exercised here.
The local-archive ingest-to-processing workflow was exercised on the synthetic
processing fixture, including output checks and resume. The full-dump processing
run, live taxonomy acquisition, and HPC execution have not been performed.
The independent taxonomy task was exercised with a tiny local tar archive,
including its handoff to processing and cached resume. Downloader contract tests
also cover ambiguous and corrupt archives without making external requests.
The direct processor integration produced byte-identical canonical tables and
selection audit to the previous adapter on this fixture. Its tests also check
that validation failures preserve an existing publication.

GMM integration was exercised locally with a real 500-row flexible smoke fit
and a separate two-batch run using two 100-row inputs (one a synthetic clone
for orchestration testing). Both produced complete checked result bundles.
Plan/result contract tests cover changed settings/inputs, missing/duplicate
batches, corrupted outputs and study failures. These bounded software checks
do not close the full-size benchmarking or independent R/Python sampling gates.

Run the offline pipeline contract tests, including durable-cache tests, with:

```bash
.venv/bin/python -m pytest data_sources/mavedb/pipeline/tests
```

The pipeline test directory's `pytest.ini` discovers the existing `test_*.py`
and `*_test.py` suites as well as the new `*_tests.py` contract suites. Tests
that launch Nextflow are opt-in and use offline synthetic inputs:

```bash
MAVEDB_NEXTFLOW_INTEGRATION=1 .venv/bin/python -m pytest \
  data_sources/mavedb/pipeline/tests/nextflow_reuse_tests.py
```

The existing processor regression suite is
`data_sources/mavedb/tests/process_mavedb_test.py`. The offline processing smoke
run above exercises Nextflow syntax and channel wiring. Repeat it without
`-resume` to check durable reuse, then exercise an explicit stage rerun before
moving to the full dump. The historical checks above do not certify the new
durable-cache path or production-scale reruns; run the current acceptance checks
on the intended execution host.
