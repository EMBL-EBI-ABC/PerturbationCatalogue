# MaveDB Nextflow workflow

This DSL2 workflow implements **raw-data ingest, canonical processing and optional GMM fitting**:

```text
pinned Zenodo ZIP -> DOWNLOAD_MAVEDB -> INGEST_MAVEDB_ARCHIVE -> raw dump + manifest
local ZIP --------------------------> INGEST_MAVEDB_ARCHIVE -> raw dump + manifest
                                                              |
local taxonomy / FETCH_MAVEDB_TAXONOMY --------------------> PROCESS_MAVEDB
                                                              |
                                    scores.tsv + metadata.tsv + audit + manifest
                                                              |
                         --run_gmm: plan -> parallel batch fits -> checked collection
```

Taxonomy is an independent `MAVEDB_TAXONOMY` workflow in `modules/taxonomy.nf`,
with its own download process and standard-library-only `scripts/download_taxonomy.py`
entry point. It runs alongside raw ingest and supplies an explicit input to
processing. It has no dependency on processor code.

Other reference acquisition, scientific validation/QC, and loading
into cloud storage or Elasticsearch are future stages. `MAVEDB_INGEST` emits
`bundle`, `dump`, and `manifest`; `MAVEDB_PROCESS` emits `bundle`, `scores`,
`metadata`, `audit`, and `manifest`. Downstream processes should consume these
channels rather than read files from the publication directory.

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

This now runs canonical processing by default. Add `--ingest_only` to stop at
the raw bundle without taxonomy acquisition or processing. Add `-resume` to
reuse prior successful tasks when their inputs and scripts still match.

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
full release URL still downloads and extracts the entire release, even locally.

To resume, repeat the same command from the same launch directory with `-resume`.
Keep both Nextflow's `.nextflow/` metadata there and the task work directories.
The default work location is `data_sources/mavedb/tmp/nextflow/work`; override it
with Nextflow's `-work-dir /absolute/path` option. Local archive caching uses
Nextflow's default file identity (path, size, modification time); changing the
expected checksum also invalidates the ingest task cache.

## Canonical processing

`MAVEDB_PROCESS` wraps the existing `scripts/process_mavedb.py`, preserving its
dataset selection and score-export rules. For every included dataset, all
original variant records remain in `scores.tsv`, including duplicate HGVS
observations, missing scores, and invalid scores retained in `score_raw`.
Nonqualifying datasets remain excluded from both canonical tables.

Taxonomy is an explicit task input. Supply `--taxonomy /path/to/taxcat.dmp` to
use a local snapshot. Otherwise `FETCH_MAVEDB_TAXONOMY` downloads NCBI's live
taxcat archive in its own task, extracts `categories.dmp`, and publishes it as
`references/taxonomy/taxcat.dmp` with a checksum manifest. Processing always
passes this file explicitly to the processor, and never uses its implicit
repository cache or downloads inside the processing task. For repeatable runs,
reuse the acquired snapshot with `--taxonomy`; a new acquisition without resume
can fetch updated NCBI content. Override the source with `--taxonomy_url URL`;
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
The raw ingest bundle is read-only. Legacy compatibility tables are not generated.

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

`--processed_dir` must contain `scores.tsv` and `metadata.tsv`; it bypasses raw
ingest, taxonomy and processing. It is useful when transferring processed data
to HPC. `--ingest_only` cannot be combined with fitting. Omitting `--run_gmm`
keeps the previous ingest/processing behavior.

The named `MAVEDB_GMM` workflow uses the existing `fit_gmm_batch.py` entry point:

1. `PLAN_GMM_BATCHES` snapshots the five GMM library modules, CLI and Stan
   sources, and writes a deterministic plan. It records input/code SHA-256
   identities, model/settings, selected studies and their whole-study batch
   assignments. Planning does not import CmdStan or fit models.
2. `FIT_GMM_BATCH` runs one positive batch index per task, with identical
   selection/settings. It verifies the frozen input/code/settings identities
   before fitting and checks the resulting assignment. Each task retains its
   native six tables, QC, and a portable result summary in its work directory.
3. `COLLECT_GMM_BATCHES` verifies every planned batch arrived exactly once,
   every assigned study succeeded, and all required tables/QC match their
   recorded hashes. It stages and publishes one complete native result tree.

Study failures make the Nextflow fit task exit nonzero, even though the legacy
standalone fitter normally returns zero with failed-study placeholders. Failure
summaries and QC remain in the failed task work directory. Successful partial
batches remain in work for resume; the final GMM bundle is emitted only after
all batches succeed. Sampler/convergence/Pareto-k warnings remain in QC and do
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
temporary disk space, retain work for resume, and benchmark full-size studies
before increasing concurrency. Slurm execution has not been exercised here.

Outputs are published beneath `<outdir>/gmm/`: `batch_plan.json`,
`completion.json`, per-batch receipts in `batches/`, and native tables/QC under
`maves/<model>/` (`smoke/maves/<model>/` for smoke). A separate planning copy is
published to `<outdir>/gmm-plan/`. Consumers should use the emitted bundle and
its manifest, rather than globbing a reused publication directory that could
contain files from older runs. Changing model/cohort should use a new outdir.

## Outputs and checks

The default output root is `data_sources/mavedb/output/nextflow`:

```text
<outdir>/ingest/sha256-<archive-digest>/
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
Archive identities keep different releases in separate publication directories.
The downloaded ZIP remains in the persistent download cache; a task-local copy
is also kept in Nextflow's work directory. Local input archives remain at their
original location. Allow disk space for the ZIP, the extracted work copy, and
the published copy. Do not edit published raw bundles.

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

Run the offline archive and processing tests with:

```bash
.venv/bin/python -m unittest discover \
  -s data_sources/mavedb/pipeline/tests -p 'test_*.py'
```

The existing processor regression suite is
`data_sources/mavedb/tests/process_mavedb_test.py`. The offline processing smoke
run above exercises Nextflow syntax and channel wiring; repeat with `-resume`
to check task reuse before moving to the full dump.
