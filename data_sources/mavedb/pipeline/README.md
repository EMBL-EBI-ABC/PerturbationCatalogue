# MaveDB Nextflow workflow

This DSL2 workflow implements raw-data ingest, canonical processing and
singles/doubles/wide tables, plus optional GMM fitting, basic QC, per-study
reports, GMM validation and UniProt residue checks. Verified completed outputs are reused by default;
regeneration requires `--rerun`.

```text
ingest + references.taxonomy -> process
ingest + references.taxonomy + process -> tables -> optional QC
process -> GMM -> merge -> reports
process + merge + references.variants -> validation
process + local ProtVar residue TSV -> residue_check
```

The public `references` stage acquires or imports taxonomy and prepares both
variant tables before downstream stages start. It runs beside ingest for a
normal source-mode run. Processing consumes its taxonomy snapshot, and
validation consumes its ClinVar and GWAS tables. Homology analysis and
cloud/database loading remain separate commands.

## Module layout

`main.nf` connects the stage workflows. Each module owns its actual commands,
staged file inputs, native outputs, resource labels, and reuse-or-execute routing:

| Module | Responsibility |
|---|---|
| `modules/ingest.nf` | Download, validate and unpack the source archive |
| `modules/references.nf` | Resolve, assemble and register the complete reference bundle |
| `modules/references_download.nf` | Acquire selected ClinVar/GWAS raw sources |
| `modules/references_prepare.nf` | Prepare or import variant tables and package references |
| `modules/taxonomy.nf` | Fetch, import or reuse the taxonomy component |
| `modules/process.nf` | Produce canonical scores, metadata and audit |
| `modules/tables.nf` | Generate the default singles, doubles and wide tables |
| `modules/gmm.nf` | Plan, fit and collect GMM batches |
| `modules/merge.nf` | Package and merge completed fits |
| `modules/qc.nf` | Generate basic QC |
| `modules/reports.nf` | Generate per-study reports |
| `modules/validation.nf` | Compare fitted scores with controls, ClinVar and GWAS |
| `modules/residue_check.nf` | Check MaveDB reference residues against a local ProtVar snapshot |

`lib/MavedbCache.groovy` provides the shared interface to the Python cache policy.
`pipeline_runtime.py` resolves reusable results, registers successful producer
outputs, and checks consistency across the saved result set. It does not execute scientific
commands. Producer scripts validate their own data contracts; registration
checks the module's required files and records the output inventory and hashes.
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

This prepares taxonomy, ClinVar and GWAS references, then runs canonical
processing and the table stage by default. Add `--ingest_only` to stop at the
raw bundle without reference acquisition or processing. Repeat the command to
reuse completed results; `-resume` is optional and does not control the durable
result policy described below.

## Reusing results and explicit regeneration

Keep the same `--outdir` to continue one saved analysis. Each stage validates
its saved output hashes and the identities of the files it actually consumes.
An unchanged invocation reuses completed outputs without fitting again, even
after Nextflow's disposable work and task-cache directories have been removed.
The pipeline generates missing default tables and reuses any compatible fits
when adding tables, QC or reports.

```bash
# Reuse completed results and generate any missing requested GMM products.
nextflow run data_sources/mavedb/pipeline -profile local \
  --outdir /absolute/path/to/results --run_gmm

# Repackage and merge saved fits without fitting again.
nextflow run data_sources/mavedb/pipeline -profile local \
  --outdir /absolute/path/to/results --rerun merge

# Add basic QC, reusing the default tables, canonical data and fits.
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
| `references` | Taxonomy and prepared ClinVar/GWAS bundle | Default except ingest-only and processed-input runs without validation |
| `process` | Canonical scores, metadata and audit | Default |
| `tables` | Raw-score singles, doubles and wide tables | Default for raw-source runs except `--ingest_only` |
| `gmm` | Frozen plan, complete batch fits and checked collection | `--run_gmm` |
| `merge` | Packaged batches and merged GMM tables | `--run_gmm` |
| `qc` | Basic source-score and dataset QC figures | `--run_qc` |
| `reports` | Per-study GMM reports | `--run_reports` (also enables GMM products) |
| `validation` | Control, ClinVar and GWAS validation tables and optional plots | `--run_validation` (also enables GMM and merge) |
| `residue_check` | UniProt reference-residue comparisons and counts | Standard ProtVar TSV present, or `--uniprot_residues PATH` |

Fitting includes its existing model diagnostics. The separate `qc` stage is
basic source-data QC, rather than a rerun of those diagnostics.

A new `--rerun gmm` starts a fresh attempt for the whole selected fitting cohort;
individual batch or study rerun selectors are not supported. To continue an
interrupted attempt, repeat its arguments **without another `--rerun gmm`**.
Its frozen code and plan, and each completed batch, are stored durably. Only
unfinished batches run, independently of the original work directories.

Changed scientific settings, model/cohort selection or upstream input content
never silently trigger regeneration of a completed stage. The pipeline stops
and identifies the required explicit reruns. References and validation compare
the code for their own stages.
A registered result with missing or corrupt files requires an explicit rerun.

When regenerated upstream files change downstream inputs, the pipeline stops
before declaring the result set consistent and lists affected selectors,
such as `--rerun gmm,merge,reports,validation`. Completed regeneration work remains
checkpointed for the next invocation. Each successful stage replaces its own
previous result; downstream directories may still contain results from earlier
inputs until their explicit reruns finish. There is no global atomic replacement
of multiple stages. Concurrent writers to one output directory are rejected.

Results live in ordinary directories: `ingest/current/`, `references/`,
`processed/`, `tables/`, `gmm-fit/`, `gmm/`, `qc/`, `reports/`, `validation/` and
`residue_check/`.
`gmm-plan/` holds the current frozen plan and code; `gmm-batches/batch<N>/`
holds completed checkpoints for that attempt. Each stage directory contains
`pipeline_manifest.json`, identifying its inputs, settings, required outputs and
hashes. Root-level `pipeline_state.json` records current stage registrations and
unfinished attempts, `pipeline_requests/` holds small execution specifications,
and `pipeline_publication.json` summarizes the last consistent result set.

```text
<outdir>/
├── ingest/current/            # raw dump and ingest provenance
├── references/                # taxonomy, source provenance and variant TSVs
├── processed/                # canonical scores, metadata and audit
├── tables/                   # raw-score singles, doubles and wide tables
├── gmm-plan/                 # current frozen plan and fitting code
├── gmm-batches/batch<N>/      # completed batches of the current attempt
├── gmm-fit/                  # verified collected fitting outputs
├── gmm/                      # packaged batches and merged tables
├── qc/                       # when requested
├── reports/                  # when requested
├── validation/               # when requested
├── residue_check/            # when a ProtVar residue TSV is available
├── pipeline_state.json
├── pipeline_requests/
└── pipeline_publication.json
```

Stage replacement uses a temporary staging directory and a rollback copy during
promotion. A failed replacement restores the previous stage. After a successful
replacement, the previous result is removed; automatic version history is not
retained. Preserve or transfer the **whole output directory**, including the
manifests and root bookkeeping. Deleting Nextflow work/cache does not remove
the reusable stage results. Use another `--outdir` or an external backup to keep
an older analysis.

The pipeline can reuse canonical processing outputs and packaged GMM results
already in the output directory after verifying manifests, receipts, hashes and
coverage. GMM reuse checks the saved batches and merged tables. Results that
fail verification require an explicit rerun.

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

To process a local ZIP entirely offline, supply taxonomy and both prepared
variant tables.

```bash
nextflow run data_sources/mavedb/pipeline -profile local \
  --archive /absolute/path/to/small-mavedb-dump.zip \
  --taxonomy /absolute/path/to/taxcat.dmp \
  --validation_clinvar /absolute/path/to/mave_clinvar_variants.tsv \
  --validation_gwas /absolute/path/to/mave_gwas_variants.tsv \
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

## Unified references stage

Every normal source-mode run resolves one complete reference bundle at
`<outdir>/references/`. It contains the saved taxonomy snapshot and prepared
ClinVar/GWAS TSVs, plus the selected raw source bundle or an external-source
descriptor. The normal source mode downloads the pinned Open Targets coding
variants release 25.06, evidence and credible-set releases 25.12, and the
checked-in MaveDB UniProt/Ensembl map. These source tables are approximately
3.4 GB. Drugs, pharmacogenomics, pockets and homology are not fetched by this
stage. Override source endpoints with `--coding_variants_url`, `--evidence_url`,
`--credible_sets_url` and `--taxonomy_url`.

Use `--taxonomy /path/to/taxcat.dmp` to import a local taxonomy snapshot. Pass
`--reference_sources_manifest /path/to/manifest.json` to prepare variants from
an existing verified version-3 source bundle without downloading or copying its
large partitions. Alternatively, pass both `--validation_clinvar` and
`--validation_gwas` to import a prepared TSV pair; a partial pair is rejected.
Those import options skip source download and variant preparation.

Completed references are reused without polling providers. `--rerun references`
regenerates the bundle from verified retained sources. Add
`--refresh_references` to redownload configured remote sources explicitly when
starting `--rerun references`. Supplied taxonomy files, raw-source manifests
and prepared TSV pairs remain fixed imports; refresh applies to the components
that are downloaded. To continue an interrupted attempt, repeat its input and
refresh options without another `--rerun references`, preserving its checkpoints.
Taxonomy reuse preserves the saved source manifest rather than silently fetching
the current NCBI file. Component checkpoints for an interrupted attempt live below
`<references_cache_dir>/<attempt>`; the default cache root is
`data_sources/mavedb/tmp/nextflow/references`. Source acquisition uses the
`reference_download` process label. With `-profile slurm` it runs locally by
default; set `--reference_download_executor slurm` when compute nodes have the
required network access. Preparation resources use the
`reference_preparation` label and can be overridden in site config.

The bundle layout is:

```text
references/
├── pipeline_manifest.json
├── manifest.json
├── taxonomy/taxcat.dmp
├── taxonomy/manifest.json
├── sources/manifest.json       # local/retained source bundle, when present
├── sources/external.json       # external-source descriptor, when needed
└── variants/
    ├── manifest.json
    ├── mave_clinvar_variants.tsv
    └── mave_gwas_variants.tsv
```

The two source manifest alternatives are mutually exclusive; imported TSV
pairs need neither. Manifests record relative product paths, hashes, sizes,
source releases and table row counts. Each reference consumer records only the
files it reads: process/tables depend on `taxonomy/taxcat.dmp`, and validation
depends on the two variant TSVs.

## Canonical processing

The `process` stage wraps the existing `scripts/process_mavedb.py`, preserving its
dataset selection and score-export rules. For every included dataset, all
original variant records remain in `scores.tsv`, including duplicate HGVS
observations, missing scores, and invalid scores retained in `score_raw`.
Nonqualifying datasets remain excluded from both canonical tables.

Taxonomy is an explicit input from `references/taxonomy/taxcat.dmp`. The
`references` stage downloads NCBI's live taxcat archive by default, or imports
the supplied `--taxonomy` snapshot. Processing passes the saved file and its
recorded source label explicitly; it never uses the processor's implicit
repository cache or downloads inside its own task. Use `--rerun references` to
regenerate reference products, and `--taxonomy_url URL` to select the taxonomy
source for that rerun. When `--taxonomy` is supplied, the URL is ignored.

The downloader can also be run independently of Nextflow and processing:

```bash
.venv/bin/python data_sources/mavedb/scripts/download_taxonomy.py \
  --output /absolute/path/to/new-taxonomy-bundle
```

Its manifest records both the downloaded archive and extracted taxonomy identities.

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
The raw ingest bundle is read-only. The separate `tables` stage generates
`<outdir>/tables/maves/{singles,doubles,wide}.tsv` by default. These
tables aggregate raw source scores and are distinct from GMM-derived summaries.
The stage invokes `scripts/generate_tables.py`, reuses the processor's target
selection, positional offsets and aggregation helpers, records its own
completion manifest, and leaves canonical
scores, metadata and their processing manifest unchanged. The canonical stage
passes `--no-tables` to `process_mavedb.py` so each product is generated only by
its owning stage. Standalone processing generates these tables by default.

The processor's `--ingest-manifest PATH` option links the raw provenance and
automatically enables the dataset audit, validation, and processing manifest.
It requires an explicit `--taxonomy` file; `--taxonomy-source URI` records its
original source. These checks run on staged outputs before publication. A new
output directory is published by rename; updating an existing directory uses
separate file replacements. Existing standalone processing arguments remain
available, and `--audit` also validates the canonical tables without requiring
an ingest manifest.

For an offline processing smoke run, use the bundled four-dataset fixture and
its header-only ClinVar and GWAS TSVs. Run from the repository root.

```bash
repo="$PWD"
fixture="$repo/data_sources/mavedb/pipeline/assets/processing-smoke"
smoke_dir="$(mktemp -d)"
(cd "$fixture" && zip -r "$smoke_dir/processing-smoke.zip" main.json csv)
nextflow run data_sources/mavedb/pipeline -profile local \
  --archive "$smoke_dir/processing-smoke.zip" \
  --taxonomy "$fixture/taxcat.dmp" \
  --validation_clinvar "$fixture/mave_clinvar_variants.tsv" \
  --validation_gwas "$fixture/mave_gwas_variants.tsv" \
  --outdir "$smoke_dir/results"
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

External processed input also supports GMM reports with `--run_reports`. This
downstream-only mode does not generate the normal source tables. Generating
singles/doubles/wide tables or basic QC requires a raw-source run with the
matching archive and taxonomy; `--processed_dir` cannot be combined with
`--run_qc` or `ingest` or `tables` reruns. `--ingest_only` cannot be
combined with fitting or other analyses. Omitting optional analysis flags keeps
the default ingest, canonical processing and table-generation behavior.

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
temporary disk space, preserve the stage directories and root bookkeeping for continuation, and benchmark full-size studies
before increasing concurrency. Slurm execution has not been exercised here.

Outputs are published beneath `<outdir>/gmm/`: `batch_plan.json`,
`completion.json`, and per-batch receipts in `batches/`. This is an ordinary
directory containing the current packaged result. Native tables and QC
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
snapshot live at `<outdir>/gmm-plan/`; current completed fitting checkpoints are
at `<outdir>/gmm-batches/batch<N>/`. The verified collected fitting bundle lives
at `<outdir>/gmm-fit/` and can regenerate `gmm/` without fitting again.
The top-level `gmm/batches/batch<N>.json` files are receipts; the
`maves/<model>/batches/batch<N>/` directories contain the batch data. The
packaging process is separate from model fitting. Use `--rerun merge` after a
packaging change to rebuild it from verified retained fits. Consumers should use
the emitted bundle and its manifest rather than combine stages whose input
identities differ. Use a new `--outdir` to retain another model/cohort as a
separate saved analysis, or explicitly rerun the affected stages.

## Optional QC and reports

Add `--run_qc` to generate basic source-score and dataset figures using the
default table-stage outputs. `analysis/mave_qc.py` requires `singles.tsv`;
doubles and wide tables are not inputs to QC. Figures are
published at `<outdir>/qc/`. GMM fitting is not a dependency of basic QC.

Add `--run_reports` to generate per-study reports from the packaged GMM model
and canonical metadata, published at `<outdir>/reports/`. It enables the GMM
products, reusing completed fits and merged tables where available. If those
products have never been generated, fitting is required to create them. Use
`--rerun qc` or `--rerun reports` to regenerate only the selected figures using
compatible saved inputs.

Add `--run_validation` to compare fitted scores with internal controls and the
prepared ClinVar and GWAS references. This enables fitting, packaging and
reference preparation in both raw-source and `--processed_dir` runs, reusing
verified stage outputs where available. To start offline, supply `--taxonomy`
and both `--validation_clinvar` and `--validation_gwas`. To prepare the variant
tables from a local raw source bundle, use `--reference_sources_manifest` in
place of the TSV pair. Results are published at `<outdir>/validation/`.

| Option | Default | Effect |
| --- | --- | --- |
| `--validation_formats` | `png,pdf` | Figure formats; comma-separated `png`, `pdf`, or both |
| `--validation_dpi` | `150` | Raster figure resolution |
| `--validation_no_plots` | `false` | Write the eight TSV tables and manifest without figures |

Validation reuses its saved result when the two reference TSVs, canonical
scores and metadata, packaged model, options and implementation are unchanged.
Changes to any of these require `--rerun validation`; use it to recompute
against the selected prepared references. The stage registers `controls.tsv`,
`clinvar_matches.tsv`, `gwas_matches.tsv`,
`ctrl_roc.tsv`, `clinvar_roc.tsv`, `gwas_roc.tsv`, `auc_summary.tsv`,
`counts.tsv` and `manifest.json`.

## UniProt reference-residue check

Processing runs automatically check reference residues when the local ProtVar
TSV exists at `data_sources/mavedb/references/protvar/protvar_2025_01_residues.tsv`
in the pipeline checkout. The TSV supplies `accession`, `protein_position` and
`protein_seq` without a database connection at runtime. The workflow logs
when it skips the stage because the standard file is missing. To use a
different snapshot, pass `--uniprot_residues PATH`.

To add this check to an existing full run, put the TSV in the standard location
and repeat that run's original command. Compatible saved stages are reused, so
the residue check is the only new task when the previous stages are complete.
Reports are saved beside them in `residue_check/`.

To create the TSV, first create the destination directory from the repository
root. Then open `psql` from the same directory using your ProtVar database
connection settings.

```bash
mkdir -p data_sources/mavedb/references/protvar
```

In `psql`, run the export on one line. `\copy` writes the TSV to the machine
running `psql`.

```sql
\copy (SELECT DISTINCT accession, protein_position, protein_seq FROM rel_2025_01_genomic_protein_mapping) TO 'data_sources/mavedb/references/protvar/protvar_2025_01_residues.tsv' WITH (FORMAT csv, DELIMITER E'\t', HEADER true)
```

Reuse the exported TSV in raw-source runs or check existing canonical tables
in a separate output directory as shown below.

```bash
nextflow run data_sources/mavedb/pipeline -profile local \
  --processed_dir data_sources/mavedb/output/nextflow/processed \
  --outdir data_sources/mavedb/output/uniprot-residue-check
```

In processed-input mode, the check runs without downloads, GMM fitting or
ClinVar/GWAS references. It uses simple protein HGVS records from datasets
with human protein-coding targets, compares the parsed reference amino acid
at `position + offset_uniprot`, and records when it defaults a missing offset
to zero. Accessions are matched exactly. Multi-target datasets and targets with
multiple UniProt identifiers are reported as ambiguous. Missing accessions or
positions and conflicting reference residues are reported separately.

The stage tries `+1` and `-1` relative to the original offset (or its zero
default). It compares the original and shifted match rates on the same distinct
checks, where both reference positions have exactly one residue. A shift must
cover at least 10 checks and 90% of the dataset's baseline checks, reach 95%
agreement, and improve agreement by at least 20 percentage points. When exactly
one direction qualifies, the stage uses the adjusted offset for all checks in
that dataset. The canonical input tables retain the original offsets.

| Output in `residue_check/` | Contents |
| --- | --- |
| `residue_checks.tsv` | One check per dataset, source position and HGVS reference residue. The original position, residue and `status` are retained alongside their `effective_*` values and `offset_applied`. |
| `dataset_summary.tsv` | One row per dataset with eligible checks. Includes original and effective offsets, `offset_shift`, `offset_applied`, baseline and effective counts, and effective agreement and check-coverage fractions. Agreement is blank when zero checks could be compared. |
| `summary.tsv` | Totals distinct checks and source observations by effective status and `offset_applied`. |
| `applied_offset_corrections.tsv` | Paired comparison evidence for each applied shift, including mismatches resolved and matches lost on the shared checks. |

Agreement is the fraction of compared checks that match. Check coverage is the
fraction of all eligible checks that could be compared. Both use the effective
offset and count distinct checks once, even when several observations share a
source position and HGVS reference residue.

The stage manifest records hashes of both canonical tables and the ProtVar TSV.
Pass `--rerun residue_check` to regenerate an existing result or to use a changed
snapshot.

## Outputs and checks

The default output root is `data_sources/mavedb/output/nextflow`. The raw bundle
is stored directly in the ordinary directory `<outdir>/ingest/current/`:

```text
<outdir>/ingest/current/
├── pipeline_manifest.json
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
represents a saved analysis; explicit regeneration replaces the selected stage
after its new result passes validation.
The downloaded ZIP remains in the persistent download cache; a task-local copy
is also kept in Nextflow's work directory. Local input archives remain at their
original location. Allow disk space for the ZIP, the extracted work copy, durable
results and temporary staging/rollback copies during stage replacement. Do not
edit completed stage outputs or their manifests.

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

## Verification

The bundled offline ingest and processing examples exercise archive checks,
reference imports, canonical outputs, dataset selection and the default MAVE
tables. Repeat the processing example with the same output directory and a
different `-work-dir` to verify durable reuse. Use `--rerun tables` to
regenerate the derived tables from the saved inputs.

The [GMM smoke run](#gmm-batch-fitting) exercises planning, CmdStan sampling,
collection and merged tables on one study. Use `--run_reports` to include a
per-study PDF. The [GMM validation record](../analysis/VALIDATION.md) documents
smoke resource use and the remaining production checks.

Full-dump processing, live acquisition of the complete reference bundle and
Slurm execution still need verification on the production execution host.
