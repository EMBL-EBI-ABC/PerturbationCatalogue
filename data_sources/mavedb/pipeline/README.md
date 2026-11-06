# MaveDB Nextflow workflow

This DSL2 workflow currently implements **raw-data ingest only**:

```text
pinned Zenodo ZIP -> DOWNLOAD_MAVEDB -> INGEST_MAVEDB_ARCHIVE -> raw dump + manifest
local ZIP --------------------------> INGEST_MAVEDB_ARCHIVE -> raw dump + manifest
```

Metadata/score processing, reference acquisition, GMM fitting, validation, QC,
and loading into cloud storage or Elasticsearch are future stages. The named
`MAVEDB_INGEST` workflow emits `bundle`, `dump`, and `manifest` channels for those
stages. Future processes should consume these channels, rather than read files
from the publication directory.

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

# Run from the repository root; only downloads need this Python dependency.
.venv/bin/python -m pip install -r data_sources/mavedb/requirements-ingest.txt
```

The workflow defaults to the repository's `.venv/bin/python`, resolved from the
pipeline directory independently of the launch directory. For another machine,
set `--python /absolute/path/to/python`. Nextflow itself is installed separately
from that Python environment. Record the exact Nextflow and Java versions used
when beginning validation.

## Running later

From the repository root, download the fixed release currently used by
`scripts/download_zenodo.py`:

```bash
nextflow run data_sources/mavedb/pipeline -profile local
```

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

Processes request one CPU and 2 GB RAM, with a 12-hour time request as initial
configuration values, not measured requirements. The local profile allows one
task at a time. Site-specific executor/resource settings can later be supplied
with `-c /path/to/site.config`; ingest has not yet been validated on HPC.

## Validation status

The pinned 1.8 GB Zenodo archive completed an ingest run on 2026-11-06 before
this checksum integration. Its SHA-256 matched the published raw-bundle
manifest, and all 3,381 extracted files matched the ZIP sizes and CRCs. The new
Zenodo API lookup, MD5-verified cache path, and metadata handoff have **not yet
been executed**. No tests were run for this update.

When execution is authorized, first run the archive contract tests:

```bash
.venv/bin/python -m unittest discover \
  -s data_sources/mavedb/pipeline/tests -p 'test_*.py'
```

Then run the local-archive workflow against a small fixture, check the published
manifest and original bytes, repeat with `-resume`, and check checksum-failure
behavior before using the full Zenodo archive. Python tests do not validate
Nextflow syntax or channel wiring; that requires the workflow exercise.
