# DepMap Co-Essentiality Feature

A self-contained pipeline and Dash web application for computing and exploring
gene co-essentiality networks from DepMap CRISPR screen data.

---

## Table of Contents

1. [Overview](#overview)
2. [Key Concepts](#key-concepts)
3. [Architecture](#architecture)
4. [Directory Layout](#directory-layout)
5. [Data I/O Schema](#data-io-schema)
6. [Pipeline](#pipeline)
   - [Step 1 — Fetch DepMap data](#step-1--fetch-depmap-data)
   - [Step 2 — GLS co-essentiality](#step-2--gls-co-essentiality)
   - [Step 3 — FDR filtering & network CSV](#step-3--fdr-filtering--network-csv)
7. [Running Locally (end-to-end)](#running-locally-end-to-end)
8. [Dash App — Standalone](#dash-app--standalone)
9. [Cloud Deployment](#cloud-deployment)
10. [Integration Guide](#integration-guide)
11. [Known Constraints & Pitfalls](#known-constraints--pitfalls)
12. [Dependencies](#dependencies)
13. [Citation](#citation)

---

## Overview

This feature calculates **gene co-essentiality** — pairs of genes whose CRISPR
knockout essentiality profiles co-vary significantly across cancer cell lines —
using a **Generalised Least Squares (GLS)** regression that accounts for
cell-line covariance structure. The resulting network is served through an
interactive Dash web application.

**Input:** DepMap CRISPR gene effect matrix (cell lines × genes).  
**Output:** A network CSV of statistically significant gene pairs and a Dash UI
for exploring them.

Scientific method: Wainberg et al. (2021), *A genome-wide atlas of co-essential
modules assigns function to uncharacterized genes*. Nature Genetics.

---

## Key Concepts
- **CRISPR screen:** An experiment where every gene in a panel of cancer cell lines is systematically knocked out (disabled) one at a time. The resulting *gene effect score* measures how essential that gene is for cell survival: a strongly negative score means the cell depends on that gene to survive.

- **Co-essentiality:** Two genes are co-essential when their essentiality scores move together across many cell lines — cell lines where gene A is critical also tend to need gene B. Co-essential genes usually function in the same biological pathway or complex, so this relationship can reveal the function of poorly characterised genes.

- **GLS (Generalised Least Squares):** A statistical regression method used here to test whether two genes' essentiality profiles are correlated. GLS is used instead of ordinary regression because cancer cell lines from the same tissue type look similar to each other; GLS corrects for this shared structure so the test is not artificially inflated.

- **FDR (False Discovery Rate):** When testing ~146 million gene pairs simultaneously, many false positives arise by chance. An FDR threshold of 5% means: among all pairs flagged as significant, we tolerate at most 5% being false positives. The Benjamini–Hochberg (BH) method converts raw p-values into FDR-adjusted p-values to enforce this guarantee. **This correction is applied globally**, across all ~146 million unique pairs at once — not per-gene as in Wainberg et al.'s original method. This is an intentional choice: global correction gives one well-defined, symmetric FDR per pair (Wainberg's per-gene method has no defined rule for symmetrizing a pair when the two genes disagree) and keeps only the strongest co-essential signals. See `REVIEW.md` for the full comparison.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        BATCH PIPELINE (monthly)                         │
│                                                                         │
│  DepMap Portal API                                                      │
│        │                                                                │
│        ▼                                                                │
│  step1_fetch_depmap.py ──► CRISPRGeneEffect_<version>.csv               │
│                    ──► depmap_version.txt                               │
│        │                                                                │
│        ▼                                                                │
│  step2_gls_coessentiality.py ──► <prefix>_GLS_p.npy                     │
│                           ──► <prefix>_GLS_sign.npy                     │
│                           ──► <prefix>_genes.txt                        │
│                           ──► <prefix>_metadata.json                    │
│        │                                                                │
│        ▼                                                                │
│  step3_fdr_coessentiality.py ──► depmap_<version>_gls_whole_            │
│                               coessential_network_FDR_10.csv            │
│        │                                                                │
│        ▼                                                                │
│   GCS Bucket  ◄──────────────────── upload (see Cloud Deployment)       │
└─────────────────────────────────────────────────────────────────────────┘
                                │
                                │ reads CSV at startup
                                ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                    ALWAYS-ON WEB APP (GCP Cloud Run)                    │
│                                                                         │
│   coessentiality_feature_pc.py  (registered as a Dash page)             │
│        └─ serves /coessentiality route on the main website              │
└─────────────────────────────────────────────────────────────────────────┘
```

The pipeline runs once per DepMap release (approximately quarterly).  
The app reads the network CSV from GCS once at startup and holds it in memory
for all subsequent requests — there are no per-request file reads.

---

## Directory Layout

```
coessentiality_feature/
├── README.md                              this file
├── requirements.txt                       Python dependencies for pipeline + app
├── pipeline/
│   ├── step1_fetch_depmap.py              downloads CRISPR data from DepMap API
│   ├── step2_gls_coessentiality.py        computes GLS p-values and sign matrix
│   └── step3_fdr_coessentiality.py        applies BH FDR correction, writes network CSV
├── app/
│   └── coessentiality_feature_pc.py       Dash web application
├── documentation/
│   ├── single_gene_workflow.md            callback flow for the single-gene explorer tab
│   └── gene_list_workflow.md              callback flow for the gene-list network tab
├── example_data/
│   └── 594_DoenchJG_A375.txt             example gene list for testing the app
└── required_data/                         local data directory (not committed to git)
    ├── depmap_version.txt
    ├── CRISPRGeneEffect_<version>.csv
    ├── <prefix>_genes.txt
    ├── <prefix>_GLS_p.npy
    ├── <prefix>_GLS_sign.npy
    ├── <prefix>_metadata.json
    ├── depmap_<version>_gls_whole_coessential_network_FDR_10.csv
    └── GO_Biological_Process_2025.gmt         GO:BP gene-set library for local enrichment
```

> **Note:** The `required_data/` directory is excluded from version control. Large
> intermediate files (`.npy`, raw CRISPR CSV) are not stored in the repo.
> Only the final network CSV and metadata files need to be accessible to the
> Dash app at runtime.

---

## Data I/O Schema

This section documents every file exchanged between pipeline steps and the app.
Use it as the contract when modifying any script or when setting up GCS.

### `depmap_version.txt`

| Property | Value |
|---|---|
| Format | Plain text, single line |
| Example content | `DepMap Public 26Q1` |
| Produced by | `step1_fetch_depmap.py` |
| Consumed by | `step2_gls_coessentiality.py` (prefix derivation), `coessentiality_feature_pc.py` (version badge) |

---

### `CRISPRGeneEffect_<version>.csv`

| Property | Value |
|---|---|
| Format | CSV |
| Rows | Cancer cell lines (1,208 in DepMap 26Q1) |
| Columns | Genes in `SYMBOL (ENTREZ_ID)` format, e.g. `A1BG (1)` |
| Index column | `ModelID` (first unnamed column) |
| Values | CRISPR gene effect score (float, typically −3 to +1) |
| Size | ~436 MB |
| Produced by | `step1_fetch_depmap.py` (downloaded from DepMap API) |
| Consumed by | `step2_gls_coessentiality.py` |

> **Pipeline note:** `step2_gls_coessentiality.py` strips the Entrez IDs from
> column names on load (e.g. `A1BG (1)` → `A1BG`). The raw file is only
> needed for the pipeline — it is not used by the Dash app.

---

### `<prefix>_genes.txt`

| Property | Value |
|---|---|
| Format | Plain text, one gene symbol per line, no header |
| Example content | `A1BG\nA1CF\nA2M\n...` |
| Row count | 17,087 (DepMap 26Q1, after NA filtering) |
| Produced by | `step2_gls_coessentiality.py` |
| Consumed by | `step3_fdr_coessentiality.py` |

---

### `<prefix>_metadata.json`

| Property | Value |
|---|---|
| Format | JSON |
| Keys | `n_cell_lines` (row count of the input CRISPR matrix), `n_genes_profiled` (gene count after NA filtering) |
| Produced by | `step2_gls_coessentiality.py` |
| Consumed by | `coessentiality_feature_pc.py` (app stat chips) |

> Avoids the app having to scan the ~440 MB `CRISPRGeneEffect_<version>.csv` or
> count lines in `<prefix>_genes.txt` just to display two numbers at startup.

---

### `<prefix>_GLS_p.npy`

| Property | Value |
|---|---|
| Format | NumPy binary array (`.npy`) |
| Shape | `(n_genes, n_genes)` — e.g. `(17087, 17087)` for DepMap 26Q1 |
| Dtype | `float64` |
| Values | Two-sided GLS p-values for each gene pair; diagonal set to 1.0 |
| Size on disk | ~2.4 GB |
| Produced by | `step2_gls_coessentiality.py` |
| Consumed by | `step3_fdr_coessentiality.py` |

---

### `<prefix>_GLS_sign.npy`

| Property | Value |
|---|---|
| Format | NumPy binary array (`.npy`) |
| Shape | `(n_genes, n_genes)` |
| Dtype | `float64` |
| Values | Sign of GLS regression coefficient: `+1.0` or `−1.0` |
| Size on disk | ~2.4 GB |
| Produced by | `step2_gls_coessentiality.py` |
| Consumed by | `step3_fdr_coessentiality.py` |

---

### `depmap_<version>_gls_whole_coessential_network_FDR_10.csv`  ← **primary app input**

| Property | Value |
|---|---|
| Format | CSV |
| Naming convention | `depmap_<version>_gls_whole_coessential_network_FDR_10.csv` |
| Example name | `depmap_26Q1_gls_whole_coessential_network_FDR_10.csv` |
| Size | ~1.7 MB (DepMap 26Q1) |
| Produced by | `step3_fdr_coessentiality.py` |
| Consumed by | `coessentiality_feature_pc.py` |

**Columns:**

| Column | Type | Description |
|---|---|---|
| `source` | string | Gene symbol (alphabetically first in the pair) |
| `target` | string | Gene symbol (alphabetically second in the pair) |
| `pvalue` | float64 | Raw two-sided GLS p-value |
| `pvalue_adj` | float64 | Benjamini-Hochberg adjusted p-value (FDR), corrected **globally** across all unique pairs (see Key Concepts above) — not per-gene as in Wainberg et al.'s original method |
| `direction` | float64 | Sign of GLS coefficient: `+1.0` = positive co-essentiality, `−1.0` = negative |

> **Important:** This file contains all pairs at FDR ≤ 10%. The app filters
> it at query time to either FDR ≤ 5% or FDR ≤ 10% based on the user's
> dropdown selection. Do not pre-filter to 5% before writing.

---

### `GO_Biological_Process_2025.gmt`

| Property | Value |
|---|---|
| Format | GMT (tab-separated: term name + GO ID, then gene symbols) |
| Size | ~1.4 MB, 5,341 terms |
| Produced by | One-time manual download (see below) — **not** part of the regular pipeline |
| Consumed by | `coessentiality_feature_pc.py` (`_GO_LIBRARY`, loaded once at startup) for local GO:BP enrichment |

**Why this exists:** GO:BP enrichment used to call the live Enrichr API per module
on every "Find modules & annotate" click — slow (~1-2s/module), and a hard
external dependency (rate limits, outages) once deployed to GCP with multiple
concurrent users sharing one outbound IP. This file lets enrichment run **fully
offline** via `gp.enrich()` (gseapy's local hypergeometric test), which is also
noticeably faster in practice than the old per-click Enrichr round trip.

**⚠️ Maintenance — this is a pinned snapshot, not auto-updated:**

Enrichr re-packages a new GO:BP library roughly every 1-2 years (available
versions as of writing: `2021`, `2023`, `2025`, `2026`). This file intentionally
does **not** auto-refresh — silently shifting enrichment results out from under
users without anyone reviewing the change would be worse than being a year
"stale." To bump it deliberately:

```bash
python3 -c "
import gseapy as gp
gp.get_library(name='GO_Biological_Process_2026', organism='Human',
               save='required_data/GO_Biological_Process_2026.gmt')
"
```
Then update `GO_GENE_SET_PATH` in `coessentiality_feature_pc.py` to point at the
new filename. Do this as a deliberate, reviewed decision (e.g. an annual check),
not as part of the routine monthly pipeline run.

**Production note:** `required_data/` is gitignored, so this file must be
uploaded to GCS once (alongside the network CSV / metadata.json — see Cloud
Deployment below). It does **not** need re-uploading on every monthly pipeline
run — only when you deliberately bump the version.

---

## Pipeline

### Step 1 — Fetch DepMap data

**Purpose:** Queries the DepMap portal files API, identifies the latest
`CRISPRGeneEffect.csv` release, and downloads it with MD5 verification.

**Usage:**
```bash
python step1_fetch_depmap.py [--output-dir DIR]
```

| Argument | Default | Description |
|---|---|---|
| `--output-dir` | `./required_data` | Directory to write all output files |

**Outputs:**

| File | Description |
|---|---|
| `CRISPRGeneEffect_<version>.csv` | CRISPR gene effect matrix |
| `depmap_version.txt` | Release name, e.g. `DepMap Public 26Q1` |

**Runtime:** ~5–10 minutes (depends on download speed; file is ~436 MB).

**API endpoint used:** `https://depmap.org/portal/api/download/files`

> **Idempotency note:** The script currently re-downloads data on every run
> regardless of version. Before using in a monthly batch job, add a version
> check to skip the download when `depmap_version.txt` already records the
> latest release.

---

### Step 2 — GLS co-essentiality

**Purpose:** Loads the CRISPR matrix, drops genes with any missing values,
computes the pseudoinverse covariance matrix (GLS), runs OLS regression in
the GLS-transformed space for all gene pairs, and saves the p-value and
sign matrices.

**Usage:**
```bash
python step2_gls_coessentiality.py <input_csv> [--output-dir DIR] [--prefix NAME]
```

| Argument | Default | Description |
|---|---|---|
| `input` | (required) | Path to `CRISPRGeneEffect_<version>.csv` |
| `--output-dir` | `./output` | Directory for output files |
| `--prefix` | input filename stem | Prefix for all output filenames |

**Outputs:**

| File | Description |
|---|---|
| `<prefix>_genes.txt` | Gene list after NA filtering, one per line |
| `<prefix>_GLS_p.npy` | `(n_genes, n_genes)` float64 p-value matrix |
| `<prefix>_GLS_sign.npy` | `(n_genes, n_genes)` float64 sign matrix |
| `<prefix>_metadata.json` | `{"n_cell_lines": ..., "n_genes_profiled": ...}` — read by the app for its stat chips |

**Runtime:** ~4–8 hours on a single core (17,087 gene iterations, each a
vectorised NumPy lstsq call over all other genes simultaneously).

**Memory requirement:** ~16 GB RAM peak. The 17,087 × 17,087 float64 matrices
each consume ~2.4 GB. Cloud Run Jobs should be configured with `--memory 16Gi`.

> **Critical constraint:** This script must run on the **full DepMap gene
> matrix** (all 17,087+ genes). Running on a gene subset causes the
> pseudoinverse covariance matrix to become rank-deficient, which breaks the
> Cholesky decomposition with a `LinAlgError`. This is expected mathematical
> behaviour, not a bug.

---

### Step 3 — FDR filtering & network CSV

**Purpose:** Loads the p-value matrix, sign matrix, and gene list. Extracts the
upper-triangle gene pairs, applies Benjamini-Hochberg FDR correction **globally**
across all ~146 million pairs at once (an intentional departure from Wainberg et
al.'s per-gene FDR — see Key Concepts above), filters at the specified FDR
threshold, and writes the significant network to a CSV file.

**Usage:**
```bash
python step3_fdr_coessentiality.py \
    --gls-p   <prefix>_GLS_p.npy \
    --gls-sign <prefix>_GLS_sign.npy \
    --genes   <prefix>_genes.txt \
    --fdr     0.10 \
    --output  depmap_<version>_gls_whole_coessential_network_FDR_10.csv
```

| Argument | Default | Description |
|---|---|---|
| `--gls-p` | (required) | Path to `_GLS_p.npy` |
| `--gls-sign` | (required) | Path to `_GLS_sign.npy` |
| `--genes` | (required) | Path to `_genes.txt` |
| `--fdr` | `0.05` | FDR threshold for filtering (use `0.10` for the app) |
| `--output` | (required) | Output CSV path |

**Outputs:**

| File | Description |
|---|---|
| `<output>` | Network CSV (see [Data I/O Schema](#data-io-schema)) |

**Runtime:** ~10–20 minutes. The stacking of a 17,087 × 17,087 matrix and
BH correction across ~146 million pairs is memory-intensive.

**Memory requirement:** ~16 GB RAM peak (stacked Pandas Series + multipletests).

> **Always use `--fdr 0.10`** when generating the file for the Dash app.
> The app loads this single file at startup and filters it at query time
> to either 5% or 10% based on the user's selection. Generating a 5%-only
> file would make the 10% view unavailable without a restart.

---

## Running Locally (end-to-end)

Prerequisites: Python 3.10+.

```bash
cd coessentiality_feature/

# One-time setup: create a virtual environment and install dependencies
# (see Dependencies below). Re-run `pip install -r requirements.txt` after
# pulling changes that touch requirements.txt.
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Step 1: download the latest DepMap data
python pipeline/step1_fetch_depmap.py --output-dir required_data/

# Step 2: compute GLS (takes several hours, needs ~16 GB RAM)
python pipeline/step2_gls_coessentiality.py \
    required_data/CRISPRGeneEffect_26Q1.csv \
    --output-dir required_data/ \
    --prefix depmap_26Q1

# Step 3: apply FDR and write the network CSV
python pipeline/step3_fdr_coessentiality.py \
    --gls-p    required_data/depmap_26Q1_GLS_p.npy \
    --gls-sign required_data/depmap_26Q1_GLS_sign.npy \
    --genes    required_data/depmap_26Q1_genes.txt \
    --fdr      0.10 \
    --output   required_data/depmap_26Q1_gls_whole_coessential_network_FDR_10.csv

# Step 4: launch the app
python app/coessentiality_feature_pc.py
# → open http://localhost:8050
```

The version tag (`26Q1`) in step 2 must match the string extracted from
`depmap_version.txt` for the app to find the correct network CSV at startup.

---

## Dash App — Standalone

```bash
python app/coessentiality_feature_pc.py
```

Opens at `http://localhost:8050`.

**What the app expects at startup:**

1. `required_data/depmap_version.txt` — determines which network CSV filename to load.
2. `required_data/depmap_<version>_gls_whole_coessential_network_FDR_10.csv` — the
   network data; loaded once into memory.
3. `required_data/depmap_<version>_metadata.json` — `n_cell_lines` and
   `n_genes_profiled` for the stat chips, written once by step 2. Optional: if
   absent, both chips show `—`.

**Two UI tabs:**

| Tab | Purpose |
|---|---|
| Single-gene explorer | Search one gene (or click a TP53/BRCA1/KRAS quick-search button); shows bar chart, partner table, Cytoscape network (top 20 partners), and a "Download all partners (CSV)" button for the full partner list |
| Gene-list network | Paste a list of genes (or click "Load example gene list"); shows a Cytoscape network of significant pairs, with a degree-of-interaction slider to expand to nearby genes, a co-essential modules table with on-demand GO:BP annotation, and CSV downloads for both the network pairs and the GO:BP terms |

**FDR filter** (top dropdown): applies globally to both tabs; filters the
in-memory DataFrame at query time — no file re-read.

For the full callback-by-callback breakdown of each tab, see
[`documentation/single_gene_workflow.md`](documentation/single_gene_workflow.md) and
[`documentation/gene_list_workflow.md`](documentation/gene_list_workflow.md).

---

## Cloud Deployment

### Data flow on GCP

```
Cloud Run Job (monthly)             Cloud Run Service (always-on)
  pipeline/ scripts                   app/ Dash app
        │                                    │
        │  writes CSV                        │  reads CSV at startup
        ▼                                    ▼
  GCS Bucket: gs://<DEPMAP_BUCKET>/<DEPMAP_BLOB>
```

### Environment variables (app)

| Variable | Example value | Description |
|---|---|---|
| `DEPMAP_BUCKET` | `my-project-depmap` | GCS bucket name |
| `DEPMAP_BLOB` | `depmap_26Q1_gls_whole_coessential_network_FDR_10.csv` | Blob path within the bucket |

The app must be updated to read the network CSV from GCS using these variables
instead of the local `_DATA_DIR` path. See [Integration Guide](#integration-guide).

### What is already built

- [x] `step1_fetch_depmap.py` — fetches from DepMap API with MD5 check
- [x] `step2_gls_coessentiality.py` — GLS computation
- [x] `step3_fdr_coessentiality.py` — FDR filtering and network CSV

### What still needs to be built for full cloud automation

- [ ] **GCS upload** — add to `step3_fdr_coessentiality.py` (or an orchestrator
  script): upload the network CSV and `depmap_version.txt` to the GCS bucket
  after writing them locally
- [ ] **Orchestrator** — a `run_pipeline.sh` or `main.py` that chains
  `step1 → step2 → step3` with correct arguments and exits non-zero on any failure
- [ ] **Separate `requirements.txt` for the pipeline container** — the repo-level
  `requirements.txt` includes app dependencies that are not needed in the pipeline
  container. Create a minimal `pipeline/requirements.txt` with only:
  `numpy`, `pandas`, `scipy`, `statsmodels`, `requests`, `google-cloud-storage`
- [ ] **Dockerfile** for the pipeline Cloud Run Job:
  ```dockerfile
  FROM python:3.11-slim
  WORKDIR /app
  COPY pipeline/ ./pipeline/
  COPY requirements.txt .
  RUN pip install --no-cache-dir -r requirements.txt
  CMD ["python", "pipeline/main.py"]
  ```
- [ ] **Cloud Run Job** — deploy with `--memory 16Gi` for the GLS step
- [ ] **Cloud Scheduler** — monthly trigger, e.g. `0 6 1 * *` (1st of each
  month, 06:00 UTC)
- [ ] **IAM** — Cloud Run Job service account needs
  `roles/storage.objectAdmin` on the GCS bucket

---

## Integration Guide

This section covers the exact changes needed to embed the Dash app as a page
in the **existing GCP Dash website** (a separate repository).

### 1. Prepare the app file

In `app/coessentiality_feature_pc.py`, make the following changes before
copying it to the main website:

**a. Register as a Dash multi-page page** (add near the top, before the layout):
```python
import dash
dash.register_page(__name__, path="/coessentiality", name="Co-Essentiality")
```

**b. Remove the standalone app initialisation** (these belong to the main app):
```python
# REMOVE these three lines:
app = dash.Dash(__name__)
app.title = "DepMap Co-Essentiality Explorer"
# and at the bottom:
if __name__ == "__main__":
    app.run(debug=True)
```

**c. Rename the layout variable:**
```python
# BEFORE:
app.layout = html.Div([...])
# AFTER:
layout = html.Div([...])
```

**d. Replace local file reads with GCS reads:**

`_DATA_DIR` (and the local `required_data/` folder it points to) is read in
**five** places — all of them need to be handled, not just the network CSV:

| Reads from `_DATA_DIR` | What for | Size |
|---|---|---|
| `_resolve_data_path()` → `df_all` | The network CSV — primary app data | ~1.7 MB |
| `_DEPMAP_VERSION` (module-level) | Version badge text, e.g. "DepMap Public 26Q1" | a few bytes |
| `_load_dataset_stats()` | `n_cell_lines` and `n_genes_profiled` from `depmap_<version>_metadata.json` for the stat chips | a few bytes |
| `_load_all_profiled_genes()` | Full gene list from `depmap_<version>_genes.txt` for the search dropdowns | a few hundred KB |
| `_GO_LIBRARY` (module-level) | Local GO:BP enrichment library from `GO_Biological_Process_2025.gmt` | ~1.4 MB |

Replace `_resolve_data_path()` with a GCS read for the network CSV. Install
`google-cloud-storage` and add:
```python
import io
from google.cloud import storage

def _load_network_from_gcs():
    bucket_name = os.environ["DEPMAP_BUCKET"]
    blob_name   = os.environ["DEPMAP_BLOB"]
    client = storage.Client()
    blob   = client.bucket(bucket_name).blob(blob_name)
    return pd.read_csv(io.BytesIO(blob.download_as_bytes()))

df_all = _load_network_from_gcs()
```

For `_DEPMAP_VERSION` and `depmap_<version>_genes.txt`, do the same small-file
GCS read (both are tiny, fine to fetch on every startup).

> **Don't** fetch the 436 MB `CRISPRGeneEffect_<version>.csv` from GCS just to
> compute `n_cell_lines` — that's a lot of bandwidth for one integer. Either
> precompute `n_cell_lines` once (e.g. write it alongside `depmap_version.txt`
> as part of the pipeline) and read that small value from GCS instead, or drop
> the cell-line stat chip in the cloud version — `_load_dataset_stats()`
> already returns `None` gracefully (chip shows `—`) if the file is missing.

**e. The gene-list Textarea fires on blur** (`n_blur`) rather than on every keystroke.
This is already implemented in the standalone app. No change needed — confirm the
callback uses `Input("multi-gene-input", "n_blur")` with `State("multi-gene-input", "value")`.

> **Note:** `dcc.Textarea` in Dash 4.x does not support a `debounce` property.
> The blur-based trigger (`n_blur`) is the correct equivalent — it fires when
> the user clicks outside the textarea.

**f. Add `server = app.server`** (for gunicorn) and **set `debug=False`** —
these are handled by the main website's `app.py`; confirm they are present there.

### 2. Add to the main website

```bash
# In the main website repository:
cp coessentiality_feature_pc.py pages/coessentiality.py
```

Add to `requirements.txt`:
```
dash-cytoscape
networkx
gseapy
google-cloud-storage
```

### 3. Redeploy with environment variables

```bash
gcloud run services update <SERVICE_NAME> \
    --update-env-vars DEPMAP_BUCKET=<bucket>,DEPMAP_BLOB=<blob>
```

### 4. Verify

- Navigate to `https://<your-domain>/coessentiality`
- Search for `TP53` — should return partners and render the Cytoscape network
- Switch the FDR dropdown between 5% and 10% — network size badge should update
- Paste `BRCA1, TP53, KRAS, PTEN` into the gene-list tab — should show edges

---

## Known Constraints & Pitfalls

| Constraint | Detail |
|---|---|
| Full gene matrix required | `step2_gls_coessentiality.py` must receive the complete DepMap matrix (~17,087 genes after NA filtering). Running on a gene subset causes a `LinAlgError` in the Cholesky step — this is expected, not a bug. |
| Memory: pipeline | `step2_gls_coessentiality.py` and `step3_fdr_coessentiality.py` each peak at ~16 GB RAM. Set `--memory 16Gi` on the Cloud Run Job. |
| Memory: app | The network CSV is ~1.7 MB and expands to ~10 MB as a DataFrame. Negligible; no concern for the app container. |
| Runtime: step 2 | ~4–8 hours on a single CPU core. Acceptable for a monthly batch job; configure Cloud Run Job `--task-timeout` accordingly (max 24h). |
| `debug=True` | The app file currently runs with `debug=True` in standalone mode. This enables Werkzeug's interactive debugger — a remote code execution vector on a public URL. Always set `debug=False` in the main website's `app.py`. |
| No per-request file I/O | The network CSV is loaded once at startup. Restarting the Cloud Run service is the mechanism for picking up newly uploaded data; there is no hot-reload. |
| FDR file must be 10% | The app expects the FDR 10% superset file. If only a 5% file is available, the 10% dropdown option will produce incorrect results (showing fewer pairs than expected). |
| `GO_Biological_Process_2025.gmt` must be present | The app loads this at startup for local GO:BP enrichment; if absent, the app fails to start. Pinned, not auto-updated — see [`GO_Biological_Process_2025.gmt`](#go_biological_process_2025gmt) in Data I/O Schema for the maintenance process. |

---

## Dependencies

### Pipeline

| Package | Purpose |
|---|---|
| `numpy` | Matrix operations, `.npy` I/O |
| `pandas` | CSV I/O, DataFrame operations |
| `scipy` | `stdtr` for p-value computation |
| `statsmodels` | Benjamini-Hochberg FDR correction |
| `requests` | DepMap API HTTP calls |
| `google-cloud-storage` | GCS upload (cloud deployment) |

### App

| Package | Purpose |
|---|---|
| `dash` | Web framework |
| `dash-cytoscape` | Interactive network graph component |
| `plotly` | Bar charts |
| `numpy` | Colour interpolation for node colours |
| `pandas` | Network CSV querying |
| `networkx` | Connected-component detection for co-essential module discovery |
| `gseapy` | Local GO Biological Process enrichment (hypergeometric test against a pinned `.gmt` file — no Enrichr API calls at runtime) |
| `google-cloud-storage` | GCS CSV read (cloud deployment only) |

---

## Citation

Wainberg, M., Kamber, R.A., Balsubramani, A. et al.  
*A genome-wide atlas of co-essential modules assigns function to uncharacterized genes.*  
Nature Genetics 53, 638–649 (2021).  
https://doi.org/10.1038/s41588-021-00840-z

DepMap data: Broad Institute DepMap Portal — https://depmap.org/portal/
