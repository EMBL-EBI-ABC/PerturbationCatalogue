# Gene-List Network View — Workflow & GO:BP Annotation

The gene-list view lets users paste a set of gene symbols and explore how they are co-essentially connected. It groups connected genes into **co-essential modules**, then annotates each module with its most enriched GO Biological Process term.

---

## 1. App startup

Before any user interaction, the app loads everything it needs from disk. This runs once when the server boots.

```mermaid
flowchart LR

classDef file  fill:#4C72B0,color:#fff,stroke:none
classDef data  fill:#2E7D52,color:#fff,stroke:none

A[depmap_version.txt]:::file --> B[FDR 10% network CSV\nsource · target · pvalue_adj · corr_genes]:::data
A --> C[CRISPRGeneEffect CSV\n→ number of cancer cell lines]:::data
A --> D[Model.csv\n→ number of cancer subtypes]:::data
B --> E[(df_all loaded into memory\nall callbacks filter this at query time)]:::data
B --> F[all_genes set\nused to validate user input]:::data
```

---

## 2. Gene-list view — full workflow

The view is driven by **three callbacks** that fire at different points.

```mermaid
flowchart TD

classDef user     fill:#4C72B0,color:#fff,stroke:none,rx:16
classDef fast     fill:#2E7D52,color:#fff,stroke:none
classDef problem  fill:#C0392B,color:#fff,stroke:#7B241C,stroke-width:3px
classDef store    fill:#6C757D,color:#fff,stroke:none
classDef out      fill:#E8F5EE,color:#212121,stroke:#2E7D52

A([Paste gene list + set FDR]):::user

subgraph CB1["Callback 1 — runs on every keystroke, no external calls"]
    B[Parse & validate gene symbols]:::fast
    C[Filter network at chosen FDR\nkeep pairs where both genes are in the list]:::fast
    D[Find connected components\nnumber by size  1 = largest\ndrop components < 2 genes]:::fast
    E[Tag each node & edge\nwith its module ID]:::fast
    B --> C --> D --> E
end

A --> B

F[/Cytoscape network rendered/]:::out
G[(Module list stored\ncluster · size · genes)]:::store

E --> F
E --> G

H([Click  Find modules & annotate]):::user
G --> H

subgraph CB2["Callback 2 — runs on button click, one API call per module"]
    I{Module\n> 3 genes?}
    J[Not annotated\ntoo small]:::out
    K[gseapy.enrichr\nHTTP call to Enrichr API\n~1–2 s per module]:::problem
    L[Top GO:BP term\nby adjusted p-value]:::fast
    I -- No --> J
    I -- Yes --> K --> L
end

H --> I

M[/Modules table\ncluster · size · GO term · p-value · adj. p-value/]:::out
J --> M
L --> M

N([Click a row in the table]):::user
M --> N

subgraph CB3["Callback 3 — runs on row click, no external calls"]
    O[Add highlight rules to Cytoscape\nselected module → orange border & edges]:::fast
end

N --> O --> F
```

---

## 3. The GO:BP annotation problem

The red box in Callback 2 is a **live HTTP call to the Enrichr web API**. This is fine for local development but causes serious problems when the app is deployed to Google Cloud Platform.

| Problem | Impact |
|---|---|
| **One API call per module, sequential** | 5 modules × ~2 s = ~10 s blocking. Cloud Run has per-request timeouts — users with larger gene lists will hit them. |
| **Rate limits** | Enrichr limits queries per IP. Multiple concurrent users share one GCP outbound IP and will get throttled or 429 errors. |
| **External dependency** | If Enrichr is unavailable (maintenance, outage), GO annotation silently fails for all users. |
| **No caching** | The same GO:BP library is re-downloaded from Enrichr on every button click, even though it never changes between requests. |

---

## 4. The fix — local enrichment (planned for GCP deployment)

The GO:BP gene-set library is **static** — it only changes when Enrichr releases a new version. The enrichment test is just a **hypergeometric test + Benjamini–Hochberg correction**, which runs in milliseconds. There is no reason to call an external API at all.

```mermaid
flowchart LR

classDef once   fill:#4C72B0,color:#fff,stroke:none
classDef boot   fill:#2E7D52,color:#fff,stroke:none
classDef live   fill:#1D5C3A,color:#fff,stroke:none

subgraph ONCE["One-time — run before deploying"]
    A[gseapy.get_library\ndownload GO_Biological_Process_2025]:::once
    B[Save as go_bp_2025.pkl\nupload to GCS with the network CSVs]:::once
    A --> B
end

subgraph BOOT["At app startup — runs once per instance"]
    C[Load go_bp_2025.pkl from disk]:::boot
    D[(GO sets in memory\nterm → gene set\nshared across all requests)]:::boot
    C --> D
end

subgraph LIVE["At annotation time — replaces the Enrichr API call"]
    E[For each module > 3 genes]:::live
    F[Hypergeometric test against each GO term\nscipy.stats.hypergeom]:::live
    G[BH correction across all tested terms\nstatsmodels.multipletests]:::live
    H[Top term by adj. p-value\n< 10 ms per module  no network call]:::live
    E --> F --> G --> H
end

ONCE -.->|bundled with app data| BOOT
BOOT --> LIVE
```

### What changes in the code

Only the `annotate_multi_modules` callback in `coessentiality_feature_pc.py` needs updating — roughly a 10-line swap:

- **Remove:** `gp.enrichr(gene_list=..., gene_sets=..., organism=..., outdir=None)`
- **Replace with:** `_local_enrichr(gene_list, _GO_SETS, _GO_BACKGROUND)` using `scipy.stats.hypergeom` + `statsmodels` BH correction
- **Add at startup:** `_GO_SETS = _load_go_sets()` — loads `go_bp_2025.pkl` once, reused for every request

### Prototype vs production at a glance

| | Prototype (now) | Production (planned) |
|---|---|---|
| GO library source | Enrichr API, live per click | `.pkl` file, loaded at startup |
| Enrichment test | Remote (Enrichr server) | Local hypergeometric + BH |
| Speed | ~1–2 s per module | < 10 ms per module |
| Fails if Enrichr is down | Yes | No |
| Rate-limit risk on GCP | Yes | None |