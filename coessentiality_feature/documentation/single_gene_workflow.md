# Single-Gene Explorer — Workflow

The single-gene explorer lets users search for one gene and see all its co-essential partners — ranked by significance, coloured by correlation, and laid out as a network graph.

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
B --> F[all_genes list\npopulates the gene search dropdown]:::data
```

---

## 2. Single-gene view — full workflow

The entire view is driven by a **single callback** that fires whenever the user picks a gene or changes the FDR threshold.

```mermaid
flowchart TD

classDef user  fill:#4C72B0,color:#fff,stroke:none,rx:16
classDef fast  fill:#2E7D52,color:#fff,stroke:none
classDef out   fill:#E8F5EE,color:#212121,stroke:#2E7D52

A([Select gene from dropdown]):::user
B([Set FDR threshold\n5% or 10%]):::user

subgraph CB["Callback — update_single  runs on gene selection or FDR change"]
    direction TB
    C[Filter df_all\nkeep rows where pvalue_adj ≤ FDR]:::fast
    D[Find all pairs containing the query gene\nbuild partner list sorted by adj. p-value]:::fast
    C --> D

    D --> E[Summary text\nGENE has N co-essential partners at FDR ≤ X%]:::out
    D --> F[Bar chart\ntop 20 partners ranked by −log₁₀ adj. p-value]:::out
    D --> G[Partner table\nPARTNER GENE · P-VALUE · ADJ. P-VALUE · CORRELATION]:::out
    D --> H[Build Cytoscape elements\nquery gene node = orange\npartner nodes coloured by corr_genes\nred = positive · green = negative correlation\nedge weight = −log₁₀ adj. p-value]:::fast
    H --> I[/Network graph rendered\ntop 20 partners + cross-edges among them/]:::out
end

A --> C
B --> C
```

---

## 3. How node colours are computed

Partner nodes are coloured on a **red → grey → green** gradient based on their `corr_genes` value — the sign of the co-essentiality relationship with the query gene.

```mermaid
flowchart LR

classDef pos fill:#E45756,color:#fff,stroke:none
classDef neu fill:#DCDCDC,color:#333,stroke:#aaa
classDef neg fill:#54A24B,color:#fff,stroke:none

A[corr_genes = +1]:::pos --> B[Red  both genes are\nco-essential together]:::pos
C[corr_genes = 0]:::neu --> D[Grey  no directional\ncorrelation]:::neu
E[corr_genes = −1]:::neg --> F[Green  genes are\nmutually essential in opposite contexts]:::neg
```

The query gene itself is always shown in **orange** to distinguish it from its partners.

---

## 4. Notes

- All computation is **local and instantaneous** — no external API calls are made in the single-gene view.
- The FDR dropdown is shared across both tabs; changing it updates the single-gene view and the gene-list view simultaneously.
- The bar chart and table both show all partners at the chosen FDR; the Cytoscape graph is capped at the **top 20** to keep the layout readable.
- For the gene-list view (Tab 2) and its GO:BP annotation workflow, see `gene_list_workflow.md`.