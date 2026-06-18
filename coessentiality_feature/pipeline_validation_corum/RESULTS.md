# CORUM Validation Results — GLS Co-essentiality Pipeline

Validation for [PR #365](https://github.com/EMBL-EBI-ABC/PerturbationCatalogue/pull/365)
Comment 4 ("we need *some* sort of benchmark which takes pipeline results as
input and produces a numerical metric"). Replicates the spirit of Wainberg et
al.'s Figure 2 / Extended Data Figure 4 ("GLS improves recall of known
functional interactions in co-essential gene pairs"), scoped to a single gold
standard (CORUM) and a single method (GLS) — see `REVIEW.md` Comment 4 for the
full scoping discussion.

## What this answers

> "After we run this pipeline, does it produce sensible results?"

## Method

1. **Gold standard:** Enrichr's `CORUM` gene-set library (1,658 human protein
   complexes, fetched via `gseapy.get_library()` — same mechanism as the
   GO:BP fix in Comment 10). Every pair of genes that co-occur in the same
   complex is treated as a "true" co-essential-like pair.
2. **Ranking:** for every gene with at least one CORUM complex-mate present
   in our 17,087-gene panel, rank all other genes by GLS p-value (ascending)
   and take its top-N partners, for N = 1 to 10.
3. **Metric ("enrichment"):**

   ```
   enrichment(N) = (% of top-N pairs that are true CORUM pairs)
                 / (% of all possible pairs that are true CORUM pairs)
   ```

   `1.0` = GLS ranking is no better than random. The paper reports
   several-fold enrichment for GLS at low N — this is the same metric, just
   scoped to one database instead of four (CORUM, hu.MAP, STRING, DoRothEA).

Data used: `depmap_26Q1_GLS_p.npy` / `depmap_26Q1_genes.txt` (current
production pipeline output, 17,087 genes, post-Comment-5 fix). No new pipeline
run was needed — this reuses the existing GLS output directly.

## Results

| N (top-N partners) | Observed hit rate | Enrichment vs. chance |
|---:|---:|---:|
| 1 | 26.93% | **1,657.6x** |
| 2 | 23.49% | 1,446.0x |
| 3 | 20.93% | 1,288.3x |
| 4 | 18.63% | 1,146.8x |
| 5 | 16.92% | 1,041.4x |
| 6 | 15.51% | 955.0x |
| 7 | 14.35% | 883.3x |
| 8 | 13.36% | 822.7x |
| 9 | 12.47% | 767.8x |
| 10 | 11.74% | 722.6x |

- **Genes evaluated:** 2,284 (genes with ≥1 CORUM complex-mate inside our panel)
- **Background rate:** 0.0162% (23,713 true CORUM pairs out of 145,974,241
  total possible pairs among our 17,087 genes)

## How to read this

At N=1, GLS's single most-significant predicted partner for a gene is its
*actual* known complex-mate **26.9% of the time** — versus a 0.016% chance of
that happening by pure luck. That's roughly **1,658x enrichment over random
chance**.

The smooth decline from N=1 to N=10 is the expected, sane pattern: your
single best-ranked prediction should be more reliable than your 10th-best
one, on average. A flat line near `1.0x` across all N would have meant GLS
isn't finding anything real; this is the opposite of that.

**Why the numbers are this large:** protein complex subunits are one of the
strongest known co-essentiality signals in the literature — losing any one
subunit of an essential complex tends to break the whole complex's function,
so its members are almost always strongly co-essential together. A few
-hundred-to-few-thousand-fold enrichment for this specific gold standard is
consistent with that, not a red flag.

## Caveats / what this does *not* show

- **Single gold standard only.** This only checks against CORUM (protein
  complexes). It doesn't yet cover STRING, hu.MAP, or DoRothEA, and doesn't
  yet test the looser "functional interaction" signal those capture beyond
  physical complex membership.
- **GLS only, no comparison method.** This doesn't show GLS is *better* than
  Pearson correlation or co-expression — only that GLS itself is sensible.
  The Comment 1 bias-correction question (does GLS need OR-gene PCA
  correction with Chronos-corrected data?) needs the with/without-correction
  comparison from Extended Data Fig 4, not this script alone.
- **Not a regression-testable benchmark yet.** This is a one-off validation
  script (`validate_corum_enrichment.py`), not wired into CI or the
  production pipeline. Re-run manually if you change something pipeline-side
  and want to check it didn't break anything.

## Files

- `validate_corum_enrichment.py` — the script that produced this
- `CORUM.gmt` — cached gold-standard data (pinned snapshot, same
  not-auto-updated philosophy as `GO_Biological_Process_2025.gmt`)
