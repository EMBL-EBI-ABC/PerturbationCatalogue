# Replogle unified pipeline vs production API

Compared on 18 September 2026 against the live production backend at
`https://perturbation-catalogue-be-asmuzum42q-nw.a.run.app` using the completed
cluster products from the unified run. Production DEA queries used the sampled
cluster target Ensembl IDs, avoiding the API's alias-only symbol resolution
problem (for example, cluster `SELK` is production `SELENOK`).

## Sample

- 249 effect genes: 2% of the 12,420-gene analysis space, selected by a
  deterministic SHA-256 ordering.
- 12 perturbation targets, selected by the same deterministic ordering.
- 2,988 cluster DEA rows and 576 cluster GSEA rows.
- 5,000-row DEA API pages; all 12 target queries returned 12,420 rows and all
  12 GSEA queries returned 48 terms.

## Key coverage

Rows are joined by `(perturbed_target_ensg, effect_gene_ensg)` for DEA and
`(perturbed_target_ensg, term)` for GSEA. Numeric metrics below are calculated
once per unique matched key.

| Product | Cluster rows | Production rows after sample filter | Matched keys | Cluster → production | Production → cluster |
|---|---:|---:|---:|---:|---:|
| DEA | 2,988 | 3,033 | 2,560 | 85.676% | 84.405% |
| GSEA | 576 | 576 | 576 | 100.000% | 100.000% |

The DEA sample contains 473 duplicate production keys after filtering, and 428
cluster keys have no exact production Ensembl-key counterpart. These are
identifier/row-coverage differences, not numerical differences on matched
rows; they should be investigated separately if complete feature-universe
identity is required.

## Numerical agreement on matched rows

| Product / metric | Pearson | Spearman | Median absolute error | Maximum absolute error | Sign agreement |
|---|---:|---:|---:|---:|---:|
| DEA log2 fold change (n=2,560) | 0.9999999999 | 0.9999998977 | 5.08e-05 | 1.25e-03 | 100% |
| DEA Wilcoxon score (n=2,560) | 0.9999999660 | 0.9999997141 | 6.15e-05 | 5.44e-04 | 100% |
| GSEA ES (n=576) | 0.9999999421 | 0.9999973627 | 3.39e-05 | 7.12e-04 | 100% |
| GSEA NES (n=576) | 0.9999566205 | 0.9997329390 | 0.00610 | 0.0399 | 100% |

For significance values, `-log10(padj)` Pearson correlation was
0.9999999984 (median absolute error 1.27e-05); GSEA `-log10(pval)` and
`-log10(fdr)` Pearson correlations were 0.99015 and 0.99466 respectively.
GSEA leading-edge overlap had median Jaccard similarity **1.000** and mean
similarity **0.99905**.

The matched DEA and GSEA values are therefore effectively identical up to
small floating-point/recalculation differences. The only material caveat in
this sample is DEA key coverage, which is reported explicitly above.
