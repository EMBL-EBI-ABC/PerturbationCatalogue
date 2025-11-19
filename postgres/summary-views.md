```sql
CREATE MATERIALIZED VIEW perturb_seq_summary_perturbation AS
SELECT
    dataset_id,
    perturbed_target_symbol,
    COUNT(*) AS n_total,
    COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
    COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up
FROM public.perturb_seq_2
WHERE padj <= 0.05
GROUP BY dataset_id, perturbed_target_symbol;
```

```sql
CREATE MATERIALIZED VIEW perturb_seq_summary_effect AS
SELECT
    dataset_id,
    gene,
    COUNT(*) AS n_total,
    COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
    COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up,
    AVG(basemean) AS base_mean
FROM public.perturb_seq_2
WHERE padj <= 0.05
GROUP BY dataset_id, gene;
```