Here's how the "public"."perturb_seq" table currently looks like:

dataset_id
perturbed_target_symbol
gene
log2foldchange
padj
basemean
max_ingested_at
adamson_2016_pilot
BHLHE40
TMSB10
0.27911572957264436
1.0606352255204349e-06
752.3006743141202
2025-09-30T13:37:59.073273Z
adamson_2016_pilot
BHLHE40
C14orf2
0.16943888923087036
0.04099857777763379
391.9803902568682
2025-09-30T13:37:59.073273Z
adamson_2016_pilot
BHLHE40
GMFG
0.8675978505270372
2.7854598007378423e-23
115.56715985676783
2025-09-30T13:37:59.073273Z
adamson_2016_pilot
BHLHE40
PFN1
0.1253097491506884
1.6602844397415236e-08
3108.9460780173954
2025-09-30T13:37:59.073273Z
adamson_2016_pilot
BHLHE40
ZNF296
3.1614764332620946
3.9957318157353094e-10
8.246859216074895
2025-09-30T13:37:59.073273Z

Based on the perturb_seq table, I need to compute two new tables:

# 1. "perturb_seq_summary_perturbation"
For each unique dataset_id + perturbed_target_symbol combination, compute:
* how many rows are present? → into n_total
* for how many of those, log2fold change is negative? → into n_down
* and for many of those, log2fold change is negative? → into n_up

# 2. "perturb_seq_summary_effect"
For each unique dataset_id + gene, compute:
* n_total, n_down, n_up using the same logic
* base_mean → use average value of all base_mean values for all rows for that given dataset_id + gene pair

Reply with the two relevant SQL blocks.
