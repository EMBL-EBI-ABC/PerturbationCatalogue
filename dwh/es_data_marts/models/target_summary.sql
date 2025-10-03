{{ config(materialized="table") }}


with
    base as (
        select
            to_hex(
                sha256(concat(dataset_id, '|', coalesce(perturbed_target_symbol, '')))
            ) as contrast_id,
            dataset_id,
            perturbed_target_symbol,
            gene,
            log2foldchange,
            padj,
            tissue_label,
            cell_type_label,
            cell_line_label,
            sex_label,
            developmental_stage_label,
            disease_label,
            significant
        from {{ ref("unified_metadata_data") }}
    ),
    sig_gene_pairs as (
        select
            perturbed_target_symbol,
            array_agg(
                struct(gene as gene, padj as padj, log2foldchange as log2foldchange)
                order by padj asc, abs(log2foldchange) desc
            ) as sig_contrasts_with_lfc
        from base
        where padj <= 0.05
        group by perturbed_target_symbol
    ),
    agg_main as (
        select
            perturbed_target_symbol,
            count(distinct contrast_id) as n_contrasts,
            countif(padj <= 0.05 or significant = "true") as n_sig_contrasts,
            countif(padj <= 0.05 and log2foldchange > 0) as n_sig_perturb_pairs_up,
            countif(padj <= 0.05 and log2foldchange < 0) as n_sig_perturb_pairs_down,
            array_agg(distinct tissue_label ignore nulls) as tissues_tested,
            array_agg(distinct cell_type_label ignore nulls) as cell_types_tested,
            array_agg(distinct cell_line_label ignore nulls) as cell_lines_tested,
            array_agg(distinct sex_label ignore nulls) as sex_tested,
            array_agg(
                distinct developmental_stage_label ignore nulls
            ) as developmental_stages_tested,
            array_agg(distinct disease_label ignore nulls) as diseases_tested,
            approx_quantiles(case when padj <= 0.05 then log2foldchange end, 5)[
                offset (2)
            ] as lfc_q25,
            approx_quantiles(case when padj <= 0.05 then log2foldchange end, 5)[
                offset (3)
            ] as lfc_median,
            approx_quantiles(case when padj <= 0.05 then log2foldchange end, 5)[
                offset (4)
            ] as lfc_q75,
            countif(significant = 'true') as n_sig_crispr_contrasts
        from base
        group by perturbed_target_symbol
    )
select
    m.perturbed_target_symbol,
    m.n_contrasts,
    m.n_sig_contrasts,
    n_sig_perturb_pairs_up,
    n_sig_perturb_pairs_down,
    m.tissues_tested,
    m.cell_types_tested,
    m.cell_lines_tested,
    m.sex_tested,
    m.developmental_stages_tested,
    m.diseases_tested,
    m.lfc_q25,
    m.lfc_median,
    m.lfc_q75,
    m.n_sig_crispr_contrasts,
    array(
        select as struct c.gene, c.padj, c.log2foldchange
        from unnest(g.sig_contrasts_with_lfc) as c
        where c.log2foldchange > 0
        order by padj asc, abs(log2foldchange) desc
        limit 10
    ) as sig_gene_pairs_up,
    array(
        select as struct c.gene, c.padj, c.log2foldchange
        from unnest(g.sig_contrasts_with_lfc) as c
        where c.log2foldchange < 0
        order by padj asc, abs(log2foldchange) desc
        limit 10
    ) as sig_gene_pairs_down,
from agg_main m
left join sig_gene_pairs g using (perturbed_target_symbol)
