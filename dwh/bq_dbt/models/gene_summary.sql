{{ config(materialized="table") }}

with
    base as (
        select
            gene,
            to_hex(
                sha256(concat(dataset_id, '|', coalesce(perturbed_target_symbol, '')))
            ) as contrast_id,
            dataset_id,
            perturbed_target_symbol,
            log2foldchange,
            padj,
            tissue_label,
            cell_type_label,
            cell_line_label,
            sex_label,
            developmental_stage_label,
            disease_label
        from {{ ref("unified_metadata_data") }}
    ),
    sig_contrasts as (
        select
            gene,
            contrast_id,
            dataset_id,
            perturbed_target_symbol,
            log2foldchange as l2fc,
            padj,
        from base
        where padj <= 0.05
    ),
    agg_main as (
        select
            gene,
            -- 1) in how many contrasts this gene was tested
            count(distinct contrast_id) as n_contrasts_tested,
            -- 2) significant contrasts (and up/down)
            count(distinct if(padj <= 0.05, contrast_id, null)) as n_contrasts_sig,
            count(
                distinct if(padj <= 0.05 and log2foldchange > 0, contrast_id, null)
            ) as n_contrasts_up_sig,
            count(
                distinct if(padj <= 0.05 and log2foldchange < 0, contrast_id, null)
            ) as n_contrasts_down_sig,
            -- 3) log2FC quantiles for significant results only
            approx_quantiles(case when padj <= 0.05 then log2foldchange end, 5)[
                offset (2)
            ] as lfc_q25,
            approx_quantiles(case when padj <= 0.05 then log2foldchange end, 5)[
                offset (3)
            ] as lfc_median,
            approx_quantiles(case when padj <= 0.05 then log2foldchange end, 5)[
                offset (4)
            ] as lfc_q75,
            -- 4) tissues where tested
            array_agg(
                distinct tissue_label ignore nulls order by tissue_label
            ) as tissues_tested,
            array_agg(
                distinct cell_type_label ignore nulls order by cell_type_label
            ) as cell_types_tested,
            array_agg(
                distinct cell_line_label ignore nulls order by cell_line_label
            ) as cell_lines_tested,
            array_agg(distinct sex_label ignore nulls order by sex_label) as sex_tested,
            array_agg(
                distinct developmental_stage_label ignore nulls
                order by developmental_stage_label
            ) as developmental_stages_tested,
            array_agg(
                distinct disease_label ignore nulls order by disease_label
            ) as diseases_tested,
        from base
        group by gene
    ),
    agg_arrays as (
        select
            gene,
            array_agg(
                struct(
                    perturbed_target_symbol as perturbed_target_symbol,
                    padj as padj,
                    l2fc as l2fc
                )
                order by padj asc, abs(l2fc) desc
            ) as sig_contrasts_with_lfc
        from sig_contrasts
        group by gene
    )

select
    m.gene,
    m.n_contrasts_tested,
    m.n_contrasts_sig,
    m.n_contrasts_up_sig,
    m.n_contrasts_down_sig,
    m.lfc_q25,
    m.lfc_median,
    m.lfc_q75,
    m.tissues_tested,
    m.cell_types_tested,
    m.cell_lines_tested,
    m.sex_tested,
    m.developmental_stages_tested,
    m.diseases_tested,
    array(
        select as struct c.perturbed_target_symbol, c.padj, c.l2fc
        from unnest(a.sig_contrasts_with_lfc) as c
        where c.l2fc > 0
        limit 10
    ) as contrasts_sig_details_up,
    array(
        select as struct c.perturbed_target_symbol, c.padj, c.l2fc
        from unnest(a.sig_contrasts_with_lfc) as c
        where c.l2fc < 0
        limit 10
    ) as contrasts_sig_details_down,
from agg_main m
left join agg_arrays a using (gene)
