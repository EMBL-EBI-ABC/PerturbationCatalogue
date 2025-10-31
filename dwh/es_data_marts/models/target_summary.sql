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
            significant,
            perturbation_name,
            score_name,
            score_value,
            data_modality
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
    mave_data as (
        select
            perturbed_target_symbol,
            array_agg(
                struct(
                    score_value as score_value, perturbation_name as perturbation_name
                )
            ) as mave_scores
        from base
        where data_modality = 'MAVE' and score_name = 'score'
        group by perturbed_target_symbol
    ),
    agg_main as (
        select
            perturbed_target_symbol,
            count(distinct contrast_id) as n_experiments,
            countif(padj <= 0.05 and log2foldchange > 0) as n_sig_perturb_pairs_up,
            countif(padj <= 0.05 and log2foldchange < 0) as n_sig_perturb_pairs_down,
            countif(significant = 'true') as n_sig_crispr,
            countif(data_modality = 'MAVE' and score_name = 'score') as n_mave,
            array_agg(distinct data_modality ignore nulls) as data_modalities,
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
            ] as lfc_q75
        from base
        group by perturbed_target_symbol
    )
select
    m.perturbed_target_symbol,
    m.n_experiments,
    m.n_sig_perturb_pairs_up,
    m.n_sig_perturb_pairs_down,
    m.n_sig_crispr,
    m.n_mave,
    m.data_modalities,
    m.tissues_tested,
    m.cell_types_tested,
    m.cell_lines_tested,
    m.sex_tested,
    m.developmental_stages_tested,
    m.diseases_tested,
    m.lfc_q25,
    m.lfc_median,
    m.lfc_q75,
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
    array(
        select as struct mv.perturbation_name, mv.score_value
        from unnest(d.mave_scores) as mv
        order by mv.score_value desc
        limit 10
    ) as mave_scores_up,
    array(
        select as struct mv.perturbation_name, mv.score_value
        from unnest(d.mave_scores) as mv
        order by mv.score_value asc
        limit 10
    ) as mave_scores_down
from agg_main m
left join sig_gene_pairs g using (perturbed_target_symbol)
left join mave_data d using (perturbed_target_symbol)
