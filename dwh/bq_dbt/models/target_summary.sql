{{ config(materialized="table") }}

with
    agg_meta as (
        select
            perturbed_target_symbol,
            count(distinct dataset_id) as n_experiments,
            countif(significant = 'true') as n_sig_crispr,
            array_agg(distinct data_modality ignore nulls) as data_modalities,
            array_agg(distinct tissue_label ignore nulls) as tissues_tested,
            array_agg(distinct cell_type_label ignore nulls) as cell_types_tested,
            array_agg(distinct cell_line_label ignore nulls) as cell_lines_tested,
            array_agg(distinct sex_label ignore nulls) as sex_tested,
            array_agg(distinct developmental_stage_label ignore nulls) as developmental_stages_tested,
            array_agg(distinct disease_label ignore nulls) as diseases_tested,
            array_agg(distinct license_label ignore nulls) as license
        from {{ ref("unified_metadata") }}
        group by perturbed_target_symbol
    ),

    agg_mave as (
        select
            perturbed_target_symbol,
            countif(score_name = 'score') as n_mave
        from {{ ref("mave_data") }}
        group by perturbed_target_symbol
    ),

    agg_ps as (
        select
            perturbed_target_symbol,
            countif(padj <= 0.05 and log2foldchange > 0) as n_sig_perturb_pairs_up,
            countif(padj <= 0.05 and log2foldchange < 0) as n_sig_perturb_pairs_down
        from {{ ref("perturb_seq_data") }}
        group by perturbed_target_symbol
    ),

    -- Supersets of symbols to ensure we don't miss any target that might exist only in data (unlikely but safe)
    symbols as (
        select perturbed_target_symbol from agg_meta
        union distinct
        select perturbed_target_symbol from agg_mave
        union distinct
        select perturbed_target_symbol from agg_ps
    )

select
    s.perturbed_target_symbol,
    coalesce(m.n_experiments, 0) as n_experiments,
    coalesce(p.n_sig_perturb_pairs_up, 0) as n_sig_perturb_pairs_up,
    coalesce(p.n_sig_perturb_pairs_down, 0) as n_sig_perturb_pairs_down,
    coalesce(m.n_sig_crispr, 0) as n_sig_crispr,
    coalesce(v.n_mave, 0) as n_mave,
    m.license,
    m.data_modalities,
    m.tissues_tested,
    m.cell_types_tested,
    m.cell_lines_tested,
    m.sex_tested,
    m.developmental_stages_tested,
    m.diseases_tested
from symbols s
left join agg_meta m on s.perturbed_target_symbol = m.perturbed_target_symbol
left join agg_mave v on s.perturbed_target_symbol = v.perturbed_target_symbol
left join agg_ps p on s.perturbed_target_symbol = p.perturbed_target_symbol
