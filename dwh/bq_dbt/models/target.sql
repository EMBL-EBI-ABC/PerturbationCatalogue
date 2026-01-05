{{ config(materialized="table") }}


with
    meta as (
        select * from {{ ref("unified_metadata") }}
    ),
    crispr_joined as (
        select
            m.*,
            d.score_name,
            cast(d.score_value as float64) as score_value,
            cast(null as float64) as log2foldchange,
            cast(null as float64) as padj
        from meta m
        join {{ source('crispr', 'data') }} d
            on m.dataset_id = d.dataset_id and m.sample_id = d.sample_id
        where m.data_modality = 'CRISPR'
    ),
    mave_joined as (
        select
            m.*,
            d.score_name,
            cast(d.score_value as float64) as score_value,
            cast(null as float64) as log2foldchange,
            cast(null as float64) as padj
        from meta m
        join {{ source('mave', 'data') }} d
            on m.dataset_id = d.dataset_id and m.sample_id = d.sample_id
        where m.data_modality = 'MAVE'
    ),
    ps_joined as (
        select
            m.*,
            cast(null as string) as score_name,
            cast(null as float64) as score_value,
            d.log2foldchange,
            d.padj
        from meta m
        join {{ source('perturb_seq', 'data') }} d
            on m.dataset_id = d.dataset_id and m.perturbed_target_symbol = d.perturbation
        where m.data_modality = 'Perturb-seq'
    ),
    base_unioned as (
        select * from crispr_joined
        union all by name
        select * from mave_joined
        union all by name
        select * from ps_joined
    ),
    base as (
        select
            to_hex(
                sha256(concat(dataset_id, '|', coalesce(perturbed_target_symbol, '')))
            ) as contrast_id,
            perturbed_target_symbol,
            log2foldchange,
            padj,
            tissue_label,
            cell_type_label,
            cell_line_label,
            sex_label,
            developmental_stage_label,
            disease_label,
            significant,
            data_modality,
            license_label,
            score_name
        from base_unioned
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
            array_agg(distinct license_label ignore nulls) as license
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
    m.license,
    m.data_modalities,
    m.tissues_tested,
    m.cell_types_tested,
    m.cell_lines_tested,
    m.sex_tested,
    m.developmental_stages_tested,
    m.diseases_tested,
from agg_main m
