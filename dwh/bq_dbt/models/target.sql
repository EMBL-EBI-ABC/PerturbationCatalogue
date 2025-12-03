{{ config(materialized="table") }}


with
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
        from {{ ref("unified_metadata_data") }}
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
