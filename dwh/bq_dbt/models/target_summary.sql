{{ config(materialized="table") }}

with
    metadata_targets as (
        select
            metadata.* except (perturbed_target_ensg),
            trim(target_ensg) as ensembl_gene_id
        from
            {{ ref("unified_metadata") }} as metadata,
            unnest(split(cast(metadata.perturbed_target_ensg as string), "|")) as target_ensg
        where
            metadata.perturbed_target_ensg is not null
            and trim(target_ensg) != ""
            and starts_with(trim(target_ensg), "ENSG")
    ),

    agg_meta as (
        select
            ensembl_gene_id,
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
        from metadata_targets
        group by ensembl_gene_id
    ),

    agg_crispr as (
        select
            perturbed_target_ensg as ensembl_gene_id,
            count(distinct dataset_id) as n_crispr,
            countif(significant = 'True') as n_sig_crispr
        from {{ ref("crispr_data") }}
        where
            perturbed_target_ensg is not null
            and starts_with(perturbed_target_ensg, "ENSG")
        group by ensembl_gene_id
    ),

    agg_mave as (
        select
            perturbed_target_ensg as ensembl_gene_id,
            count(distinct dataset_id) as n_mave
        from {{ ref("mave_data") }}
        where
            perturbed_target_ensg is not null
            and starts_with(perturbed_target_ensg, "ENSG")
        group by ensembl_gene_id
    ),

    agg_ps as (
        select
            perturbed_target_ensg as ensembl_gene_id,
            count(distinct dataset_id) as n_perturb_seq,
            countif(padj <= 0.05 and log2foldchange > 0) as n_sig_perturb_pairs_up,
            countif(padj <= 0.05 and log2foldchange < 0) as n_sig_perturb_pairs_down
        from {{ source("perturb_seq", "pertpy_dea") }}
        where
            perturbed_target_ensg is not null
            and starts_with(perturbed_target_ensg, "ENSG")
        group by ensembl_gene_id
    ),

    gsea_ranked as (
        select
            perturbed_target_ensg as ensembl_gene_id,
            term,
            sidak,
            row_number() over (
                partition by perturbed_target_ensg order by sidak asc
            ) as rn
        from {{ source("perturb_seq", "pertpy_gsea") }}
        where
            perturbed_target_ensg is not null
            and starts_with(perturbed_target_ensg, "ENSG")
            and sidak <= 0.05
    ),

    agg_gsea as (
        select
            ensembl_gene_id,
            array_agg(term order by sidak asc) as top_gsea_terms
        from gsea_ranked
        where rn <= 5
        group by ensembl_gene_id
    ),

    effect_targets as (
        select distinct effect_gene_ensg as ensembl_gene_id
        from {{ source("perturb_seq", "pertpy_dea") }}
        where
            effect_gene_ensg is not null
            and starts_with(effect_gene_ensg, "ENSG")
    ),

    targets as (
        select ensembl_gene_id
        from agg_meta
        union distinct
        select ensembl_gene_id
        from agg_mave
        union distinct
        select ensembl_gene_id
        from agg_ps
        union distinct
        select ensembl_gene_id
        from agg_crispr
        union distinct
        select ensembl_gene_id
        from agg_gsea
        union distinct
        select ensembl_gene_id
        from effect_targets
    )

select
    targets.ensembl_gene_id,
    opentargets.approved_symbol,
    opentargets.approved_name,
    opentargets.exact_aliases,
    opentargets.search_keywords,
    coalesce(perturb_seq.n_perturb_seq, 0) as n_perturb_seq,
    coalesce(perturb_seq.n_sig_perturb_pairs_up, 0) as n_sig_perturb_pairs_up,
    coalesce(perturb_seq.n_sig_perturb_pairs_down, 0) as n_sig_perturb_pairs_down,
    coalesce(crispr.n_crispr, 0) as n_crispr,
    coalesce(crispr.n_sig_crispr, 0) as n_sig_crispr,
    coalesce(mave.n_mave, 0) as n_mave,
    gsea.top_gsea_terms,
    metadata.license,
    metadata.data_modalities,
    metadata.tissues_tested,
    metadata.cell_types_tested,
    metadata.cell_lines_tested,
    metadata.sex_tested,
    metadata.developmental_stages_tested,
    metadata.diseases_tested
from targets
inner join {{ ref("stg_opentargets_targets") }} as opentargets
    on targets.ensembl_gene_id = opentargets.ensembl_gene_id
left join agg_meta as metadata on targets.ensembl_gene_id = metadata.ensembl_gene_id
left join agg_mave as mave on targets.ensembl_gene_id = mave.ensembl_gene_id
left join agg_ps as perturb_seq on targets.ensembl_gene_id = perturb_seq.ensembl_gene_id
left join agg_crispr as crispr on targets.ensembl_gene_id = crispr.ensembl_gene_id
left join agg_gsea as gsea on targets.ensembl_gene_id = gsea.ensembl_gene_id
