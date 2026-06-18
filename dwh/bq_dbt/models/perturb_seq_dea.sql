{{ config(materialized="table") }}

with
    source_data as (
        select
            dataset_id,
            perturbed_target_symbol,
            gene as effect_gene_symbol,
            padj,
            log2FoldChange,
            score_name,
            score_value,
            cell_type,
            ingested_at
        from {{ source("perturb_seq", "pertpy_dea") }}
    )

select
    source_data.dataset_id,
    target_map.ensembl_gene_id as perturbed_target_ensg,
    effect_map.ensembl_gene_id as effect_gene_ensg,
    source_data.padj,
    source_data.log2FoldChange as log2foldchange,
    source_data.score_name,
    source_data.score_value,
    source_data.cell_type,
    source_data.ingested_at as max_ingested_at
from source_data
left join {{ ref("stg_opentargets_alias_resolution") }} as target_map
    on lower(trim(source_data.perturbed_target_symbol)) = target_map.alias_key
left join {{ ref("stg_opentargets_alias_resolution") }} as effect_map
    on lower(trim(source_data.effect_gene_symbol)) = effect_map.alias_key
