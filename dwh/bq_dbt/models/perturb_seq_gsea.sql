{{ config(materialized="table") }}

with
    source_data as (
        select
            dataset_id,
            term,
            perturbed_target_symbol,
            es,
            nes,
            pval,
            sidak,
            fdr,
            geneset_size,
            leading_edge,
            cell_type,
            ingested_at
        from {{ source("perturb_seq", "pertpy_gsea") }}
    )

select
    source_data.dataset_id,
    source_data.term,
    target_map.ensembl_gene_id as perturbed_target_ensg,
    source_data.es,
    source_data.nes,
    source_data.pval,
    source_data.sidak,
    source_data.fdr,
    source_data.geneset_size,
    source_data.leading_edge,
    source_data.cell_type,
    source_data.ingested_at as max_ingested_at
from source_data
left join {{ ref("stg_opentargets_alias_resolution") }} as target_map
    on lower(trim(source_data.perturbed_target_symbol)) = target_map.alias_key
