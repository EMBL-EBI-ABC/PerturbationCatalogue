{{ config(materialized="table") }}

with
    mapped as (
        select
            source_data.perturbed_target_symbol as source_value,
            alias_resolution.n_ensembl_matches,
            alias_resolution.ensembl_gene_ids
        from {{ source("perturb_seq", "pertpy_gsea") }} as source_data
        left join {{ ref("stg_opentargets_alias_resolution") }} as alias_resolution
            on lower(trim(source_data.perturbed_target_symbol)) = alias_resolution.alias_key
    )

select
    'perturbed_target_symbol' as source_field,
    source_value,
    case
        when source_value is null or trim(source_value) = '' then 'blank'
        when n_ensembl_matches is null then 'unmapped'
        when n_ensembl_matches = 1 then 'mapped'
        else 'ambiguous'
    end as mapping_status,
    coalesce(n_ensembl_matches, 0) as n_ensembl_matches,
    ensembl_gene_ids,
    count(*) as row_count
from mapped
group by source_value, mapping_status, n_ensembl_matches, ensembl_gene_ids
