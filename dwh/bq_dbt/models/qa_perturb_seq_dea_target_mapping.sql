{{ config(materialized="table") }}

with
    source_symbols as (
        select 'perturbed_target_symbol' as source_field, perturbed_target_symbol as source_value
        from {{ source("perturb_seq", "pertpy_dea") }}

        union all

        select 'effect_gene_symbol' as source_field, gene as source_value
        from {{ source("perturb_seq", "pertpy_dea") }}
    ),

    mapped as (
        select
            source_symbols.source_field,
            source_symbols.source_value,
            alias_resolution.n_ensembl_matches,
            alias_resolution.ensembl_gene_ids
        from source_symbols
        left join {{ ref("stg_opentargets_alias_resolution") }} as alias_resolution
            on lower(trim(source_symbols.source_value)) = alias_resolution.alias_key
    )

select
    source_field,
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
group by source_field, source_value, mapping_status, n_ensembl_matches, ensembl_gene_ids
