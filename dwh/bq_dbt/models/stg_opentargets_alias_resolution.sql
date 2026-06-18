{{ config(materialized="table") }}

with
    alias_rows as (
        select
            lower(trim(alias)) as alias_key,
            trim(alias) as alias,
            ensembl_gene_id
        from {{ ref("stg_opentargets_targets") }}, unnest(exact_aliases) as alias
        where alias is not null and trim(alias) != ''

        union distinct

        select
            lower(trim(ensembl_gene_id)) as alias_key,
            ensembl_gene_id as alias,
            ensembl_gene_id
        from {{ ref("stg_opentargets_targets") }}
        where ensembl_gene_id is not null and trim(ensembl_gene_id) != ''
    ),

    grouped as (
        select
            alias_key,
            array_agg(distinct alias order by alias) as aliases,
            array_agg(distinct ensembl_gene_id order by ensembl_gene_id) as ensembl_gene_ids,
            count(distinct ensembl_gene_id) as n_ensembl_matches
        from alias_rows
        group by alias_key
    )

select
    alias_key,
    aliases,
    ensembl_gene_ids,
    n_ensembl_matches,
    if(n_ensembl_matches = 1, ensembl_gene_ids[offset(0)], null) as ensembl_gene_id
from grouped
