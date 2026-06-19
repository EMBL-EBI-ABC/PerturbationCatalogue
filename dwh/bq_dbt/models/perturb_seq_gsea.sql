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
    ),

    distinct_targets as (
        select distinct perturbed_target_symbol
        from source_data
        where perturbed_target_symbol is not null
    ),

    split_targets as (
        select
            targets.perturbed_target_symbol,
            offset_pos,
            trim(symbol_part) as symbol_part,
            lower(trim(symbol_part)) as alias_key
        from
            distinct_targets as targets,
            unnest(split(targets.perturbed_target_symbol, "|")) as symbol_part
            with offset as offset_pos
    ),

    mapped_target_parts as (
        select
            split_targets.perturbed_target_symbol,
            split_targets.offset_pos,
            case
                when starts_with(split_targets.alias_key, "control_")
                then split_targets.alias_key
                else alias_resolution.ensembl_gene_id
            end as mapped_symbol
        from split_targets
        left join {{ ref("stg_opentargets_alias_resolution") }} as alias_resolution
            on split_targets.alias_key = alias_resolution.alias_key
    ),

    mapped_targets as (
        select
            perturbed_target_symbol,
            string_agg(mapped_symbol, "|" order by offset_pos) as perturbed_target_ensg
        from mapped_target_parts
        group by perturbed_target_symbol
        having
            count(*) > 0
            and countif(mapped_symbol is null) = 0
            and countif(
                not (
                    starts_with(mapped_symbol, "ENSG")
                    or starts_with(mapped_symbol, "control_")
                )
            )
            = 0
    )

select
    source_data.dataset_id,
    source_data.term,
    mapped_targets.perturbed_target_ensg,
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
inner join mapped_targets
    on source_data.perturbed_target_symbol = mapped_targets.perturbed_target_symbol
