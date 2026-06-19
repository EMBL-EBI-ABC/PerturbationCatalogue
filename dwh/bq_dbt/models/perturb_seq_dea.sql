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
    mapped_targets.perturbed_target_ensg,
    effect_map.ensembl_gene_id as effect_gene_ensg,
    source_data.padj,
    source_data.log2FoldChange as log2foldchange,
    source_data.score_name,
    source_data.score_value,
    source_data.cell_type,
    source_data.ingested_at as max_ingested_at
from source_data
inner join mapped_targets
    on source_data.perturbed_target_symbol = mapped_targets.perturbed_target_symbol
left join {{ ref("stg_opentargets_alias_resolution") }} as effect_map
    on lower(trim(source_data.effect_gene_symbol)) = effect_map.alias_key
