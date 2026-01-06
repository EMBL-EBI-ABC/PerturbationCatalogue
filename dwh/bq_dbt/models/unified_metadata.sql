{{
    config(
        materialized="table",
        partition_by={
            "field": "max_ingested_at",
            "data_type": "timestamp",
            "granularity": "day",
        },
        cluster_by=["dataset_id", "sample_id", "perturbed_target_symbol"],
    )
}}

with

    -- Base metadata from sources
    crispr_base as (
        select
            *,
            ingested_at as max_ingested_at
        from {{ source('crispr', 'metadata') }}

    ),
    mave_base as (
        select
            *,
            ingested_at as max_ingested_at
        from {{ source('mave', 'metadata') }}

    ),
    ps_base as (
        select
            *,
            ingested_at as max_ingested_at
        from (
            select distinct * except (sample_id)
            from {{ source('perturb_seq', 'metadata') }}
            where
                perturbed_target_symbol not like 'control%'
                and perturbed_target_symbol not like '%None%'

        )
    ),
    
    -- Normalize columns for UNION ALL BY NAME
    -- We must ensure ALL branches have the SAME columns to satisfy the specific BigQuery constraint reported.
    
    crispr as (
        select
            * except (ingested_at)
        from crispr_base
    ),
    mave as (
        select
            * except (ingested_at)
        from mave_base
    ),
    ps as (
        select
            * except (
                ingested_at,
                significant,
                significance_criteria,
                number_of_perturbed_targets,
                number_of_perturbed_samples,
                library_total_grnas
            ),
            -- Add sample_id which is missing in PS
            cast(null as string) as sample_id,
            -- Override these to match previous logic (nulling and casting to string)
            cast(null as string) as significant,
            cast(null as string) as significance_criteria,
            cast(number_of_perturbed_targets as string) as number_of_perturbed_targets,
            cast(number_of_perturbed_samples as string) as number_of_perturbed_samples,
            cast(library_total_grnas as string) as library_total_grnas
        from ps_base
    ),
    
    unified as (
        select * from crispr
        full outer union all by name
        select * from mave
        full outer union all by name
        select * from ps
    )
    
select * from unified
