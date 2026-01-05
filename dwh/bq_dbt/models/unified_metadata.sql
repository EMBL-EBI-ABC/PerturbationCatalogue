{{
    config(
        materialized="incremental",
        incremental_strategy="insert_overwrite",
        partition_by={
            "field": "max_ingested_at",
            "data_type": "timestamp",
            "granularity": "day",
        },
        cluster_by=["dataset_id", "sample_id", "perturbed_target_symbol"],
    )
}}

with
    latest_loaded_partition as (
        select parse_date('%Y%m%d', max(partition_id)) as pdate
        from `{{ this.database }}`.`{{ this.schema }}.INFORMATION_SCHEMA.PARTITIONS`
        where
            table_name = '{{ this.identifier }}'
            and partition_id not in ('__NULL__', '__UNPARTITIONED__')
    ),
    -- Base metadata from sources, ensuring common ingested_at column name
    crispr_base as (
        select
            *,
            ingested_at as max_ingested_at
        from {{ source('crispr', 'metadata') }}
        {% if is_incremental() %}
            where timestamp_trunc(ingested_at, day) > (select timestamp(pdate) from latest_loaded_partition)
        {% endif %}
    ),
    mave_base as (
        select
            *,
            ingested_at as max_ingested_at
        from {{ source('mave', 'metadata') }}
        {% if is_incremental() %}
            where timestamp_trunc(ingested_at, day) > (select timestamp(pdate) from latest_loaded_partition)
        {% endif %}
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
            {% if is_incremental() %}
                and timestamp_trunc(ingested_at, day) > (select timestamp(pdate) from latest_loaded_partition)
            {% endif %}
        )
    ),
    
    -- Identify the superset of columns across all sources
    -- Note: BigQuery's UNION ALL BY NAME requires all branches to have the SAME columns.
    -- We'll explicitly select and cast to ensure alignment.
    
    crispr as (
        select
            * except (ingested_at),
            cast(null as string) as perturbation_name,
            cast(null as string) as guide_sequence
        from crispr_base
    ),
    mave as (
        select
            * except (ingested_at),
            cast(null as string) as guide_sequence
        from mave_base
    ),
    ps as (
        select
            * except (ingested_at),
            cast(null as string) as sample_id,
            cast(null as string) as significant,
            cast(null as string) as significance_criteria,
            cast(number_of_perturbed_targets as string) as number_of_perturbed_targets,
            cast(number_of_perturbed_samples as string) as number_of_perturbed_samples,
            cast(library_total_grnas as string) as library_total_grnas
        from ps_base
    ),
    
    unified as (
        select * from crispr
        union all by name
        select * from mave
        union all by name
        select * from ps
    )
    
select * from unified
