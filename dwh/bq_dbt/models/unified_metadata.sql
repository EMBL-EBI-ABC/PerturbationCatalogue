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
    crispr as (
        select
            * except (ingested_at),
            ingested_at as max_ingested_at
        from {{ source('crispr', 'metadata') }}
        {% if is_incremental() %}
            where timestamp_trunc(ingested_at, day) > (select timestamp(pdate) from latest_loaded_partition)
        {% endif %}
    ),
    mave as (
        select
            * except (ingested_at),
            ingested_at as max_ingested_at
        from {{ source('mave', 'metadata') }}
        {% if is_incremental() %}
            where timestamp_trunc(ingested_at, day) > (select timestamp(pdate) from latest_loaded_partition)
        {% endif %}
    ),
    ps as (
        select
            * except (ingested_at),
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
    unified as (
        select * from crispr
        union all by name
        select * from mave
        union all by name
        select
            * except (
                significant,
                significance_criteria,
                number_of_perturbed_targets,
                number_of_perturbed_samples,
                library_total_grnas
            ),
            null as significant,
            null as significance_criteria,
            cast(number_of_perturbed_targets as string) as number_of_perturbed_targets,
            cast(number_of_perturbed_samples as string) as number_of_perturbed_samples,
            cast(library_total_grnas as string) as library_total_grnas
        from ps
    )
select * from unified
