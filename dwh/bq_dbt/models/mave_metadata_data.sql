{{
    config(
        materialized="incremental",
        incremental_strategy="insert_overwrite",
        unique_key=["dataset_id", "sample_id"],
        on_schema_change="sync_all_columns",
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
            and partition_id not in ('__NULL__', '__UNPARTITIONED__')  -- ignore null + unpartitioned pseudo-ids
    ),
    base as (
        select
            m.* except (ingested_at),
            d.score_name,
            d.score_value,
            -- Keep as DATETIME
            greatest(m.ingested_at, d.ingested_at) as max_ingested_at
        from {{ source("mave", "metadata") }} as m
        join
            {{ source("mave", "data") }} as d
            on m.dataset_id = d.dataset_id
            and m.sample_id = d.sample_id

        {% if is_incremental() %}
            where
                timestamp_trunc(m.ingested_at, day)
                > (select timestamp(pdate) from latest_loaded_partition)
                and timestamp_trunc(d.ingested_at, day)
                > (select timestamp(pdate) from latest_loaded_partition)
        {% endif %}
    )

select *
from base
