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
    )
select
    dataset_id,
    sample_id,
    perturbation_name,
    perturbed_target_symbol,
    score_name,
    score_value,
    max_ingested_at
from {{ ref("mave_metadata_data") }}
{% if is_incremental() %}
    where
        timestamp_trunc(max_ingested_at, day)
        > (select timestamp(pdate) from latest_loaded_partition)
{% endif %}
