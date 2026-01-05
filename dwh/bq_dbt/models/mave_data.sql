{{
    config(
        materialized="incremental",
        incremental_strategy="insert_overwrite",
        unique_key=["dataset_id", "sample_id", "score_name"],
        partition_by={
            "field": "max_ingested_at",
            "data_type": "timestamp",
            "granularity": "day",
        },
        cluster_by=["dataset_id", "sample_id", "perturbed_target_symbol"],
    )
}}

with
    meta as (
        select dataset_id, sample_id, perturbed_target_symbol, ingested_at
        from {{ source('mave', 'metadata') }}
    ),
    data as (
        select dataset_id, sample_id, score_name, score_value, ingested_at
        from {{ source('mave', 'data') }}
    )
select
    d.dataset_id,
    d.sample_id,
    m.perturbed_target_symbol,
    d.score_name,
    d.score_value,
    greatest(d.ingested_at, m.ingested_at) as max_ingested_at
from data d
join meta m
    on d.dataset_id = m.dataset_id and d.sample_id = m.sample_id
{% if is_incremental() %}
    where timestamp_trunc(d.ingested_at, day) >= (
        select timestamp_sub(timestamp_trunc(max(max_ingested_at), day), interval 1 day)
        from {{ this }}
    )
{% endif %}
