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
    meta as (
        select dataset_id, sample_id, perturbed_target_symbol, significant, significance_criteria, ingested_at
        from {{ source('crispr', 'metadata') }}
    ),
    data as (
        select dataset_id, sample_id, score_name, score_value, ingested_at
        from {{ source('crispr', 'data') }}
    )
select
    d.dataset_id,
    d.sample_id,
    m.perturbed_target_symbol,
    d.score_name,
    d.score_value,
    m.significant,
    m.significance_criteria,
    greatest(d.ingested_at, m.ingested_at) as max_ingested_at
from data d
join meta m
    on d.dataset_id = m.dataset_id and d.sample_id = m.sample_id

