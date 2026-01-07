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

select
    dataset_id,
    sample_id,
    perturbed_target_symbol,
    score_name,
    score_value,
    significant,
    significance_criteria,
    ingested_at as max_ingested_at
from {{ source('crispr', 'data') }}

