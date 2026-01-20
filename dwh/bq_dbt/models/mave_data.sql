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
    perturbation_name,
    cast(
        regexp_extract(perturbation_name, r'p\.[a-zA-Z]+(\d+)') as int64
    ) as perturbation_position,
    regexp_extract(perturbation_name, r'p\.([a-zA-Z]+)\d+') as perturbation_aa_wt,
    regexp_extract(
        perturbation_name, r'p\.[a-zA-Z]+\d+([a-zA-Z=]+)'
    ) as perturbation_aa_change,
    ingested_at as max_ingested_at
from {{ source("mave", "data") }}
