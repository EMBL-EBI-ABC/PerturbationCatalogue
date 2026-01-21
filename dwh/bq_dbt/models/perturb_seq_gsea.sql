{{
    config(
        materialized="table",
        partition_by={
            "field": "max_ingested_at",
            "data_type": "timestamp",
            "granularity": "day",
        },
        cluster_by=["dataset_id", "perturbed_target_symbol"],
    )
}}

select
    dataset_id,
    term,
    perturbed_target_symbol,
    es,
    nes,
    pval,
    sidak,
    fdr,
    geneset_size,
    leading_edge,
    cell_type,
    current_timestamp() as max_ingested_at
from {{ source("perturb_seq", "pertpy_gsea") }}
