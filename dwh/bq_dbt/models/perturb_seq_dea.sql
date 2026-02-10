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
    perturbed_target_symbol,
    gene,
    padj,
    log2FoldChange as log2foldchange,
    score_name,
    score_value,
    cell_type,
    ingested_at as max_ingested_at
from {{ source("perturb_seq", "pertpy_dea") }}
