{{
    config(
        materialized="incremental",
        incremental_strategy="merge",
        unique_key=["dataset_id", "perturbed_target_symbol", "gene"],
        on_schema_change="sync_all_columns",
        partition_by={
            "field": "max_ingested_at",
            "data_type": "datetime",
            "granularity": "day",
        },
    )
}}

select
    dataset_id,
    perturbed_target_symbol,
    gene,
    log2foldchange,
    padj,
    basemean,
    max_ingested_at
from {{ ref("unified_metadata_data") }}
where
    data_modality = "Perturb-seq"
    {% if is_incremental() %}
        and datetime(max_ingested_at) > (
            select coalesce(max(max_ingested_at), datetime '1970-01-01 00:00:00')
            from {{ this }}
        )
    {% endif %}
