{{
    config(
        materialized="incremental",
        incremental_strategy="merge",
        unique_key=["dataset_id", "sample_id", "score_name"],
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
    sample_id,
    perturbed_target_symbol,
    score_name,
    score_value,
    significant,
    significance_criteria,
    max_ingested_at
from {{ ref("unified_metadata_data") }}
where
    data_modality = "CRISPR screen"
    {% if is_incremental() %}
        and datetime(max_ingested_at) > (
            select coalesce(max(max_ingested_at), datetime '1970-01-01 00:00:00')
            from {{ this }}
        )
    {% endif %}
