{{
    config(
        materialized="incremental",
        incremental_strategy="merge",
        unique_key=["dataset_id", "sample_id", "perturbed_target_symbol", "gene"],
        on_schema_change="sync_all_columns",
        partition_by={
            "field": "max_ingested_at",
            "data_type": "datetime",
            "granularity": "day",
        },
    )
}}

with
    base as (

        {% if is_incremental() %}

            -- Only pull rows newer than what's already loaded into
            -- unified_metadata_data
            select *
            from {{ ref("crispr_metadata_data") }}
            where
                datetime(max_ingested_at) > (
                    select
                        coalesce(max(max_ingested_at), datetime '1970-01-01 00:00:00')
                    from {{ this }}
                )
                full outer
            union all by name
            select
                * except (significant, significance_criteria),
                null as significant,
                null as significance_criteria
            from {{ ref("perturb_seq_metadata_data") }}
            where
                datetime(max_ingested_at) > (
                    select
                        coalesce(max(max_ingested_at), datetime '1970-01-01 00:00:00')
                    from {{ this }}
                )

        {% else %}

            -- First run: take everything
            select *
            from {{ ref("crispr_metadata_data") }} full outer
            union all by name
            select
                * except (significant, significance_criteria),
                null as significant,
                null as significance_criteria
            from {{ ref("perturb_seq_metadata_data") }}

        {% endif %}

    )

select *
from base
