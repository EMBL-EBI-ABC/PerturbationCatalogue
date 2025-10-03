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

with
    base as (
        select
            m.* except (ingested_at),
            d.score_name,
            d.score_value,
            -- Keep as DATETIME
            greatest(
                datetime(m.ingested_at), datetime(d.ingested_at)
            ) as max_ingested_at
        from {{ source("crispr", "metadata") }} as m
        join
            {{ source("crispr", "data") }} as d
            on m.dataset_id = d.dataset_id
            and m.sample_id = d.sample_id

        {% if is_incremental() %}
            where
                datetime(m.ingested_at) > (
                    select
                        coalesce(max(max_ingested_at), datetime '1970-01-01 00:00:00')
                    from {{ this }}
                )
                and datetime(d.ingested_at) > (
                    select
                        coalesce(max(max_ingested_at), datetime '1970-01-01 00:00:00')
                    from {{ this }}
                )
        {% endif %}
    )

select *
from base
