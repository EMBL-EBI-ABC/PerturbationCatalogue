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

with
    base as (
        select
            d.* except (dataset_id, perturbation, ingested_at),
            m.* except (ingested_at),
            greatest(
                datetime(m.ingested_at), datetime(d.ingested_at)
            ) as max_ingested_at
        from {{ source("perturb_seq", "data") }} as d
        right join
            (
                select distinct * except (sample_id)
                from {{ source("perturb_seq", "metadata") }}
                where
                    perturbed_target_symbol not like 'control%'
                    and perturbed_target_symbol not like '%None%'
            ) as m
            on m.dataset_id = d.dataset_id
            and m.perturbed_target_symbol = d.perturbation

        {% if is_incremental() %}
            -- No late arrivals: only load rows newer than what we've already loaded.
            -- Using a strict greater-than avoids reprocessing the last batch.
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
