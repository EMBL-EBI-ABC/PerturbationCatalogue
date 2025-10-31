{{
    config(
        materialized="incremental",
        incremental_strategy="insert_overwrite",
        unique_key=["dataset_id", "sample_id", "perturbed_target_symbol", "gene"],
        on_schema_change="sync_all_columns",
        partition_by={
            "field": "max_ingested_at",
            "data_type": "timestamp",
            "granularity": "day",
        },
        cluster_by=["dataset_id", "sample_id", "perturbed_target_symbol"],
    )
}}

with
    latest_loaded_partition as (
        select parse_date('%Y%m%d', max(partition_id)) as pdate
        from `{{ this.database }}`.`{{ this.schema }}.INFORMATION_SCHEMA.PARTITIONS`
        where
            table_name = '{{ this.identifier }}'
            and partition_id not in ('__NULL__', '__UNPARTITIONED__')  -- ignore null + unpartitioned pseudo-ids
    ),
    base as (
        {% if is_incremental() %}
            -- Only pull rows newer than what's already loaded into
            -- unified_metadata_data
            select *
            from {{ ref("crispr_metadata_data") }}
            where
                timestamp_trunc(max_ingested_at, day)
                > (select timestamp(pdate) from latest_loaded_partition) full outer
            union all by name
            select
                * except (
                    significant,
                    significance_criteria,
                    number_of_perturbed_targets,
                    number_of_perturbed_samples,
                    library_total_grnas
                ),
                null as significant,
                null as significance_criteria,
                cast(
                    number_of_perturbed_targets as string
                ) as number_of_perturbed_targets,
                cast(
                    number_of_perturbed_samples as string
                ) as number_of_perturbed_samples,
                cast(library_total_grnas as string) as library_total_grnas

            from {{ ref("perturb_seq_metadata_data") }}
            where
                timestamp_trunc(max_ingested_at, day)
                > (select timestamp(pdate) from latest_loaded_partition) full outer
            union all by name
            select *
            from {{ ref("mave_metadata_data") }}
            where
                timestamp_trunc(max_ingested_at, day)
                > (select timestamp(pdate) from latest_loaded_partition)

        {% else %}

            -- First run: take everything
            select *
            from {{ ref("crispr_metadata_data") }} full outer
            union all by name
            select
                * except (
                    significant,
                    significance_criteria,
                    number_of_perturbed_targets,
                    number_of_perturbed_samples,
                    library_total_grnas
                ),
                null as significant,
                null as significance_criteria,
                cast(
                    number_of_perturbed_targets as string
                ) as number_of_perturbed_targets,
                cast(
                    number_of_perturbed_samples as string
                ) as number_of_perturbed_samples,
                cast(library_total_grnas as string) as library_total_grnas
            from {{ ref("perturb_seq_metadata_data") }} full outer
            union all by name
            select *
            from {{ ref("mave_metadata_data") }}

        {% endif %}

    )

select *
from base
