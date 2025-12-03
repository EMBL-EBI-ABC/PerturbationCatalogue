{{
    config(
        materialized="incremental",
        incremental_strategy="merge",
        unique_key=["dataset_id"],
        on_schema_change="sync_all_columns",
        partition_by={
            "field": "max_ingested_at",
            "data_type": "timestamp",
            "granularity": "day",
        },
        cluster_by=["dataset_id"],
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
        select
            dataset_id,
            array_agg(distinct data_modality ignore nulls) as data_modalities,
            array_agg(distinct timepoint ignore nulls) as timepoints,
            array_agg(distinct treatment_label ignore nulls) as treatment_labels,
            array_agg(distinct treatment_id ignore nulls) as treatment_ids,
            array_agg(distinct model_system_label ignore nulls) as model_system_labels,
            array_agg(distinct model_system_id ignore nulls) as model_system_ids,
            array_agg(distinct tissue_label ignore nulls) as tissue_labels,
            array_agg(distinct tissue_id ignore nulls) as tissue_ids,
            array_agg(distinct cell_type_label ignore nulls) as cell_type_labels,
            array_agg(distinct cell_type_id ignore nulls) as cell_type_ids,
            array_agg(distinct cell_line_label ignore nulls) as cell_line_labels,
            array_agg(distinct cell_line_id ignore nulls) as cell_line_ids,
            array_agg(distinct sex_label ignore nulls) as sex_labels,
            array_agg(distinct sex_id ignore nulls) as sex_ids,
            array_agg(
                distinct developmental_stage_label ignore nulls
            ) as developmental_stage_labels,
            array_agg(
                distinct developmental_stage_id ignore nulls
            ) as developmental_stage_ids,
            array_agg(distinct disease_label ignore nulls) as disease_labels,
            array_agg(distinct disease_id ignore nulls) as disease_ids,
            any_value(study_title) as study_title,
            any_value(study_uri) as study_uri,
            any_value(study_year) as study_year,
            any_value(first_author) as first_author,
            any_value(last_author) as last_author,
            any_value(experiment_title) as experiment_title,
            any_value(experiment_summary) as experiment_summary,
            array_agg(
                distinct number_of_perturbed_targets ignore nulls
            ) as number_of_perturbed_targets,
            array_agg(
                distinct number_of_perturbed_samples ignore nulls
            ) as number_of_perturbed_samples,
            array_agg(
                distinct library_generation_type_label ignore nulls
            ) as library_generation_type_labels,
            array_agg(
                distinct library_generation_type_id ignore nulls
            ) as library_generation_type_ids,
            array_agg(
                distinct library_generation_method_label ignore nulls
            ) as library_generation_method_labels,
            array_agg(
                distinct library_generation_method_id ignore nulls
            ) as library_generation_method_ids,
            array_agg(
                distinct enzyme_delivery_method_label ignore nulls
            ) as enzyme_generation_delivery_labels,
            array_agg(
                distinct enzyme_delivery_method_id ignore nulls
            ) as enzyme_generation_delivery_ids,
            array_agg(
                distinct library_delivery_method_label ignore nulls
            ) as library_generation_delivery_labels,
            array_agg(
                distinct library_delivery_method_id ignore nulls
            ) as library_generation_delivery_ids,
            array_agg(
                distinct enzyme_integration_state_label ignore nulls
            ) as enzyme_integration_state_labels,
            array_agg(
                distinct enzyme_integration_state_id ignore nulls
            ) as enzyme_integration_state_ids,
            array_agg(
                distinct library_integration_state_label ignore nulls
            ) as library_integration_state_labels,
            array_agg(
                distinct library_integration_state_id ignore nulls
            ) as library_integration_state_ids,
            array_agg(
                distinct enzyme_expression_control_label ignore nulls
            ) as enzyme_expression_control_labels,
            array_agg(
                distinct enzyme_expression_control_id ignore nulls
            ) as enzyme_expression_control_ids,
            array_agg(
                distinct library_expression_control_label ignore nulls
            ) as library_expression_control_labels,
            array_agg(
                distinct library_expression_control_id ignore nulls
            ) as library_expression_control_ids,
            array_agg(distinct library_name ignore nulls) as library_names,
            array_agg(distinct library_uri ignore nulls) as library_uris,
            array_agg(
                distinct library_format_label ignore nulls
            ) as library_format_labels,
            array_agg(distinct library_format_id ignore nulls) as library_format_ids,
            array_agg(
                distinct library_scope_label ignore nulls
            ) as library_scope_labels,
            array_agg(distinct library_scope_id ignore nulls) as library_scope_ids,
            array_agg(
                distinct library_perturbation_type_label ignore nulls
            ) as library_perturbation_type_labels,
            array_agg(
                distinct library_perturbation_type_id ignore nulls
            ) as library_perturbation_type_ids,
            array_agg(
                distinct library_manufacturer ignore nulls
            ) as library_manufacturers,
            array_agg(
                distinct library_lentiviral_generation ignore nulls
            ) as library_lentivaral_generations,
            array_agg(
                distinct library_grnas_per_target ignore nulls
            ) as library_grnas_per_targets,
            array_agg(distinct library_total_grnas ignore nulls) as library_total_grnas,
            array_agg(
                distinct library_total_variants ignore nulls
            ) as library_total_variants,
            array_agg(
                distinct readout_dimensionality_label ignore nulls
            ) as readout_dimensionality_labels,
            array_agg(
                distinct readout_dimensionality_id ignore nulls
            ) as readout_dimensionality_ids,
            array_agg(distinct readout_type_label ignore nulls) as readout_type_labels,
            array_agg(distinct readout_type_id ignore nulls) as readout_type_ids,
            array_agg(
                distinct readout_technology_label ignore nulls
            ) as readout_technology_labels,
            array_agg(
                distinct readout_technology_id ignore nulls
            ) as readout_technology_ids,
            array_agg(distinct method_name_label ignore nulls) as method_name_labels,
            array_agg(distinct method_name_id ignore nulls) as method_name_ids,
            array_agg(distinct method_uri ignore nulls) as method_uri,
            array_agg(
                distinct sequencing_library_kit_label ignore nulls
            ) as sequencing_library_kit_labels,
            array_agg(
                distinct sequencing_library_kit_id ignore nulls
            ) as sequencing_library_ki_ids,
            array_agg(
                distinct sequencing_platform_label ignore nulls
            ) as sequencing_platform_labels,
            array_agg(
                distinct sequencing_platform_id ignore nulls
            ) as sequencing_platform_ids,
            array_agg(
                distinct sequencing_strategy_label ignore nulls
            ) as sequencing_strategy_labels,
            array_agg(
                distinct sequencing_strategy_id ignore nulls
            ) as sequencing_strategy_ids,
            array_agg(
                distinct software_counts_label ignore nulls
            ) as software_counts_labels,
            array_agg(distinct software_counts_id ignore nulls) as software_counts_ids,
            array_agg(
                distinct software_analysis_label ignore nulls
            ) as software_analysis_labels,
            array_agg(
                distinct software_analysis_id ignore nulls
            ) as software_analysis_ids,
            array_agg(
                distinct reference_genome_label ignore nulls
            ) as reference_genome_labels,
            array_agg(
                distinct reference_genome_id ignore nulls
            ) as reference_genome_ids,
            array_agg(distinct license_label ignore nulls) as license_labels,
            array_agg(distinct license_id ignore nulls) as license_ids,
            array_agg(distinct associated_datasets ignore nulls) as associated_datasets,
            score_interpretation,
            max(max_ingested_at) as max_ingested_at
        from {{ ref("unified_metadata_data") }}

        {% if is_incremental() %}
            -- No late arrivals: only load rows newer than what we've already loaded.
            -- Using a strict greater-than avoids reprocessing the last batch.
            where
                timestamp_trunc(max_ingested_at, day)
                > (select timestamp(pdate) from latest_loaded_partition)
        {% endif %}
        group by dataset_id
    )

select *
from base
