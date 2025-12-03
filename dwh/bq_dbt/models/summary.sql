{{ config(materialized="table") }}


with
    base as (select * from {{ ref("unified_metadata_data") }}),
    -- Top-k helpers (by #datasets using that attribute)
    modalities as (
        select data_modality as value, count(distinct dataset_id) as n_datasets
        from base
        where data_modality is not null
        group by 1
        order by n_datasets desc
        limit 10
    ),
    tissues as (
        select tissue_label as value, count(distinct dataset_id) as n_datasets
        from base
        where tissue_label is not null
        group by 1
        order by n_datasets desc
        limit 10
    ),
    cell_types as (
        select cell_type_label as value, count(distinct dataset_id) as n_datasets
        from base
        where cell_type_label is not null
        group by 1
        order by n_datasets desc
        limit 10
    ),
    cell_lines as (
        select cell_line_label as value, count(distinct dataset_id) as n_datasets
        from base
        where cell_line_label is not null
        group by 1
        order by n_datasets desc
        limit 10
    ),
    sexes as (
        select sex_label as value, count(distinct dataset_id) as n_datasets
        from base
        where sex_label is not null
        group by 1
        order by n_datasets desc
        limit 10
    ),
    dev_stages as (
        select
            developmental_stage_label as value, count(distinct dataset_id) as n_datasets
        from base
        where developmental_stage_label is not null
        group by 1
        order by n_datasets desc
        limit 10
    ),
    diseases as (
        select disease_label as value, count(distinct dataset_id) as n_datasets
        from base
        where disease_label is not null
        group by 1
        order by n_datasets desc
        limit 10
    ),
    perturbation_types as (
        select
            perturbation_type_label as value, count(distinct dataset_id) as n_datasets
        from base
        where perturbation_type_label is not null
        group by 1
        order by n_datasets desc
        limit 10
    ),

    summary as (
        select
            count(distinct dataset_id) as n_datasets,
            count(distinct experiment_title) as n_experiments,
            min(study_year) as min_year,
            max(study_year) as max_year,
            count(distinct perturbed_target_ensg) as n_targets,
            count(distinct tissue_label) as n_tissues,
            count(distinct cell_type_label) as n_cell_types,
            count(distinct cell_line_label) as n_cell_lines,
            count(distinct disease_label) as n_diseases,
        from base
    )

select
    s.*,
    -- structured arrays for API use
    (select array_agg(struct(value, n_datasets)) from modalities) as top_modalities,
    (select array_agg(struct(value, n_datasets)) from tissues) as top_tissues,
    (select array_agg(struct(value, n_datasets)) from cell_types) as top_cell_types,
    (select array_agg(struct(value, n_datasets)) from cell_lines) as top_cell_lines,
    (
        select array_agg(struct(value, n_datasets)) from perturbation_types
    ) as top_perturbation_types,
    (select array_agg(struct(value, n_datasets)) from diseases) as top_diseases,
    (select array_agg(struct(value, n_datasets)) from sexes) as top_sexes,
    (select array_agg(struct(value, n_datasets)) from dev_stages) as top_dev_stages,
from summary s