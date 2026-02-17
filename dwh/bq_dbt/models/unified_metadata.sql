{{ config(materialized="table") }}

with
    unified as (
        select *
        from {{ source("crispr", "metadata") }} full outer
        union all by name
        select *
        from {{ source("mave", "metadata") }} full outer
        union all by name
        select distinct * except (sample_id)
        from {{ source("perturb_seq", "metadata") }}
        where not (
            -- Exclude rows where EVERY '|' separated part contains "control"
            not exists (
                select 1
                from unnest(split(perturbed_target_symbol, '|')) as part
                where lower(part) not like '%control%'
            )
        )
    )

select *
from unified
{% if var('suppress_datasets', None) %}
    where dataset_id not in (
        {% set suppressed_ids = var('suppress_datasets').split(',') %}
        {% for id in suppressed_ids %}
            '{{ id | trim }}'{% if not loop.last %}, {% endif %}
        {% endfor %}
    )
{% endif %}
