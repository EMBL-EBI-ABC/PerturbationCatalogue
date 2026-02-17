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
            -- Exclude rows where every '|' separated part contains 'control',
            -- meaning that the sample *only* contains controls.
            array_size(
                array_filter(
                    split(perturbed_target_symbol, '|'),
                    x -> lower(x) not like '%control%'
                )
            ) = 0
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
