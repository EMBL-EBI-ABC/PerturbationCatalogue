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
        where
            perturbed_target_symbol not like 'control%'
            and perturbed_target_symbol not like '%None%'
    )

select *
from unified
