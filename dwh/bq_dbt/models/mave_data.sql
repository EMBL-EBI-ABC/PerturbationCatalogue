{{ config(materialized="table") }}

select
    dataset_id,
    sample_id,
    perturbed_target_ensg,
    score_name,
    score_value,
    perturbation_name,
    cast(regexp_extract(perturbation_name, r'p\.[a-zA-Z]+(\d+)') as int64) as perturbation_position,
    regexp_extract(perturbation_name, r'p\.([a-zA-Z]+)\d+') as perturbation_aa_wt,
    regexp_extract(perturbation_name, r'p\.[a-zA-Z]+\d+([a-zA-Z=]+)') as perturbation_aa_change,
    ingested_at as max_ingested_at
from {{ source("mave", "data") }}
