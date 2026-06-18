{{ config(materialized="table") }}

select
    dataset_id,
    sample_id,
    perturbed_target_ensg,
    score_name,
    score_value,
    significant,
    significance_criteria,
    ingested_at as max_ingested_at
from {{ source("crispr", "data") }}
