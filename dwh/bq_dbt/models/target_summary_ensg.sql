{{ config(materialized="table") }}

with
    target_aliases as (
        select
            ot.ensembl_gene_id,
            lower(trim(alias)) as alias_normalized
        from {{ ref("stg_opentargets_targets") }} as ot, unnest(ot.exact_aliases) as alias
        where alias is not null and alias != ''
    ),

    unique_aliases as (
        select
            alias_normalized,
            any_value(ensembl_gene_id) as ensembl_gene_id
        from target_aliases
        group by alias_normalized
        having count(distinct ensembl_gene_id) = 1
    ),

    catalogue_symbols_raw as (
        select perturbed_target_symbol
        from {{ ref("unified_metadata") }}
        union distinct
        select perturbed_target_symbol
        from {{ source("mave", "data") }}
        union distinct
        select perturbed_target_symbol
        from {{ source("perturb_seq", "pertpy_dea") }}
        union distinct
        select perturbed_target_symbol
        from {{ source("crispr", "data") }}
        union distinct
        select perturbed_target_symbol
        from {{ source("perturb_seq", "pertpy_gsea") }}
    ),

    mapped_symbols as (
        select
            s.perturbed_target_symbol,
            a.ensembl_gene_id
        from catalogue_symbols_raw as s
        inner join unique_aliases as a
            on lower(trim(s.perturbed_target_symbol)) = a.alias_normalized
        where s.perturbed_target_symbol is not null
    ),

    agg_meta as (
        select
            m.ensembl_gene_id,
            array_agg(distinct metadata.data_modality ignore nulls) as data_modalities,
            array_agg(distinct metadata.tissue_label ignore nulls) as tissues_tested,
            array_agg(distinct metadata.cell_type_label ignore nulls) as cell_types_tested,
            array_agg(distinct metadata.cell_line_label ignore nulls) as cell_lines_tested,
            array_agg(distinct metadata.sex_label ignore nulls) as sex_tested,
            array_agg(
                distinct metadata.developmental_stage_label ignore nulls
            ) as developmental_stages_tested,
            array_agg(distinct metadata.disease_label ignore nulls) as diseases_tested,
            array_agg(distinct metadata.license_label ignore nulls) as license
        from {{ ref("unified_metadata") }} as metadata
        inner join mapped_symbols as m
            on metadata.perturbed_target_symbol = m.perturbed_target_symbol
        group by m.ensembl_gene_id
    ),

    agg_crispr as (
        select
            m.ensembl_gene_id,
            count(distinct crispr.dataset_id) as n_crispr,
            countif(crispr.significant = 'True') as n_sig_crispr
        from {{ source("crispr", "data") }} as crispr
        inner join mapped_symbols as m
            on crispr.perturbed_target_symbol = m.perturbed_target_symbol
        group by m.ensembl_gene_id
    ),

    agg_mave as (
        select
            m.ensembl_gene_id,
            count(distinct mave.dataset_id) as n_mave
        from {{ source("mave", "data") }} as mave
        inner join mapped_symbols as m
            on mave.perturbed_target_symbol = m.perturbed_target_symbol
        group by m.ensembl_gene_id
    ),

    agg_ps as (
        select
            m.ensembl_gene_id,
            count(distinct dea.dataset_id) as n_perturb_seq,
            countif(dea.padj <= 0.05 and dea.log2FoldChange > 0) as n_sig_perturb_pairs_up,
            countif(dea.padj <= 0.05 and dea.log2FoldChange < 0) as n_sig_perturb_pairs_down
        from {{ source("perturb_seq", "pertpy_dea") }} as dea
        inner join mapped_symbols as m
            on dea.perturbed_target_symbol = m.perturbed_target_symbol
        group by m.ensembl_gene_id
    ),

    gsea_ranked as (
        select
            m.ensembl_gene_id,
            gsea.term,
            gsea.sidak,
            row_number() over (
                partition by m.ensembl_gene_id order by gsea.sidak asc
            ) as rn
        from {{ source("perturb_seq", "pertpy_gsea") }} as gsea
        inner join mapped_symbols as m
            on gsea.perturbed_target_symbol = m.perturbed_target_symbol
        where gsea.sidak <= 0.05
    ),

    agg_gsea as (
        select
            ensembl_gene_id,
            array_agg(term order by sidak asc) as top_gsea_terms
        from gsea_ranked
        where rn <= 5
        group by ensembl_gene_id
    ),

    targets as (
        select ensembl_gene_id
        from agg_meta
        union distinct
        select ensembl_gene_id
        from agg_mave
        union distinct
        select ensembl_gene_id
        from agg_ps
        union distinct
        select ensembl_gene_id
        from agg_crispr
        union distinct
        select ensembl_gene_id
        from agg_gsea
    )

select
    s.ensembl_gene_id,
    ot.approved_symbol,
    ot.approved_name,
    ot.exact_aliases,
    ot.search_keywords,
    coalesce(p.n_perturb_seq, 0) as n_perturb_seq,
    coalesce(p.n_sig_perturb_pairs_up, 0) as n_sig_perturb_pairs_up,
    coalesce(p.n_sig_perturb_pairs_down, 0) as n_sig_perturb_pairs_down,
    coalesce(c.n_crispr, 0) as n_crispr,
    coalesce(c.n_sig_crispr, 0) as n_sig_crispr,
    coalesce(v.n_mave, 0) as n_mave,
    g.top_gsea_terms,
    m.license,
    m.data_modalities,
    m.tissues_tested,
    m.cell_types_tested,
    m.cell_lines_tested,
    m.sex_tested,
    m.developmental_stages_tested,
    m.diseases_tested
from targets as s
inner join {{ ref("stg_opentargets_targets") }} as ot
    on s.ensembl_gene_id = ot.ensembl_gene_id
left join agg_meta as m on s.ensembl_gene_id = m.ensembl_gene_id
left join agg_mave as v on s.ensembl_gene_id = v.ensembl_gene_id
left join agg_ps as p on s.ensembl_gene_id = p.ensembl_gene_id
left join agg_crispr as c on s.ensembl_gene_id = c.ensembl_gene_id
left join agg_gsea as g on s.ensembl_gene_id = g.ensembl_gene_id
