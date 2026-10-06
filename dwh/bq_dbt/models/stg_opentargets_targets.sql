{{ config(materialized="table") }}

with flattened as (
    select
        id as ensembl_gene_id,
        approvedSymbol as approved_symbol,
        approvedName as approved_name,
        biotype,
        genomicLocation as genomic_location,
        array(select element.label from unnest(symbolSynonyms.list)) as symbol_synonyms,
        array(select element.label from unnest(obsoleteSymbols.list)) as obsolete_symbols,
        array(select element.label from unnest(nameSynonyms.list)) as name_synonyms,
        array(select element.label from unnest(obsoleteNames.list)) as obsolete_names,
        array(select element.label from unnest(synonyms.list)) as synonyms
    from {{ source("reference", "opentargets_targets") }}
)

select
    ensembl_gene_id,
    approved_symbol,
    approved_name,
    biotype,
    genomic_location,
    array(
        select distinct label
        from unnest(array_concat([approved_symbol], symbol_synonyms, obsolete_symbols)) as label
        where label is not null and label != ''
        order by lower(label)
    ) as exact_aliases,
    array(
        select distinct label
        from unnest(array_concat([ensembl_gene_id, approved_symbol, approved_name], symbol_synonyms, obsolete_symbols, name_synonyms, obsolete_names, synonyms)) as label
        where label is not null and label != ''
        order by lower(label)
    ) as search_keywords
from flattened
