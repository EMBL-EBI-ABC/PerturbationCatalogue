{{ config(materialized="table") }}

with
    targets as (
        select
            id as ensembl_gene_id,
            approvedSymbol as approved_symbol,
            approvedName as approved_name,
            biotype,
            genomicLocation as genomic_location,
            symbolSynonyms as symbol_synonyms,
            nameSynonyms as name_synonyms,
            obsoleteSymbols as obsolete_symbols,
            obsoleteNames as obsolete_names,
            synonyms
        from {{ source("reference", "opentargets_targets") }}
    ),

    exact_alias_rows as (
        select ensembl_gene_id, approved_symbol as alias
        from targets
        where approved_symbol is not null and approved_symbol != ''

        union distinct

        select t.ensembl_gene_id, synonym.label as alias
        from targets t, unnest(t.symbol_synonyms) as synonym
        where synonym.label is not null and synonym.label != ''

        union distinct

        select t.ensembl_gene_id, obsolete_symbol.label as alias
        from targets t, unnest(t.obsolete_symbols) as obsolete_symbol
        where obsolete_symbol.label is not null and obsolete_symbol.label != ''
    ),

    search_keyword_rows as (
        select ensembl_gene_id, ensembl_gene_id as keyword
        from targets
        where ensembl_gene_id is not null and ensembl_gene_id != ''

        union distinct

        select ensembl_gene_id, approved_symbol as keyword
        from targets
        where approved_symbol is not null and approved_symbol != ''

        union distinct

        select ensembl_gene_id, approved_name as keyword
        from targets
        where approved_name is not null and approved_name != ''

        union distinct

        select t.ensembl_gene_id, synonym.label as keyword
        from targets t, unnest(t.symbol_synonyms) as synonym
        where synonym.label is not null and synonym.label != ''

        union distinct

        select t.ensembl_gene_id, obsolete_symbol.label as keyword
        from targets t, unnest(t.obsolete_symbols) as obsolete_symbol
        where obsolete_symbol.label is not null and obsolete_symbol.label != ''

        union distinct

        select t.ensembl_gene_id, name_synonym.label as keyword
        from targets t, unnest(t.name_synonyms) as name_synonym
        where name_synonym.label is not null and name_synonym.label != ''

        union distinct

        select t.ensembl_gene_id, obsolete_name.label as keyword
        from targets t, unnest(t.obsolete_names) as obsolete_name
        where obsolete_name.label is not null and obsolete_name.label != ''

        union distinct

        select t.ensembl_gene_id, synonym.label as keyword
        from targets t, unnest(t.synonyms) as synonym
        where synonym.label is not null and synonym.label != ''
    )

select
    t.ensembl_gene_id,
    t.approved_symbol,
    t.approved_name,
    t.biotype,
    t.genomic_location,
    array(
        select alias
        from exact_alias_rows e
        where e.ensembl_gene_id = t.ensembl_gene_id
        order by lower(alias)
    ) as exact_aliases,
    array(
        select keyword
        from search_keyword_rows k
        where k.ensembl_gene_id = t.ensembl_gene_id
        order by lower(keyword)
    ) as search_keywords
from targets t
