{{ config(materialized="table") }}

with
    targets as (
        select
            id as ensembl_gene_id,
            approvedSymbol as approved_symbol,
            approvedName as approved_name,
            biotype,
            genomicLocation as genomic_location,
            symbolSynonyms.list as symbol_synonyms,
            nameSynonyms.list as name_synonyms,
            obsoleteSymbols.list as obsolete_symbols,
            obsoleteNames.list as obsolete_names,
            synonyms.list as synonyms
        from {{ source("reference", "opentargets_targets") }}
    ),

    exact_alias_rows as (
        select ensembl_gene_id, approved_symbol as alias
        from targets
        where approved_symbol is not null and approved_symbol != ''

        union distinct

        select t.ensembl_gene_id, synonym.element.label as alias
        from targets t, unnest(t.symbol_synonyms) as synonym
        where synonym.element.label is not null and synonym.element.label != ''

        union distinct

        select t.ensembl_gene_id, obsolete_symbol.element.label as alias
        from targets t, unnest(t.obsolete_symbols) as obsolete_symbol
        where
            obsolete_symbol.element.label is not null
            and obsolete_symbol.element.label != ''
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

        select t.ensembl_gene_id, synonym.element.label as keyword
        from targets t, unnest(t.symbol_synonyms) as synonym
        where synonym.element.label is not null and synonym.element.label != ''

        union distinct

        select t.ensembl_gene_id, obsolete_symbol.element.label as keyword
        from targets t, unnest(t.obsolete_symbols) as obsolete_symbol
        where
            obsolete_symbol.element.label is not null
            and obsolete_symbol.element.label != ''

        union distinct

        select t.ensembl_gene_id, name_synonym.element.label as keyword
        from targets t, unnest(t.name_synonyms) as name_synonym
        where name_synonym.element.label is not null and name_synonym.element.label != ''

        union distinct

        select t.ensembl_gene_id, obsolete_name.element.label as keyword
        from targets t, unnest(t.obsolete_names) as obsolete_name
        where
            obsolete_name.element.label is not null
            and obsolete_name.element.label != ''

        union distinct

        select t.ensembl_gene_id, synonym.element.label as keyword
        from targets t, unnest(t.synonyms) as synonym
        where synonym.element.label is not null and synonym.element.label != ''
    ),

    exact_aliases as (
        select ensembl_gene_id, array_agg(alias order by lower(alias)) as exact_aliases
        from (select distinct ensembl_gene_id, alias from exact_alias_rows)
        group by ensembl_gene_id
    ),

    search_keywords as (
        select
            ensembl_gene_id,
            array_agg(keyword order by lower(keyword)) as search_keywords
        from (select distinct ensembl_gene_id, keyword from search_keyword_rows)
        group by ensembl_gene_id
    )

select
    t.ensembl_gene_id,
    t.approved_symbol,
    t.approved_name,
    t.biotype,
    t.genomic_location,
    coalesce(e.exact_aliases, []) as exact_aliases,
    coalesce(k.search_keywords, []) as search_keywords
from targets t
left join exact_aliases e using (ensembl_gene_id)
left join search_keywords k using (ensembl_gene_id)
