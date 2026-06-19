import re
from typing import Any, Dict, List, Optional


TARGET_EXACT_SEARCH_FIELDS = {
    "ensembl_gene_id": 10.0,
    "approved_symbol": 8.0,
    "exact_aliases": 6.0,
}

TARGET_TEXT_SEARCH_FIELDS = [
    "ensembl_gene_id^10.0",
    "ensembl_gene_id.text^4.0",
    "approved_symbol^8.0",
    "approved_symbol.text^4.0",
    "exact_aliases^6.0",
    "exact_aliases.text^3.0",
    "approved_name^3.0",
    "approved_name.text^2.0",
    "search_keywords^2.0",
    "search_keywords.text^1.5",
]

TARGET_SEARCHABLE_FIELDS = [
    "ensembl_gene_id",
    "approved_symbol",
    "exact_aliases",
    "approved_name",
    "search_keywords",
]


def escape_wildcard(value: str) -> str:
    """Escape characters that have special meaning in wildcard queries."""
    return re.sub(r"([\\*?])", r"\\\1", value)


def _target_filter_clauses(
    filters: Optional[Dict[str, List[str]]], facet_fields: List[str]
) -> List[Dict[str, Any]]:
    if not filters:
        return []

    return [
        {"terms": {field: values}}
        for field, values in filters.items()
        if field in facet_fields and values
    ]


def _bool_or_match_all(
    should_clauses: List[Dict[str, Any]], filter_clauses: List[Dict[str, Any]]
) -> Dict[str, Any]:
    bool_query: Dict[str, Any] = {}

    if filter_clauses:
        bool_query["filter"] = filter_clauses

    if should_clauses:
        bool_query["should"] = should_clauses
        bool_query["minimum_should_match"] = 1

    if bool_query:
        return {"bool": bool_query}

    return {"match_all": {}}


def build_target_exact_query(
    query: str,
    filters: Optional[Dict[str, List[str]]] = None,
    facet_fields: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Build a target query that only matches exact ENSGs, symbols, and aliases."""
    cleaned_query = query.strip()
    filter_clauses = _target_filter_clauses(filters, facet_fields or [])
    should_clauses = []

    if cleaned_query:
        for field, boost in TARGET_EXACT_SEARCH_FIELDS.items():
            should_clauses.append(
                {
                    "term": {
                        field: {
                            "value": cleaned_query,
                            "case_insensitive": True,
                            "boost": boost,
                        }
                    }
                }
            )

    return _bool_or_match_all(should_clauses, filter_clauses)


def build_target_fuzzy_query(
    query: Optional[str],
    filters: Optional[Dict[str, List[str]]] = None,
    facet_fields: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Build the full target query with exact, fuzzy, prefix, and wildcard matches."""
    filter_clauses = _target_filter_clauses(filters, facet_fields or [])
    should_clauses = []

    if query:
        cleaned_query = query.strip()
        if cleaned_query:
            for field, boost in TARGET_EXACT_SEARCH_FIELDS.items():
                should_clauses.append(
                    {
                        "term": {
                            field: {
                                "value": cleaned_query,
                                "case_insensitive": True,
                                "boost": boost,
                            }
                        }
                    }
                )

            should_clauses.append(
                {
                    "multi_match": {
                        "query": cleaned_query,
                        "fields": TARGET_TEXT_SEARCH_FIELDS,
                        "type": "best_fields",
                        "fuzziness": "AUTO:5,8",
                    }
                }
            )

            for field in TARGET_SEARCHABLE_FIELDS:
                should_clauses.append(
                    {
                        "match_phrase_prefix": {
                            f"{field}.text": {
                                "query": cleaned_query,
                                "slop": 1,
                                "boost": 1.2,
                            }
                        }
                    }
                )

            wildcard_terms = []
            for term in cleaned_query.split():
                safe_term = escape_wildcard(term.lower())
                if safe_term:
                    wildcard_terms.append(f"*{safe_term}*")

            if not wildcard_terms:
                safe_term = escape_wildcard(cleaned_query.lower())
                if safe_term:
                    wildcard_terms.append(f"*{safe_term}*")

            for wildcard_value in wildcard_terms:
                for field in TARGET_SEARCHABLE_FIELDS:
                    should_clauses.append(
                        {
                            "wildcard": {
                                field: {
                                    "value": wildcard_value,
                                    "case_insensitive": True,
                                    "boost": 0.3,
                                }
                            }
                        }
                    )

    return _bool_or_match_all(should_clauses, filter_clauses)

