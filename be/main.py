from fastapi import FastAPI, Query, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional, List, Dict, Any
from elasticsearch import Elasticsearch
import os
from dotenv import load_dotenv
import re
from urllib.parse import urlparse

try:
    from .models import (  # type: ignore
        SearchRequest,
        SearchResponse,
        FacetValue,
        Facets,
        LandingPageSummary,
    )
except ImportError:  # pragma: no cover - fallback for running as a script
    from models import (  # type: ignore
        SearchRequest,
        SearchResponse,
        FacetValue,
        Facets,
        LandingPageSummary,
    )

# Import data query APIs.
from data_query import lifespan as data_query_lifespan, router as data_query_router

load_dotenv()

app = FastAPI(title="Search API", version="1.0.0", lifespan=data_query_lifespan)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(data_query_router)
ES_URL = os.getenv("ES_URL")
ES_USERNAME = os.getenv("ES_USERNAME")
ES_PASSWORD = os.getenv("ES_PASSWORD")
ES_INDEX = "target-summary-v2"


es = None
if ES_URL and ES_USERNAME and ES_PASSWORD:
    es_config = {"hosts": [ES_URL], "basic_auth": (ES_USERNAME, ES_PASSWORD)}
    es = Elasticsearch(**es_config)
else:
    print(
        "WARNING: Missing one or more Elasticsearch environment variables (ES_URL, ES_USERNAME, ES_PASSWORD). Elasticsearch client not created."
    )


# Facet fields
FACET_FIELDS = [
    "data_modalities",
    "tissues_tested",
    "cell_types_tested",
    "cell_lines_tested",
    "sex_tested",
    "developmental_stages_tested",
    "diseases_tested",
]


# Elasticsearch helper functions
def _escape_wildcard(value: str) -> str:
    """Escape characters that have special meaning in wildcard queries."""
    return re.sub(r"([\\*?])", r"\\\1", value)


def build_elasticsearch_query(
    query: Optional[str], filters: Optional[Dict[str, List[str]]]
) -> Dict[str, Any]:
    """Build Elasticsearch query with search and filters"""
    filter_clauses = []
    should_clauses = []

    # Text search on perturbed_target_symbol
    if query:
        cleaned_query = query.strip()
        if cleaned_query:
            # Exact/fuzzy matches with boost
            should_clauses.append(
                {
                    "multi_match": {
                        "query": cleaned_query,
                        "fields": [
                            "perturbed_target_symbol^2",  # Boost exact match
                            "perturbed_target_symbol.text^1.5",  # Boost text field
                        ],
                        "type": "best_fields",
                        "fuzziness": "AUTO",
                    }
                }
            )

            # Prefix support for token beginnings (e.g. "SU" -> "SUMO1")
            should_clauses.append(
                {
                    "match_phrase_prefix": {
                        "perturbed_target_symbol.text": {
                            "query": cleaned_query,
                            "slop": 1,
                            "boost": 1.2,
                        }
                    }
                }
            )

            # Wildcard for partial/infix search (case-insensitive)
            wildcard_terms = []
            for term in cleaned_query.split():
                safe_term = _escape_wildcard(term.lower())
                if safe_term:
                    wildcard_terms.append(f"*{safe_term}*")

            # Include whole query if no spaces
            if not wildcard_terms:
                safe_term = _escape_wildcard(cleaned_query.lower())
                if safe_term:
                    wildcard_terms.append(f"*{safe_term}*")

            for wildcard_value in wildcard_terms:
                should_clauses.append(
                    {
                        "wildcard": {
                            "perturbed_target_symbol.keyword": {
                                "value": wildcard_value,
                                "case_insensitive": True,
                                "boost": 0.8,
                            }
                        }
                    }
                )
                should_clauses.append(
                    {
                        "wildcard": {
                            "perturbed_target_symbol": {
                                "value": wildcard_value,
                                "case_insensitive": True,
                                "boost": 0.6,
                            }
                        }
                    }
                )

    # Filters for facet fields
    if filters:
        for field, values in filters.items():
            if field in FACET_FIELDS and values:
                filter_clauses.append({"terms": {field: values}})

    bool_query: Dict[str, Any] = {}

    if filter_clauses:
        bool_query["filter"] = filter_clauses

    if should_clauses:
        bool_query["should"] = should_clauses
        bool_query["minimum_should_match"] = 1

    if bool_query:
        return {"bool": bool_query}

    return {"match_all": {}}


def build_aggregations() -> Dict[str, Any]:
    """Build aggregations for all facet fields"""
    aggs = {}
    for field in FACET_FIELDS:
        aggs[field] = {
            "terms": {
                "field": field,
                "size": 100,  # Adjust if you need more facet values
            }
        }
    return aggs


def parse_filters_from_params(
    data_modalities: Optional[str] = None,
    tissues_tested: Optional[str] = None,
    cell_types_tested: Optional[str] = None,
    cell_lines_tested: Optional[str] = None,
    sex_tested: Optional[str] = None,
    developmental_stages_tested: Optional[str] = None,
    diseases_tested: Optional[str] = None,
) -> Optional[Dict[str, List[str]]]:
    """Parse filters from query parameters"""
    filters = {}
    if data_modalities:
        filters["data_modalities"] = [v.strip() for v in data_modalities.split(",")]
    if tissues_tested:
        filters["tissues_tested"] = [v.strip() for v in tissues_tested.split(",")]
    if cell_types_tested:
        filters["cell_types_tested"] = [v.strip() for v in cell_types_tested.split(",")]
    if cell_lines_tested:
        filters["cell_lines_tested"] = [v.strip() for v in cell_lines_tested.split(",")]
    if sex_tested:
        filters["sex_tested"] = [v.strip() for v in sex_tested.split(",")]
    if developmental_stages_tested:
        filters["developmental_stages_tested"] = [
            v.strip() for v in developmental_stages_tested.split(",")
        ]
    if diseases_tested:
        filters["diseases_tested"] = [v.strip() for v in diseases_tested.split(",")]

    return filters if filters else None


async def perform_search(
    query: Optional[str], filters: Optional[Dict[str, List[str]]], page: int, size: int
) -> SearchResponse:
    """Perform the actual search operation"""
    # Build Elasticsearch query
    es_query = build_elasticsearch_query(query, filters)
    aggs = build_aggregations()

    # Calculate from/size for pagination
    from_ = (page - 1) * size

    # Execute search
    try:
        response = es.search(
            index=ES_INDEX, query=es_query, aggs=aggs, from_=from_, size=size
        )
    except Exception as e:
        error_detail = str(e)
        safe_config = {
            k: v for k, v in es_config.items() if k not in ["basic_auth", "api_key"]
        }
        safe_config["has_auth"] = "basic_auth" in es_config or "api_key" in es_config
        print(f"Elasticsearch config: {safe_config}")
        print(f"Elasticsearch error: {error_detail}")
        raise HTTPException(
            status_code=500, detail=f"Elasticsearch error: {error_detail}"
        )

    # Process results
    hits = response.get("hits", {})
    total = hits.get("total", {}).get("value", 0)
    results = [hit["_source"] for hit in hits.get("hits", [])]

    # Process facets
    facets_dict = {}
    aggregations = response.get("aggregations", {})
    for field in FACET_FIELDS:
        buckets = aggregations.get(field, {}).get("buckets", [])
        facets_dict[field] = [
            FacetValue(value=bucket["key"], count=bucket["doc_count"])
            for bucket in buckets
        ]

    facets = Facets(**facets_dict)

    # Calculate total pages
    total_pages = (total + size - 1) // size if total > 0 else 0

    return SearchResponse(
        total=total,
        page=page,
        size=size,
        total_pages=total_pages,
        results=results,
        facets=facets,
    )


@app.get("/")
async def root():
    return {"message": "Search API", "version": "1.0.0"}


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    es_status = "unknown"
    es_error = None

    # Check Elasticsearch
    try:
        if es.ping():
            es_status = "connected"
        else:
            es_status = "not_connected"
    except Exception as e:
        es_status = "not_connected"
        es_error = str(e)

    overall_status = "healthy" if es_status == "connected" else "unhealthy"

    return {
        "status": overall_status,
        "elasticsearch": {
            "status": es_status,
            "host": urlparse(ES_URL).hostname,
            "index": ES_INDEX,
            "error": es_error,
        },
    }


@app.get("/summary", response_model=LandingPageSummary)
async def get_landing_page_summary():
    """
    Retrieve the landing page summary document from Elasticsearch.
    """
    try:
        response = es.get(index="landing-page-summary", id="summary")
    except Exception as exc:
        raise HTTPException(
            status_code=500, detail=f"Elasticsearch error: {exc}"
        ) from exc

    source = response.get("_source")
    if not source:
        raise HTTPException(status_code=404, detail="Summary document not found")

    return LandingPageSummary(**source)


@app.get("/search", response_model=SearchResponse)
async def search_get(
    query: Optional[str] = Query(
        None, description="Search query for perturbed_target_symbol"
    ),
    data_modalities: Optional[str] = Query(
        None, description="Comma-separated list of data_modalities"
    ),
    tissues_tested: Optional[str] = Query(
        None, description="Comma-separated list of tissues_tested"
    ),
    cell_types_tested: Optional[str] = Query(
        None, description="Comma-separated list of cell_types_tested"
    ),
    cell_lines_tested: Optional[str] = Query(
        None, description="Comma-separated list of cell_lines_tested"
    ),
    sex_tested: Optional[str] = Query(
        None, description="Comma-separated list of sex_tested"
    ),
    developmental_stages_tested: Optional[str] = Query(
        None, description="Comma-separated list of developmental_stages_tested"
    ),
    diseases_tested: Optional[str] = Query(
        None, description="Comma-separated list of diseases_tested"
    ),
    page: int = Query(1, ge=1, description="Page number (1-indexed)"),
    size: int = Query(6, ge=1, le=100, description="Number of results per page"),
):
    """
    Search endpoint (GET) with pagination, facets, filtering, and text search.
    """
    filters = parse_filters_from_params(
        data_modalities,
        tissues_tested,
        cell_types_tested,
        cell_lines_tested,
        sex_tested,
        developmental_stages_tested,
        diseases_tested,
    )
    return await perform_search(query, filters, page, size)


@app.post("/search", response_model=SearchResponse)
async def search_post(request: SearchRequest):
    """
    Search endpoint (POST) with pagination, facets, filtering, and text search.
    Accepts JSON body with search parameters.
    """
    return await perform_search(
        request.query, request.filters, request.page, request.size
    )


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
