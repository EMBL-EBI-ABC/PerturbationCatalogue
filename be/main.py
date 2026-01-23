from fastapi import FastAPI, Query, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional, List, Dict, Any
from elasticsearch import AsyncElasticsearch
import asyncpg
from dotenv import load_dotenv
import re
from urllib.parse import urlparse
from contextlib import asynccontextmanager
from pydantic_settings import BaseSettings

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
from data_query import router as data_query_router, db_pools

load_dotenv()


# Elastic indexes to use.
ES_LANDING_PAGE_SUMMARY = "landing-page-summary"
ES_TARGET_SUMMARY = "target-summary"
ES_DATASET_SUMMARY = "dataset-summary"


# Configuration
class Settings(BaseSettings):
    pg_host: str
    pg_port: int
    pg_user: str
    pg_password: str
    pg_db: str
    es_url: str
    es_username: str
    es_password: str


settings = Settings()


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup: Initialize connections
    db_pools["pg"] = await asyncpg.create_pool(
        user=settings.pg_user,
        password=settings.pg_password,
        database=settings.pg_db,
        host=settings.pg_host,
        port=settings.pg_port,
    )
    db_pools["es"] = AsyncElasticsearch(
        [settings.es_url], basic_auth=(settings.es_username, settings.es_password)
    )
    yield
    # Shutdown: Close connections
    await db_pools["pg"].close()
    await db_pools["es"].close()


app = FastAPI(title="Search API", version="1.0.0", lifespan=lifespan)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(data_query_router)


# Facet fields
FACET_FIELDS = [
    "license",
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

    # Text search across all searchable fields
    if query:
        cleaned_query = query.strip()
        if cleaned_query:
            # Exact/fuzzy matches with equal boost across all fields
            should_clauses.append(
                {
                    "multi_match": {
                        "query": cleaned_query,
                        "fields": [
                            "perturbed_target_symbol^1.5",
                            "perturbed_target_symbol.text^1.0",
                            "license^1.5",
                            "license.text^1.0",
                            "data_modalities^1.5",
                            "data_modalities.text^1.0",
                            "tissues_tested^1.5",
                            "tissues_tested.text^1.0",
                            "cell_types_tested^1.5",
                            "cell_types_tested.text^1.0",
                            "cell_lines_tested^1.5",
                            "cell_lines_tested.text^1.0",
                            "sex_tested^1.5",
                            "sex_tested.text^1.0",
                            "developmental_stages_tested^1.5",
                            "developmental_stages_tested.text^1.0",
                            "diseases_tested^1.5",
                            "diseases_tested.text^1.0",
                        ],
                        "type": "best_fields",
                        "fuzziness": "AUTO",
                    }
                }
            )

            # Prefix support for token beginnings (e.g. "SU" -> "SUMO1")
            searchable_fields = [
                "perturbed_target_symbol",
                "license",
                "data_modalities",
                "tissues_tested",
                "cell_types_tested",
                "cell_lines_tested",
                "sex_tested",
                "developmental_stages_tested",
                "diseases_tested",
            ]
            for field in searchable_fields:
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
                for field in searchable_fields:
                    should_clauses.append(
                        {
                            "wildcard": {
                                f"{field}.keyword": {
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
                                field: {
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


def _calculate_facets_from_results(results: List[Dict[str, Any]]) -> Facets:
    """Calculate facet counts from displayed results only.

    Used for search queries to ensure facets match the limited result set shown.
    """
    from collections import Counter

    facets_dict = {}
    for field in FACET_FIELDS:
        counter = Counter()
        for result in results:
            value = result.get(field)
            if value is None:
                continue
            if isinstance(value, list):
                for v in value:
                    if v:
                        counter[v] += 1
            elif value:
                counter[value] += 1

        facets_dict[field] = [
            FacetValue(value=val, count=count)
            for val, count in counter.most_common(100)
        ]
    return Facets(**facets_dict)


def parse_filters_from_params(
    license: Optional[str] = None,
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
    if license:
        filters["license"] = [v.strip() for v in license.split(",")]
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

    # When searching, limit to 20 best hits to show most relevant results
    # When browsing (no query), use normal pagination
    max_search_results = 20
    has_query = query and query.strip()

    if has_query:
        # For search queries, always return top 20 results (no pagination)
        effective_size = max_search_results
        from_ = 0
    else:
        # Normal pagination for browsing
        effective_size = size
        from_ = (page - 1) * size

    # Execute search with aggregations
    try:
        response = await db_pools["es"].search(
            index=ES_TARGET_SUMMARY,
            query=es_query,
            aggs=aggs,
            from_=from_,
            size=effective_size,
        )
    except Exception as e:
        error_detail = str(e)
        raise HTTPException(
            status_code=500, detail=f"Elasticsearch error: {error_detail}"
        )

    # Process results
    hits = response.get("hits", {})
    results = [hit["_source"] for hit in hits.get("hits", [])]

    if has_query:
        # For search queries: facets from displayed results only, capped total
        total = min(hits.get("total", {}).get("value", 0), max_search_results)
        total_pages = (total + effective_size - 1) // effective_size if total > 0 and effective_size > 0 else 0
        facets = _calculate_facets_from_results(results)
    else:
        # For browsing: use ES aggregations for full dataset facets
        total = hits.get("total", {}).get("value", 0)
        total_pages = (total + size - 1) // size if total > 0 else 0
        aggregations = response.get("aggregations", {})
        facets_dict = {}
        for field in FACET_FIELDS:
            buckets = aggregations.get(field, {}).get("buckets", [])
            facets_dict[field] = [
                FacetValue(value=bucket["key"], count=bucket["doc_count"])
                for bucket in buckets
            ]
        facets = Facets(**facets_dict)

    return SearchResponse(
        total=total,
        page=page,
        size=effective_size if has_query else size,
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
        if await db_pools["es"].ping():
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
            "host": urlparse(settings.es_url).hostname,
            "error": es_error,
        },
    }


@app.get("/summary", response_model=LandingPageSummary)
async def get_landing_page_summary():
    """
    Retrieve the landing page summary document from Elasticsearch.
    """
    try:
        response = await db_pools["es"].get(index=ES_LANDING_PAGE_SUMMARY, id="summary")
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
    license: Optional[str] = Query(
        None, description="Comma-separated list of licenses"
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
        license,
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


@app.get("/dataset/{dataset_id}")
async def get_dataset(dataset_id: str):
    """
    Retrieve a specific dataset record from Elasticsearch by dataset_id.
    """
    try:
        # Search for the dataset by dataset_id field
        response = await db_pools["es"].search(
            index=ES_DATASET_SUMMARY,
            query={"term": {"dataset_id": dataset_id}},
            size=1,
        )
    except Exception as exc:
        raise HTTPException(
            status_code=500, detail=f"Elasticsearch error: {exc}"
        ) from exc

    hits = response.get("hits", {}).get("hits", [])
    if not hits:
        raise HTTPException(
            status_code=404, detail=f"Dataset with id '{dataset_id}' not found"
        )

    return hits[0]["_source"]


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
