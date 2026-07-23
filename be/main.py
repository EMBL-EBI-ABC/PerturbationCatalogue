from fastapi import FastAPI, Query, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional, List, Dict, Any
from elasticsearch import AsyncElasticsearch
import asyncpg
from dotenv import load_dotenv
import json
import logging
import os
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
from target_search import build_target_fuzzy_query

load_dotenv()


# Elastic indexes to use.
ES_INDEX_SET = os.getenv("ES_INDEX_SET", "")
ES_TARGET_SUMMARY = f"target-summary{ES_INDEX_SET}"
ES_DATASET_SUMMARY = f"dataset-summary{ES_INDEX_SET}"
ES_LANDING_PAGE_SUMMARY = f"landing-page-summary{ES_INDEX_SET}"


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
logger = logging.getLogger(__name__)


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


# Facet fields for target-summary index
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

# Facet fields for dataset-summary index (actual ES field names)
DATASET_FACET_FIELDS = [
    "license_labels",
    "data_modalities",
    "tissue_labels",
    "cell_type_labels",
    "cell_line_labels",
    "sex_labels",
    "developmental_stage_labels",
    "disease_labels",
    "perturb_seq_reprocessed",
]

# Facet fields backed by an Elasticsearch boolean. Their aggregation buckets carry the
# value in `key_as_string` ("true"/"false") rather than a string `key`.
BOOLEAN_FACET_FIELDS = {"perturb_seq_reprocessed"}

# Mapping from canonical (target) field names to dataset index field names
TARGET_TO_DATASET_FIELD = {
    "license": "license_labels",
    "data_modalities": "data_modalities",
    "tissues_tested": "tissue_labels",
    "cell_types_tested": "cell_type_labels",
    "cell_lines_tested": "cell_line_labels",
    "sex_tested": "sex_labels",
    "developmental_stages_tested": "developmental_stage_labels",
    "diseases_tested": "disease_labels",
}
DATASET_TO_TARGET_FIELD = {v: k for k, v in TARGET_TO_DATASET_FIELD.items()}

# Searchable fields for dataset-summary index
DATASET_SEARCHABLE_FIELDS = [
    "dataset_id",
    "study_title",
    "experiment_title",
    "experiment_summary",
    "first_author",
    "last_author",
    "tissue_labels",
    "cell_type_labels",
    "cell_line_labels",
    "sex_labels",
    "developmental_stage_labels",
    "disease_labels",
    "license_labels",
    "library_perturbation_type_labels",
]


# Elasticsearch helper functions
def _escape_wildcard(value: str) -> str:
    """Escape characters that have special meaning in wildcard queries."""
    return re.sub(r"([\\*?])", r"\\\1", value)


def build_aggregations(facet_fields: Optional[List[str]] = None) -> Dict[str, Any]:
    """Build aggregations for all facet fields."""
    fields = facet_fields or FACET_FIELDS
    aggs = {}
    for field in fields:
        aggs[field] = {
            "terms": {
                "field": field,
                "size": 100,
            }
        }
    return aggs


def build_dataset_elasticsearch_query(
    query: Optional[str], filters: Optional[Dict[str, List[str]]]
) -> Dict[str, Any]:
    """Build Elasticsearch query for searching the dataset-summary index."""
    filter_clauses = []
    should_clauses = []

    if query:
        cleaned_query = query.strip()
        if cleaned_query:
            # Multi-match with boosted fields
            text_fields = []
            for field in DATASET_SEARCHABLE_FIELDS:
                text_fields.append(f"{field}^1.5")
                text_fields.append(f"{field}.text^1.0")
            should_clauses.append(
                {
                    "multi_match": {
                        "query": cleaned_query,
                        "fields": text_fields,
                        "type": "best_fields",
                        "fuzziness": "AUTO:5,8",
                    }
                }
            )

            # Prefix support
            for field in DATASET_SEARCHABLE_FIELDS:
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

            # Wildcard for partial/infix search
            wildcard_terms = []
            for term in cleaned_query.split():
                safe_term = _escape_wildcard(term.lower())
                if safe_term:
                    wildcard_terms.append(f"*{safe_term}*")

            if not wildcard_terms:
                safe_term = _escape_wildcard(cleaned_query.lower())
                if safe_term:
                    wildcard_terms.append(f"*{safe_term}*")

            for wildcard_value in wildcard_terms:
                for field in DATASET_SEARCHABLE_FIELDS:
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

    # Filters for dataset facet fields (map canonical target names to dataset names)
    if filters:
        for field, values in filters.items():
            dataset_field = TARGET_TO_DATASET_FIELD.get(field, field)
            if dataset_field in DATASET_FACET_FIELDS and values:
                filter_clauses.append({"terms": {dataset_field: values}})

    bool_query: Dict[str, Any] = {}

    if filter_clauses:
        bool_query["filter"] = filter_clauses

    if should_clauses:
        bool_query["should"] = should_clauses
        bool_query["minimum_should_match"] = 1

    if bool_query:
        return {"bool": bool_query}

    return {"match_all": {}}


def parse_filters_from_params(
    license: Optional[str] = None,
    data_modalities: Optional[str] = None,
    tissues_tested: Optional[str] = None,
    cell_types_tested: Optional[str] = None,
    cell_lines_tested: Optional[str] = None,
    sex_tested: Optional[str] = None,
    developmental_stages_tested: Optional[str] = None,
    diseases_tested: Optional[str] = None,
    # Dataset-mode filters
    license_labels: Optional[str] = None,
    library_perturbation_type_labels: Optional[str] = None,
    tissue_labels: Optional[str] = None,
    cell_type_labels: Optional[str] = None,
    cell_line_labels: Optional[str] = None,
    sex_labels: Optional[str] = None,
    developmental_stage_labels: Optional[str] = None,
    disease_labels: Optional[str] = None,
    perturb_seq_reprocessed: Optional[str] = None,
) -> Optional[Dict[str, List[str]]]:
    """Parse filters from query parameters"""
    filters = {}
    # Target-mode filters
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
    # Dataset-mode filters
    if license_labels:
        filters["license_labels"] = [v.strip() for v in license_labels.split(",")]
    if library_perturbation_type_labels:
        filters["library_perturbation_type_labels"] = [
            v.strip() for v in library_perturbation_type_labels.split(",")
        ]
    if tissue_labels:
        filters["tissue_labels"] = [v.strip() for v in tissue_labels.split(",")]
    if cell_type_labels:
        filters["cell_type_labels"] = [v.strip() for v in cell_type_labels.split(",")]
    if cell_line_labels:
        filters["cell_line_labels"] = [v.strip() for v in cell_line_labels.split(",")]
    if sex_labels:
        filters["sex_labels"] = [v.strip() for v in sex_labels.split(",")]
    if developmental_stage_labels:
        filters["developmental_stage_labels"] = [
            v.strip() for v in developmental_stage_labels.split(",")
        ]
    if disease_labels:
        filters["disease_labels"] = [v.strip() for v in disease_labels.split(",")]
    if perturb_seq_reprocessed:
        filters["perturb_seq_reprocessed"] = [
            v.strip() for v in perturb_seq_reprocessed.split(",")
        ]

    return filters if filters else None


async def perform_search(
    query: Optional[str],
    filters: Optional[Dict[str, List[str]]],
    page: int,
    size: int,
    search_mode: str = "targets",
    search_after: Optional[List[Any]] = None,
) -> SearchResponse:
    """Perform the actual search operation"""
    is_dataset_mode = search_mode == "datasets"

    # Build Elasticsearch query based on mode
    if is_dataset_mode:
        es_query = build_dataset_elasticsearch_query(query, filters)
        facet_fields = DATASET_FACET_FIELDS
        es_index = ES_DATASET_SUMMARY
        sort_field = "dataset_id"
    else:
        es_query = build_target_fuzzy_query(query, filters, FACET_FIELDS)
        facet_fields = FACET_FIELDS
        es_index = ES_TARGET_SUMMARY
        sort_field = "approved_symbol"

    aggs = build_aggregations(facet_fields)

    # Deterministic sort for search_after pagination support.
    # Always include _score first (1.0 for match_all) and a unique tiebreaker.
    sort = [{"_score": "desc"}, {sort_field: "asc"}]

    search_kwargs = {
        "index": es_index,
        "query": es_query,
        "aggs": aggs,
        "size": size,
        "sort": sort,
        "track_total_hits": True,
    }

    if search_after is not None:
        search_kwargs["search_after"] = search_after
    else:
        search_kwargs["from_"] = (page - 1) * size

    # Execute search with aggregations
    try:
        response = await db_pools["es"].search(**search_kwargs)
    except Exception as e:
        error_detail = str(e)
        raise HTTPException(
            status_code=500, detail=f"Elasticsearch error: {error_detail}"
        )

    # Process results
    hits = response.get("hits", {})
    hit_list = hits.get("hits", [])
    results = [hit["_source"] for hit in hit_list]

    total = hits.get("total", {}).get("value", 0)
    total_pages = (total + size - 1) // size if total > 0 else 0
    aggregations = response.get("aggregations", {})

    # Extract search_after cursor from last hit for deep pagination
    last_sort = hit_list[-1]["sort"] if hit_list else None

    # Build a lookup of lowercase → original-case values from result _source
    # fields, so we can restore proper display case for aggregation keys
    # (the lc_ascii normalizer lowercases all aggregation values).
    original_case: Dict[str, Dict[str, str]] = {}
    for field in facet_fields:
        field_map: Dict[str, str] = {}
        for result in results:
            val = result.get(field)
            if val is None:
                continue
            vals = val if isinstance(val, list) else [val]
            for v in vals:
                if v and isinstance(v, str):
                    field_map.setdefault(v.lower(), v)
        original_case[field] = field_map

    facets_dict = {}
    for field in facet_fields:
        buckets = aggregations.get(field, {}).get("buckets", [])
        case_map = original_case.get(field, {})
        if field in BOOLEAN_FACET_FIELDS:
            # Boolean aggregation buckets expose "true"/"false" via key_as_string.
            facets_dict[field] = [
                FacetValue(
                    value=bucket.get("key_as_string", str(bucket["key"])),
                    count=bucket["doc_count"],
                )
                for bucket in buckets
            ]
        else:
            facets_dict[field] = [
                FacetValue(
                    value=case_map.get(bucket["key"], bucket["key"]),
                    count=bucket["doc_count"],
                )
                for bucket in buckets
            ]

    # Remap dataset field names to canonical target field names
    if is_dataset_mode:
        facets_dict = {
            DATASET_TO_TARGET_FIELD.get(k, k): v for k, v in facets_dict.items()
        }

    facets = Facets(**facets_dict)

    return SearchResponse(
        total=total,
        page=page,
        size=size,
        total_pages=total_pages,
        results=results,
        facets=facets,
        search_after=last_sort,
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
        es_error = "Elasticsearch health check failed"
        logger.exception(
            "Elasticsearch health check failed: %s",
            {
                "status": "unhealthy",
                "elasticsearch": {
                    "status": es_status,
                    "host": urlparse(settings.es_url).hostname,
                    "error": str(e),
                },
            },
        )

    overall_status = "healthy" if es_status == "connected" else "unhealthy"

    return {
        "status": (
            overall_status + ". Check the logs"
            if overall_status == "unhealthy"
            else overall_status
        ),
        "elasticsearch": {
            "status": es_status,
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
        None, description="Search query for target symbol, synonym, name, or Ensembl ID"
    ),
    search_mode: str = Query(
        "targets", description="Search mode: 'targets' or 'datasets'"
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
    # Dataset-mode filter parameters
    license_labels: Optional[str] = Query(
        None, description="Comma-separated list of license labels (dataset mode)"
    ),
    library_perturbation_type_labels: Optional[str] = Query(
        None,
        description="Comma-separated list of perturbation type labels (dataset mode)",
    ),
    tissue_labels: Optional[str] = Query(
        None, description="Comma-separated list of tissue labels (dataset mode)"
    ),
    cell_type_labels: Optional[str] = Query(
        None, description="Comma-separated list of cell type labels (dataset mode)"
    ),
    cell_line_labels: Optional[str] = Query(
        None, description="Comma-separated list of cell line labels (dataset mode)"
    ),
    sex_labels: Optional[str] = Query(
        None, description="Comma-separated list of sex labels (dataset mode)"
    ),
    developmental_stage_labels: Optional[str] = Query(
        None,
        description="Comma-separated list of developmental stage labels (dataset mode)",
    ),
    disease_labels: Optional[str] = Query(
        None, description="Comma-separated list of disease labels (dataset mode)"
    ),
    perturb_seq_reprocessed: Optional[str] = Query(
        None,
        description="Comma-separated 'true'/'false' to filter by perturb-seq re-processing (dataset mode)",
    ),
    page: int = Query(1, ge=1, description="Page number (1-indexed)"),
    size: int = Query(6, ge=1, le=100, description="Number of results per page"),
    search_after: Optional[str] = Query(
        None,
        description="JSON-encoded search_after cursor for deep pagination beyond 10k results",
    ),
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
        license_labels,
        library_perturbation_type_labels,
        tissue_labels,
        cell_type_labels,
        cell_line_labels,
        sex_labels,
        developmental_stage_labels,
        disease_labels,
        perturb_seq_reprocessed,
    )
    parsed_search_after = None
    if search_after:
        try:
            parsed_search_after = json.loads(search_after)
        except (json.JSONDecodeError, TypeError):
            parsed_search_after = None
    return await perform_search(
        query, filters, page, size, search_mode, parsed_search_after
    )


@app.post("/search", response_model=SearchResponse)
async def search_post(request: SearchRequest):
    """
    Search endpoint (POST) with pagination, facets, filtering, and text search.
    Accepts JSON body with search parameters.
    """
    return await perform_search(
        request.query,
        request.filters,
        request.page,
        request.size,
        getattr(request, "search_mode", "targets"),
        request.search_after,
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
