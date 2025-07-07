import os
import urllib.parse
from contextlib import asynccontextmanager
import base64
import json
import asyncio

from fastapi import FastAPI, HTTPException, Query, Path
from elasticsearch import AsyncElasticsearch
from fastapi.middleware.cors import CORSMiddleware
from collections import defaultdict
from typing import Annotated, List

from models import (
    get_list_of_aggregations,
    ElasticResponse,
    ElasticDetailsResponse,
    MaveDBData,
    MaveDBSearchParams,
    MaveDBAggregationResponse,
    DepMapData,
    DepMapSearchParams,
    DepMapAggregationResponse,
    PerturbSeqData,
    PerturbSeqSearchParams,
    PerturbSeqAggregationResponse,
    Aggregation,
    RangeAggregation,
)


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Initialize AsyncElasticsearch.
    es_client = AsyncElasticsearch(
        [os.getenv("ES_URL")],
        http_auth=(os.getenv("ES_USERNAME"), os.getenv("ES_PASSWORD")),
        verify_certs=True,
    )
    # Pass the client to the app's state so it's accessible in routes.
    app.state.es_client = es_client
    yield
    # Clean up by closing the Elasticsearch client.
    await es_client.close()


# Initialize FastAPI with lifespan manager.
app = FastAPI(
    lifespan=lifespan,
    title="Perturbation Catalogue API",
    version="0.0.1",
    contact={
        "name": "Perturbation Catalogue Helpdesk",
        "email": "perturbation-catalogue-help@ebi.ac.uk",
    },
    license_info={
        "name": "Apache 2.0",
        "url": "https://www.apache.org/licenses/LICENSE-2.0.html",
    },
)

# Allow all origins.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allows all origins
    allow_credentials=True,
    allow_methods=["*"],  # Allows all HTTP methods
    allow_headers=["*"],  # Allows all headers
)


# Generic search methods.


async def elastic_search(
    index_name, params, data_class, aggregation_class, nested_query=None
):
    # Build the query body based on whether there is full text search.
    if params.q:
        # The base query will always search top-level fields.
        should_clauses = [{"multi_match": {"query": params.q, "fields": ["*"]}}]
        # If nested query info is provided, add a nested search clause.
        if nested_query and nested_query.get("path") and nested_query.get("field"):
            nested_path = nested_query["path"]
            nested_field = nested_query["field"]
            should_clauses.append(
                {
                    "nested": {
                        "path": nested_path,
                        "query": {
                            "match": {
                                # Construct the full field path: path.field
                                f"{nested_path}.{nested_field}": params.q
                            }
                        },
                    }
                }
            )
        # Combine clauses: a document matches if the term is in the top level OR the nested field.
        query_body = {"bool": {"should": should_clauses, "minimum_should_match": 1}}
    else:
        query_body = {"match_all": {}}

    # Adding filters.
    filters = []
    for field_name, model_field in aggregation_class.model_fields.items():
        filter_value = getattr(params, field_name, None)
        if not filter_value:
            continue

        # Check the type of aggregation to determine the filter logic.
        # The annotation will be either the Aggregation or RangeAggregation class.
        if model_field.annotation is Aggregation:  # List filter
            filters.append({"terms": {field_name: [filter_value]}})
        elif model_field.annotation is RangeAggregation:  # Slider/range filter
            try:
                min_val, max_val = filter_value.split("-", 1)
                range_query = {"gte": float(min_val), "lte": float(max_val)}
                filters.append({"range": {field_name: range_query}})
            except:
                raise HTTPException(
                    status_code=400,
                    detail=f"Invalid range format for '{field_name}'. Expected 'min-max' with numbers.",
                )

    # Base search body, without pagination fields yet.
    search_body = {
        "query": {"bool": {"must": query_body, "filter": filters}},
        "aggs": defaultdict(dict),
    }

    # Adding aggregation fields.
    for field_name, model_field in aggregation_class.model_fields.items():
        if model_field.annotation is Aggregation:
            search_body["aggs"][field_name] = {
                "terms": {"field": field_name, "size": 1000000}
            }
        elif model_field.annotation is RangeAggregation:
            # Use a stats aggregation to get min and max for the range.
            search_body["aggs"][field_name] = {"stats": {"field": field_name}}

    # Adding sort field and sort order. This must be deterministic for search_after.
    sort_config = [{params.sort_field: {"order": params.sort_order}}]
    # Add a tie-breaker for consistent sorting, required for search_after.
    if params.sort_field != "_id":
        sort_config.append({"_id": "asc"})
    search_body["sort"] = sort_config

    # Pagination Logic
    MAX_WINDOW = 10000

    if params.start + params.size <= MAX_WINDOW:
        # Use efficient 'from' for shallow pages
        search_body["from"] = params.start
    else:
        # Deep pagination using search_after.
        docs_to_skip = params.start
        last_hit_sort_values = None
        while docs_to_skip > 0:
            batch_size = min(docs_to_skip, MAX_WINDOW)
            preflight_body = {
                "query": search_body["query"],
                "sort": search_body["sort"],
                "size": batch_size,
                "_source": False,
                "track_total_hits": False,
            }
            if last_hit_sort_values:
                preflight_body["search_after"] = last_hit_sort_values
            try:
                preflight_response = await app.state.es_client.search(
                    index=index_name, body=preflight_body
                )
            except Exception as e:
                raise HTTPException(
                    status_code=500,
                    detail=f"Search error during deep pagination: {str(e)}",
                )

            preflight_hits = preflight_response["hits"]["hits"]
            if not preflight_hits:
                # The requested page is out of bounds. To return a valid response,
                # we run one final query to get total count and aggregations.
                agg_query_body = {
                    "query": search_body["query"],
                    "aggs": search_body["aggs"],
                    "size": 0,
                    "track_total_hits": True,
                }
                final_response = await app.state.es_client.search(
                    index=index_name, body=agg_query_body
                )
                return ElasticResponse[data_class, aggregation_class](
                    total=final_response["hits"]["total"]["value"],
                    start=params.start,
                    size=params.size,
                    results=[],
                    aggregations=final_response["aggregations"],
                )
            last_hit_sort_values = preflight_hits[-1].get("sort")
            docs_to_skip -= len(preflight_hits)
        if last_hit_sort_values:
            search_body["search_after"] = last_hit_sort_values

    search_body["size"] = params.size

    # Performing the search.
    try:
        # Execute the async search request.
        response = await app.state.es_client.search(
            index=index_name, body=search_body, track_total_hits=True
        )
        # Extract total count and hits.
        total = response["hits"]["total"]["value"]
        hits = response["hits"]["hits"]
        aggregations = response["aggregations"]

        # Return the results.
        return ElasticResponse[data_class, aggregation_class](
            total=total,
            start=params.start,
            size=params.size,
            results=[r["_source"] for r in hits],
            aggregations=aggregations,
        )

    except Exception as e:
        # Handle Elasticsearch errors.
        raise HTTPException(status_code=500, detail=f"Search error: {str(e)}")


async def elastic_details(index_name, record_id, data_class):
    try:
        # Quote the request, except for colons; they might appear in IDs.
        quoted_id = urllib.parse.quote(record_id).replace("%3A", ":")
        response = await app.state.es_client.search(
            index=index_name, query={"term": {"_id": quoted_id}}
        )
        hits = [r["_source"] for r in response["hits"]["hits"]]
        return ElasticDetailsResponse[data_class](results=hits)
    except Exception as e:
        # Handle Elasticsearch errors.
        raise HTTPException(status_code=500, detail=f"Search error: {str(e)}")


async def get_all_unique_terms_paginated(
    index_name: str,
    field: str,
    nested_path: str | None = None,
) -> set[str]:
    """Performs a paginated composite aggregation to retrieve all unique terms from a field.
    Can handle a single top-level or a single nested field."""
    unique_terms = set()
    after_key = None
    es_client = app.state.es_client
    sources = [{"term": {"terms": {"field": field}}}]
    while True:
        composite_agg = {"sources": sources, "size": 1000}
        if after_key:
            composite_agg["after"] = after_key
        agg_name = "paginated_terms"
        aggs_query = {agg_name: {"composite": composite_agg}}
        if nested_path:
            nested_agg_name = "nested_agg"
            aggs_query = {
                nested_agg_name: {
                    "nested": {"path": nested_path},
                    "aggs": aggs_query,
                }
            }
        search_body = {"size": 0, "aggs": aggs_query}
        try:
            response = await es_client.search(
                index=index_name, body=search_body, track_total_hits=False
            )
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Aggregation error: {str(e)}")
        agg_result = response["aggregations"]
        if nested_path:
            agg_result = agg_result.get(nested_agg_name, {})
        agg_result = agg_result.get(agg_name, {})
        buckets = agg_result.get("buckets", [])
        if not buckets:
            break
        for bucket in buckets:
            # The key in a composite agg is a dictionary of values.
            # Since we have one source, it will have one value.
            for value in bucket["key"].values():
                if value is not None:
                    unique_terms.add(str(value))
        if "after_key" in agg_result and agg_result["buckets"]:
            after_key = agg_result["after_key"]
        else:
            # No more pages
            break
    return unique_terms


# MaveDB.


@app.get("/mavedb/search")
async def mavedb_search(
    params: Annotated[MaveDBSearchParams, Query()],
) -> ElasticResponse[MaveDBData, MaveDBAggregationResponse]:
    return await elastic_search(
        index_name="mavedb",
        params=params,
        data_class=MaveDBData,
        aggregation_class=MaveDBAggregationResponse,
    )


@app.get("/mavedb/search/{record_id}")
async def mavedb_details(
    record_id: Annotated[str, Path(description="Record ID")],
) -> ElasticDetailsResponse[MaveDBData]:
    return await elastic_details(
        index_name="mavedb",
        record_id=record_id,
        data_class=MaveDBData,
    )


@app.get("/mavedb/genes", response_model=list[str])
async def mavedb_genes():
    """Returns the complete set of all genes which have any information for MaveDB."""
    terms = await get_all_unique_terms_paginated(
        index_name="mavedb", field="normalisedGeneName"
    )
    return sorted(list(terms))


# DepMap.


@app.get("/depmap/search")
async def depmap_search(
    params: Annotated[DepMapSearchParams, Query()],
) -> ElasticResponse[DepMapData, DepMapAggregationResponse]:
    nested_query_config = {"path": "high_dependency_genes", "field": "name"}
    return await elastic_search(
        index_name="depmap",
        params=params,
        data_class=DepMapData,
        aggregation_class=DepMapAggregationResponse,
        nested_query=nested_query_config,
    )


@app.get("/depmap/search/{record_id}")
async def depmap_details(
    record_id: Annotated[str, Path(description="Record ID")],
) -> ElasticDetailsResponse[DepMapData]:
    return await elastic_details(
        index_name="depmap",
        record_id=record_id,
        data_class=DepMapData,
    )


@app.get("/depmap/genes", response_model=list[str])
async def depmap_genes():
    """Returns the complete set of all genes which have any information for DepMap."""
    terms = await get_all_unique_terms_paginated(
        index_name="depmap",
        field="high_dependency_genes.name",
        nested_path="high_dependency_genes",
    )
    return sorted(list(terms))


# Perturb-Seq.


@app.get("/perturb-seq/search")
async def perturb_seq_search(
    params: Annotated[PerturbSeqSearchParams, Query()],
) -> ElasticResponse[PerturbSeqData, PerturbSeqAggregationResponse]:
    return await elastic_search(
        index_name="perturb-seq",
        params=params,
        data_class=PerturbSeqData,
        aggregation_class=PerturbSeqAggregationResponse,
    )


@app.get("/perturb-seq/search/{record_id}")
async def perturb_seq_details(
    record_id: Annotated[str, Path(description="Record ID")],
) -> ElasticDetailsResponse[PerturbSeqData]:
    return await elastic_details(
        index_name="perturb-seq",
        record_id=record_id,
        data_class=PerturbSeqData,
    )


@app.get("/perturb-seq/genes", response_model=list[str])
async def perturb_seq_genes():
    """Returns the complete set of all genes which have any information for Perturb-Seq."""
    perturbation_task = get_all_unique_terms_paginated(
        index_name="perturb-seq", field="perturbation"
    )
    gene_task = get_all_unique_terms_paginated(index_name="perturb-seq", field="gene")
    perturbation_terms, gene_terms = await asyncio.gather(perturbation_task, gene_task)
    all_terms = perturbation_terms.union(gene_terms)
    return sorted(list(all_terms))
