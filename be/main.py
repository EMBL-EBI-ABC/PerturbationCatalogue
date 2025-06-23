import os
import urllib.parse
from contextlib import asynccontextmanager

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


async def elastic_search(index_name, params, data_class, aggregation_class):

    # Build the query body based on whether there is full text search.
    if params.q:
        query_body = {"multi_match": {"query": params.q, "fields": ["*"]}}
    else:
        query_body = {"match_all": {}}

    # Adding filters.
    filters = []
    aggregation_fields = get_list_of_aggregations(aggregation_class)
    if aggregation_fields:
        for aggregation_field in aggregation_fields:
            filter_value = getattr(params, aggregation_field)
            if filter_value:
                filters.append({"terms": {aggregation_field: [filter_value]}})

    # Combine query with filters.
    search_body = {
        "from": params.start,
        "size": params.size,
        "query": {
            "bool": {
                "must": query_body,
                "filter": filters,
            }
        },
        "aggs": defaultdict(dict),
    }

    # Adding aggregation fields.
    if aggregation_fields:
        for aggregation_field in aggregation_fields:
            search_body["aggs"][aggregation_field] = {
                "terms": {"field": aggregation_field, "size": 1000000}
            }

    # Adding sort field and sort order
    search_body["sort"] = [{params.sort_field: {"order": params.sort_order}}]

    # Performing the search.
    try:
        # Execute the async search request.
        response = await app.state.es_client.search(
            index=index_name, body=search_body, track_total_hits=True
        )
        # Extract total count and hits.
        total = response["hits"]["total"]["value"]
        hits = [r["_source"] for r in response["hits"]["hits"]]
        aggregations = response["aggregations"]
        # Return the results.
        return ElasticResponse[data_class, aggregation_class](
            total=total,
            start=params.start,
            size=params.size,
            results=hits,
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


async def get_unique_terms(index_name: str, aggs_body: dict) -> list[str]:
    """Performs an Elasticsearch aggregation to get unique terms from one or more fields."""
    search_body = {"size": 0, "aggs": aggs_body}
    try:
        response = await app.state.es_client.search(
            index=index_name, body=search_body, track_total_hits=False
        )

        unique_terms = set()
        aggregations = response["aggregations"]

        def extract_keys_from_buckets(agg_result):
            """Recursively traverses aggregation results to find and extract keys from buckets."""
            if "buckets" in agg_result:
                for bucket in agg_result["buckets"]:
                    unique_terms.add(bucket["key"])

            # Recurse into nested aggregations or sub-aggregations
            for value in agg_result.values():
                if isinstance(value, dict):
                    extract_keys_from_buckets(value)

        extract_keys_from_buckets(aggregations)

        return sorted(list(unique_terms))

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Aggregation error: {str(e)}")


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


# DepMap.


@app.get("/depmap/search")
async def depmap_search(
    params: Annotated[DepMapSearchParams, Query()],
) -> ElasticResponse[DepMapData, DepMapAggregationResponse]:
    return await elastic_search(
        index_name="depmap",
        params=params,
        data_class=DepMapData,
        aggregation_class=DepMapAggregationResponse,
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
