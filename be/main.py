import os
import urllib.parse
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Query, Path
from elasticsearch import AsyncElasticsearch
from fastapi.middleware.cors import CORSMiddleware
from collections import defaultdict
from typing import Annotated

from models import (
    get_list_of_aggregations,
    ElasticResponse,
    ElasticDetailsResponse,
    MaveDBData,
    MaveDBSearchParams,
    MaveDBAggregationResponse
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
    }
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

async def elastic_search(params, index_name, data_class, aggregation_class):

    # Build the query body based on whether there is full text search.
    if params.q:
        query_body = {"multi_match": {"query": params.q, "fields": ["*"]}}
    else:
        query_body = {"match_all": {}}

    # Adding filters.
    filters = []
    aggregation_fields=get_list_of_aggregations(MaveDBAggregationResponse)
    if aggregation_fields:
        for aggregation_field in aggregation_fields:
            filter_value = getattr(params, aggregation_field)
            print(f"{aggregation_field}={filter_value}")
            if filter_value:
                filters.append(
                    {"terms": {aggregation_field: [filter_value]}})

    # Combine query with filters.
    search_body = {"from": params.start, "size": params.size, "query": {
        "bool": {
            "must": query_body,
            "filter": filters,
        }
    }, "aggs": defaultdict(dict)}

    # Adding aggregation fields.
    if aggregation_fields:
        for aggregation_field in aggregation_fields:
            search_body["aggs"][aggregation_field] = {
                "terms": {"field": aggregation_field}
            }

    # Adding sort field and sort order
    search_body['sort'] = [{params.sort_field: {"order": params.sort_order}}]

    # Performing the search.
    try:
        # Execute the async search request.
        response = await app.state.es_client.search(index=index_name, body=search_body)
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

async def elastic_details(index_name, record_id):
    try:
        response = await app.state.es_client.search(
            index=index_name,
            q=f"_id:{urllib.parse.quote(record_id)}"
        )
        hits = [r["_source"] for r in response["hits"]["hits"]]
        return hits
    except Exception as e:
        # Handle Elasticsearch errors.
        raise HTTPException(status_code=500, detail=f"Search error: {str(e)}")


# MaveDB.

@app.get("/mavedb/search")
async def mavedb_search(
    params: Annotated[MaveDBSearchParams, Query()]
    ) -> ElasticResponse[MaveDBData, MaveDBAggregationResponse]:
    return await elastic_search(
        params=params,
        index_name="mavedb",
        data_class=MaveDBData,
        aggregation_class=MaveDBAggregationResponse,
    )

@app.get("/mavedb/search/{record_id}")
async def mavedb_details(
    record_id: Annotated[str, Path(description="Record ID")]
    ) -> ElasticDetailsResponse[MaveDBData]:
    hits = await elastic_details("mavedb", record_id)
    return ElasticDetailsResponse[MaveDBData](results=hits)
