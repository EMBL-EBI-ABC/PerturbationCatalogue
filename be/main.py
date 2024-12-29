import os
import urllib.parse
from contextlib import asynccontextmanager
from fastapi import FastAPI, HTTPException, Query
from elasticsearch import AsyncElasticsearch
from fastapi.middleware.cors import CORSMiddleware
from collections import defaultdict

from constants import AGGREGATION_FIELDS


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
app = FastAPI(lifespan=lifespan)

# Allow all origins.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allows all origins
    allow_credentials=True,
    allow_methods=["*"],  # Allows all HTTP methods
    allow_headers=["*"],  # Allows all headers
)


@app.get("/search")
async def search(
    q: str = Query(None, description="Search query string"),
    data_filter: str = Query(
        None,
        description="Filter query in the format field1:value1,value2+field2:value3,value4",
    ),
    start: int = Query(0, description="Starting point of the results"),
    size: int = Query(10, description="Number of results per page"),
    search_field: str = Query('publicationYear', description="Search field"),
    search_order: str = Query('desc', description="Search order"),
):
    # Access the Elasticsearch client from the app's state.
    es = app.state.es_client

    # Parse filters from the filter query parameter.
    filters = []
    if data_filter:
        try:
            for f in data_filter.split("+"):
                field, values = f.split(":")
                values_list = values.split(",")
                filters.append({"terms": {f"{field}": values_list}})
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid filter format")

    # Build the query body based on whether there is full text search.
    if q:
        query_body = {"multi_match": {"query": q, "fields": ["*"]}}
    else:
        query_body = {"match_all": {}}


    # Combine query with filters.
    search_body = {"from": start, "size": size, "query": {
        "bool": {
            "must": query_body,
            "filter": filters,
        }
    }, "aggs": defaultdict(dict)}

    # Adding aggregation fields
    for aggregation_field in AGGREGATION_FIELDS:
        search_body["aggs"][aggregation_field] = {
            "terms": {"field": aggregation_field}
        }

    search_body['sort'] = [{search_field: {"order": search_order}}]

    try:
        # Execute the async search request.
        response = await es.search(index="mavedb", body=search_body)
        # Extract total count and hits.
        total = response["hits"]["total"]["value"]
        hits = [r["_source"] for r in response["hits"]["hits"]]
        # Return the response with pagination info.
        return {
            "total": total,
            "start": start,
            "size": size,
            "results": hits,
            "aggregations": response["aggregations"]
        }
    except Exception as e:
        # Handle Elasticsearch errors.
        raise HTTPException(status_code=500, detail=f"Search error: {str(e)}")


@app.get("/search/{record_id}")
async def search(record_id: str):
    es = app.state.es_client
    try:
        response = await es.search(index="mavedb",
                                   q=f"_id:{urllib.parse.quote(record_id)}")
        hits = [r["_source"] for r in response["hits"]["hits"]]
        return {
            "results": hits
        }
    except Exception as e:
        # Handle Elasticsearch errors.
        raise HTTPException(status_code=500, detail=f"Search error: {str(e)}")