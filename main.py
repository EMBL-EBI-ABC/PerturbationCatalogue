import os
from contextlib import asynccontextmanager
from fastapi import FastAPI, HTTPException, Query
from elasticsearch import AsyncElasticsearch


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


@app.get("/search")
async def search(
    q: str = Query(None, description="Search query string"),
    start: int = Query(0, description="Starting point of the results"),
    size: int = Query(10, description="Number of results per page"),
):
    # Access the Elasticsearch client from the app's state.
    es = app.state.es_client
    # Build query body based on whether there is full text search.
    if q:
        query_body = {"multi_match": {"query": q, "fields": ["*"]}}
    else:
        query_body = {"match_all": {}}
    # Build the search query.
    search_body = {
        "from": start,
        "size": size,
        "query": query_body,
    }

    try:
        # Execute the async search request.
        response = await es.search(index="mavedb", body=search_body)
        # Extract total count and hits.
        total = response["hits"]["total"]["value"]
        hits = [r["_source"] for r in response["hits"]["hits"]]
        # Return the response with pagination info.
        return {"total": total, "start": start, "size": size, "results": hits}
    except Exception as e:
        # Handle Elasticsearch errors.
        raise HTTPException(status_code=500, detail=f"Search error: {str(e)}")
