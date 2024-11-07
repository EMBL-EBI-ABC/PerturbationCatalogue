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
async def search(q: str = Query(None, description="Search query string")):
    # Access the Elasticsearch client from the app's state.
    es = app.state.es_client
    # Build the search query.
    if q:
        search_body = {
            "query": {"multi_match": {"query": q, "fields": ["*"]}},
            "size": 10,
        }
    else:
        search_body = {"query": {"match_all": {}}, "size": 10}

    try:
        # Execute the async search request.
        response = await es.search(index="mavedb", body=search_body)
        # Return the hits.
        return [r["_source"] for r in response["hits"]["hits"]]
    except Exception as e:
        # Handle ElasticSearch errors.
        raise HTTPException(status_code=500, detail=f"Search error: {str(e)}")
