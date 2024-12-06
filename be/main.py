import os
from contextlib import asynccontextmanager
from fastapi import FastAPI, HTTPException, Query
from elasticsearch import AsyncElasticsearch
from fastapi.middleware.cors import CORSMiddleware


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
    filter: str = Query(
        None,
        description="Filter query in the format field1:value1,value2+field2:value3,value4",
    ),
    start: int = Query(0, description="Starting point of the results"),
    size: int = Query(10, description="Number of results per page"),
):
    # Access the Elasticsearch client from the app's state.
    es = app.state.es_client

    # Parse filters from the filter query parameter.
    filters = []
    if filter:
        try:
            for f in filter.split("+"):
                field, values = f.split(":")
                values_list = values.split(",")
                filters.append({"terms": {f"{field}.keyword": values_list}})
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid filter format")

    # Build the query body based on whether there is full text search.
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
