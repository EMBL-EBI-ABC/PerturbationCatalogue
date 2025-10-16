import asyncio
import os
from contextlib import asynccontextmanager

import asyncpg
from elasticsearch import AsyncElasticsearch
from fastapi import FastAPI, HTTPException, Request

# --- Database and Elasticsearch Connection ---

es_client = None
db_pool = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    """
    Asynchronous context manager to handle startup and shutdown of database connections.
    """
    global es_client, db_pool
    # Startup: Initialize connections
    es_client = AsyncElasticsearch(
        os.getenv("ES_URL"),
        basic_auth=(os.getenv("ES_USERNAME"), os.getenv("ES_PASSWORD")),
    )
    try:
        db_pool = await asyncpg.create_pool(
            host=os.getenv("PS_HOST"),
            port=os.getenv("PS_PORT"),
            user=os.getenv("PS_USER"),
            password=os.getenv("PS_PASSWORD"),
            database=os.getenv("PS_DB"),
            min_size=1,
            max_size=10,
        )
    except asyncpg.PostgresError as e:
        raise RuntimeError(f"Failed to connect to PostgreSQL: {e}") from e

    yield

    # Shutdown: Close connections
    await es_client.close()
    if db_pool:
        await db_pool.close()


app = FastAPI(lifespan=lifespan)

# --- API Endpoint ---

VALID_PARAMS = {
    "dataset_metadata",
    "perturbation_gene_name",
    "change_direction",
    "phenotype_gene_name",
    "group_by",
}

METADATA_FIELDS = [
    "tissue_labels",
    "cell_type_labels",
    "cell_line_labels",
    "sex_labels",
    "developmental_stage_labels",
    "disease_labels",
    "library_perturbation_type_labels",
]


@app.get("/v1/search")
async def search(request: Request):
    """
    Search endpoint to query perturbation data.
    """
    params = dict(request.query_params)
    validate_params(params)

    group_by = params.get("group_by")
    if not group_by:
        raise HTTPException(status_code=400, detail="group_by parameter is mandatory")

    # 1. Fetch and filter datasets from Elasticsearch
    datasets = await fetch_datasets(params)
    if not datasets:
        return []

    dataset_ids = [d["dataset_id"] for d in datasets]
    dataset_map = {d["dataset_id"]: d for d in datasets}

    # 2. Fetch data from Postgres based on grouping
    if group_by == "perturbation_gene_name":
        results = await get_by_perturbation(params, dataset_ids, dataset_map)
    elif group_by == "phenotype_gene_name":
        results = await get_by_phenotype(params, dataset_ids, dataset_map)
    else:
        raise HTTPException(status_code=400, detail="Invalid group_by value")

    return results


# --- Parameter Validation ---


def validate_params(params: dict):
    """
    Validates that only recognized query parameters are used.
    """
    unrecognized_params = set(params.keys()) - VALID_PARAMS
    if unrecognized_params:
        raise HTTPException(
            status_code=400,
            detail=f"Unrecognized parameters: {', '.join(unrecognized_params)}",
        )


# --- Elasticsearch Querying ---


async def fetch_datasets(params: dict) -> list[dict]:
    """
    Fetches dataset metadata from Elasticsearch, applying filters.
    """
    es_query = {"bool": {"must": []}}
    if query_str := params.get("dataset_metadata"):
        es_query["bool"]["must"].append(
            {"query_string": {"query": query_str, "fields": METADATA_FIELDS}}
        )

    resp = await es_client.search(
        index="dataset-summary-v1",
        query=es_query,
        source=["dataset_id", *METADATA_FIELDS],
        size=100,  # Max datasets to process
    )

    datasets = []
    for hit in resp["hits"]["hits"]:
        source = hit["_source"]
        dataset = {"dataset_id": source["dataset_id"]}
        for field in METADATA_FIELDS:
            if field in source and source[field]:
                # Flatten list and rename field as per spec
                assert len(source[field]) == 1, f"Expected single entry for {field}"
                dataset[field.removesuffix("_labels")] = source[field][0]
        datasets.append(dataset)
    return datasets


# --- PostgreSQL Querying ---


async def get_by_perturbation(params: dict, dataset_ids: list[str], dataset_map: dict):
    """
    Handles the case where data is grouped by perturbation.
    """
    summary_view = "perturb_seq_summary_perturbation"
    group_key = "perturbed_target_symbol"

    # Build filters for the summary view query
    filters = ["dataset_id = ANY($1)"]
    args = [dataset_ids]
    if p_gene := params.get("perturbation_gene_name"):
        filters.append(f"{group_key} = ${len(args) + 1}")
        args.append(p_gene)

    # Fetch top 3 perturbations per dataset
    query = f"""
        SELECT * FROM (
            SELECT dataset_id, {group_key} AS perturbation_gene_name, n_total, n_up, n_down,
                   ROW_NUMBER() OVER(PARTITION BY dataset_id ORDER BY n_total DESC) as rn
            FROM {summary_view}
            WHERE {" AND ".join(filters)}
        ) t WHERE rn <= 3
    """
    top_entities = await db_pool.fetch(query, *args)

    if not top_entities:
        return []

    # Fetch underlying change/phenotype data for each top perturbation
    details_tasks = [
        fetch_perturbation_details(entity, params) for entity in top_entities
    ]
    details_results = await asyncio.gather(*details_tasks)

    # Structure the final response
    return structure_response(
        top_entities,
        details_results,
        dataset_map,
        "by_perturbation",
        "perturbation",
        "change_phenotype",
    )


async def get_by_phenotype(params: dict, dataset_ids: list[str], dataset_map: dict):
    """
    Handles the case where data is grouped by phenotype.
    """
    summary_view = "perturb_seq_summary_effect"
    group_key = "gene"

    # Build filters for the summary view query
    filters = ["dataset_id = ANY($1)"]
    args = [dataset_ids]
    if ph_gene := params.get("phenotype_gene_name"):
        filters.append(f"{group_key} = ${len(args) + 1}")
        args.append(ph_gene)

    # Fetch top 3 phenotypes per dataset
    query = f"""
        SELECT * FROM (
            SELECT dataset_id, {group_key} AS phenotype_gene_name, n_total, n_up, n_down, base_mean,
                   ROW_NUMBER() OVER(PARTITION BY dataset_id ORDER BY n_total DESC) as rn
            FROM {summary_view}
            WHERE {" AND ".join(filters)}
        ) t WHERE rn <= 3
    """
    top_entities = await db_pool.fetch(query, *args)

    if not top_entities:
        return []

    # Fetch underlying perturbation/change data for each top phenotype
    details_tasks = [fetch_phenotype_details(entity, params) for entity in top_entities]
    details_results = await asyncio.gather(*details_tasks)

    return structure_response(
        top_entities,
        details_results,
        dataset_map,
        "by_phenotype",
        "phenotype",
        "perturbation_change",
    )


async def fetch_perturbation_details(entity, params):
    """
    Fetches the top 10 phenotype changes for a given perturbation.
    """
    filters = ["dataset_id = $1", "perturbed_target_symbol = $2"]
    args = [entity["dataset_id"], entity["perturbation_gene_name"]]

    if ph_gene := params.get("phenotype_gene_name"):
        filters.append(f"gene = ${len(args) + 1}")
        args.append(ph_gene)
    if direction := params.get("change_direction"):
        op = ">" if direction == "increased" else "<"
        filters.append(f"log2foldchange {op} 0")

    query = f"""
        SELECT log2foldchange, padj, gene as phenotype_gene_name, basemean
        FROM perturb_seq
        WHERE {" AND ".join(filters)}
        ORDER BY padj ASC
        LIMIT 10
    """
    async with db_pool.acquire() as conn:
        return await conn.fetch(query, *args)


async def fetch_phenotype_details(entity, params):
    """
    Fetches the top 10 perturbations for a given phenotype change.
    """
    filters = ["dataset_id = $1", "gene = $2"]
    args = [entity["dataset_id"], entity["phenotype_gene_name"]]

    if p_gene := params.get("perturbation_gene_name"):
        filters.append(f"perturbed_target_symbol = ${len(args) + 1}")
        args.append(p_gene)
    if direction := params.get("change_direction"):
        op = ">" if direction == "increased" else "<"
        filters.append(f"log2foldchange {op} 0")

    query = f"""
        SELECT perturbed_target_symbol as perturbation_gene_name, log2foldchange, padj
        FROM perturb_seq
        WHERE {" AND ".join(filters)}
        ORDER BY padj ASC
        LIMIT 10
    """
    async with db_pool.acquire() as conn:
        return await conn.fetch(query, *args)


# --- Response Formatting ---


def structure_response(
    top_entities, details_results, dataset_map, group_name, entity_name, details_name
):
    """
    Assembles the final nested JSON response from the query results.
    """
    # Organize details by their parent entity
    entity_details_map = {}
    for entity, details in zip(top_entities, details_results):
        key = (entity["dataset_id"], entity[f"{entity_name}_gene_name"])
        if key not in entity_details_map:
            entity_details_map[key] = {
                entity_name: dict(entity),
                details_name: format_details(details, entity_name),
            }

    # Group entities by dataset
    dataset_groups = {}
    for (dataset_id, _), data in entity_details_map.items():
        if dataset_id not in dataset_groups:
            dataset_groups[dataset_id] = []
        dataset_groups[dataset_id].append(data)

    # Build final list
    final_results = []
    for dataset_id, items in dataset_groups.items():
        final_results.append({"dataset": dataset_map[dataset_id], group_name: items})
    return final_results


def format_details(details, parent_entity_name):
    """
    Formats the innermost data rows into the required object structure.
    """
    formatted = []
    for row in details:
        if parent_entity_name == "perturbation":
            formatted.append(
                {
                    "change": {
                        "direction": (
                            "increased" if row["log2foldchange"] > 0 else "decreased"
                        ),
                        "log2fc": row["log2foldchange"],
                        "padj": row["padj"],
                    },
                    "phenotype": {
                        "phenotype_gene_name": row["phenotype_gene_name"],
                        "base_mean": row["basemean"],
                    },
                }
            )
        else:  # parent is phenotype
            formatted.append(
                {
                    "perturbation": {
                        "perturbation_gene_name": row["perturbation_gene_name"],
                    },
                    "change": {
                        "direction": (
                            "increased" if row["log2foldchange"] > 0 else "decreased"
                        ),
                        "log2fc": row["log2foldchange"],
                        "padj": row["padj"],
                    },
                }
            )
    return formatted
