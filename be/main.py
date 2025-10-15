import asyncio
import os
from functools import lru_cache
from typing import List, Optional

import asyncpg
from elasticsearch import AsyncElasticsearch
from fastapi import FastAPI, HTTPException, Query
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    es_url: str = "http://localhost:9200"
    es_username: str = "elastic"
    es_password: str = "changeme"
    ps_host: str = "localhost"
    ps_port: int = 5432
    ps_user: str = "user"
    ps_password: str = "password"
    ps_db: str = "perturb_seq"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


@lru_cache()
def get_settings():
    return Settings()


app = FastAPI()


@app.on_event("startup")
async def startup_event():
    settings = get_settings()
    app.state.es_client = AsyncElasticsearch(
        settings.es_url,
        basic_auth=(settings.es_username, settings.es_password),
        verify_certs=False,
    )
    try:
        app.state.pg_pool = await asyncpg.create_pool(
            host=settings.ps_host,
            port=settings.ps_port,
            user=settings.ps_user,
            password=settings.ps_password,
            database=settings.ps_db,
        )
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to connect to Postgres: {e}"
        )


@app.on_event("shutdown")
async def shutdown_event():
    await app.state.es_client.close()
    await app.state.pg_pool.close()


async def query_elasticsearch(
    dataset_metadata: Optional[str],
) -> List[str]:
    es_client = app.state.es_client
    if not dataset_metadata:
        return []

    search_fields = [
        "dataset_id",
        "tissue_labels",
        "cell_type_labels",
        "cell_line_labels",
        "sex_labels",
        "developmental_stage_labels",
        "disease_labels",
        "library_perturbation_type_labels",
    ]

    query = {
        "query": {
            "multi_match": {
                "query": dataset_metadata,
                "fields": search_fields,
            }
        },
        "_source": ["dataset_id"],
        "size": 1000,
    }

    response = await es_client.search(index="dataset-summary-v1", body=query)
    return [hit["_source"]["dataset_id"] for hit in response["hits"]["hits"]]


async def get_dataset_details(dataset_ids: List[str]) -> dict:
    es_client = app.state.es_client
    if not dataset_ids:
        return {}

    query = {"query": {"terms": {"dataset_id": dataset_ids}}}
    fields_to_return = [
        "dataset_id",
        "tissue_labels",
        "cell_type_labels",
        "cell_line_labels",
        "sex_labels",
        "developmental_stage_labels",
        "disease_labels",
        "library_perturbation_type_labels",
    ]
    query["_source"] = fields_to_return
    query["size"] = len(dataset_ids)

    response = await es_client.search(index="dataset-summary-v1", body=query)

    results = {}
    for hit in response["hits"]["hits"]:
        source = hit["_source"]
        dataset_id = source["dataset_id"]
        details = {"dataset_id": dataset_id}
        for key, value in source.items():
            if key.endswith("_labels"):
                singular_key = key.replace("_labels", "_label")
                if value:
                    assert (
                        len(value) == 1
                    ), f"Expected 1 value for {key}, got {len(value)}"
                    details[singular_key] = value[0]
                else:
                    details[singular_key] = None
        results[dataset_id] = details
    return results


@app.get("/v1/search")
async def search(
    dataset_metadata: Optional[str] = Query(None, alias="dataset"),
    perturbation_gene_name: Optional[str] = None,
    effect_gene_name: Optional[str] = None,
    group_by: Optional[str] = None,
):
    if group_by and group_by not in ["perturbation_gene_name", "effect_gene_name"]:
        raise HTTPException(status_code=400, detail="Invalid group_by value")

    es_dataset_ids = []
    if dataset_metadata:
        es_dataset_ids = await query_elasticsearch(dataset_metadata)
        if not es_dataset_ids:
            return []

    if group_by is None:
        async with app.state.pg_pool.acquire() as connection:
            results = await no_grouping(
                connection, es_dataset_ids, perturbation_gene_name, effect_gene_name
            )
    elif group_by == "perturbation_gene_name":
        results = await group_by_perturbation(
            app.state.pg_pool, es_dataset_ids, perturbation_gene_name, effect_gene_name
        )
    else:  # effect_gene_name
        results = await group_by_effect(
            app.state.pg_pool, es_dataset_ids, perturbation_gene_name, effect_gene_name
        )

    return results


async def no_grouping(
    connection, es_dataset_ids, perturbation_gene_name, effect_gene_name
):
    query = """
        SELECT dataset_id, perturbed_target_symbol as perturbation_gene_name, gene as effect_gene_name, basemean, log2foldchange, padj
        FROM perturb_seq
        WHERE 1=1
    """
    params = []
    if es_dataset_ids:
        params.append(tuple(es_dataset_ids))
        query += f" AND dataset_id = ANY(${len(params)})"
    if perturbation_gene_name:
        params.append(perturbation_gene_name)
        query += f" AND perturbed_target_symbol = ${len(params)}"
    if effect_gene_name:
        params.append(effect_gene_name)
        query += f" AND gene = ${len(params)}"

    # This is tricky. We need top 20 per dataset.
    # We can't just do a simple LIMIT. We need to use a window function.
    base_query = query
    query = f"""
        WITH ranked_data AS (
            {base_query.replace("SELECT", "SELECT ROW_NUMBER() OVER(PARTITION BY dataset_id ORDER BY padj ASC) as rn,")}
        )
        SELECT dataset_id, perturbation_gene_name, effect_gene_name, basemean, log2foldchange, padj
        FROM ranked_data
        WHERE rn <= 20
    """

    rows = await connection.fetch(query, *params)

    dataset_ids = list(set(row["dataset_id"] for row in rows))
    dataset_details = await get_dataset_details(dataset_ids)

    output = []
    for ds_id in dataset_ids:
        ds_details = dataset_details.get(ds_id)
        if not ds_details:
            continue

        data = [
            {
                "perturbation_gene_name": r["perturbation_gene_name"],
                "effect_gene_name": r["effect_gene_name"],
                "base_mean": r["basemean"],
                "log2foldchange": r["log2foldchange"],
                "padj": r["padj"],
            }
            for r in rows
            if r["dataset_id"] == ds_id
        ]
        if data:
            ds_details["data"] = data
            output.append(ds_details)

    return output


async def group_by_perturbation(
    pg_pool, es_dataset_ids, perturbation_gene_name, effect_gene_name
):
    # Get top 3 perturbations per dataset
    summary_query = """
        SELECT dataset_id, perturbed_target_symbol
        FROM (
            SELECT *, ROW_NUMBER() OVER(PARTITION BY dataset_id ORDER BY n_total DESC) as rn
            FROM perturb_seq_summary_perturbation
            WHERE 1=1
    """
    params = []
    if es_dataset_ids:
        params.append(tuple(es_dataset_ids))
        summary_query += f" AND dataset_id = ANY(${len(params)})"
    if perturbation_gene_name:
        params.append(perturbation_gene_name)
        summary_query += f" AND perturbed_target_symbol = ${len(params)}"

    summary_query += ") as ranked WHERE rn <= 3"

    async with pg_pool.acquire() as connection:
        top_perturbations = await connection.fetch(summary_query, *params)
    if not top_perturbations:
        return []

    dataset_to_perturbations = {}
    for row in top_perturbations:
        if row["dataset_id"] not in dataset_to_perturbations:
            dataset_to_perturbations[row["dataset_id"]] = []
        dataset_to_perturbations[row["dataset_id"]].append(
            row["perturbed_target_symbol"]
        )

    # Get top 10 effects for each of those perturbations
    tasks = []
    for dataset_id, perturbations in dataset_to_perturbations.items():
        task = get_effects_for_perturbations(
            pg_pool, dataset_id, perturbations, effect_gene_name
        )
        tasks.append(task)

    results_from_tasks = await asyncio.gather(*tasks)

    # results_from_tasks is a list of tuples (dataset_id, data)
    all_data = {}
    for dataset_id, data in results_from_tasks:
        if data:
            all_data[dataset_id] = data

    dataset_ids = list(all_data.keys())
    dataset_details = await get_dataset_details(dataset_ids)

    output = []
    for ds_id, data in all_data.items():
        ds_details = dataset_details.get(ds_id)
        if ds_details:
            ds_details["data"] = data
            output.append(ds_details)

    return output


async def get_effects_for_perturbations(
    pg_pool, dataset_id, perturbations, effect_gene_name
):
    query = """
        WITH ranked_effects AS (
            SELECT
                perturbed_target_symbol,
                gene as effect_gene_name,
                basemean,
                log2foldchange,
                padj,
                ROW_NUMBER() OVER(PARTITION BY perturbed_target_symbol ORDER BY padj ASC) as rn
            FROM perturb_seq
            WHERE dataset_id = $1 AND perturbed_target_symbol = ANY($2)
    """
    params = [dataset_id, perturbations]
    if effect_gene_name:
        query += " AND gene = $3"
        params.append(effect_gene_name)

    query += """
        )
        SELECT * FROM ranked_effects WHERE rn <= 10
    """

    async with pg_pool.acquire() as connection:
        effects_rows = await connection.fetch(query, *params)

        # Need summary data too
        summary_query = "SELECT perturbed_target_symbol, n_total, n_up, n_down FROM perturb_seq_summary_perturbation WHERE dataset_id = $1 AND perturbed_target_symbol = ANY($2)"
        summary_rows = await connection.fetch(summary_query, dataset_id, perturbations)

    summary_map = {r["perturbed_target_symbol"]: r for r in summary_rows}

    data = []
    for p_gene in perturbations:
        summary = summary_map.get(p_gene)
        if not summary:
            continue

        effects = [
            {
                "effect_gene_name": r["effect_gene_name"],
                "base_mean": r["basemean"],
                "log2foldchange": r["log2foldchange"],
                "padj": r["padj"],
            }
            for r in effects_rows
            if r["perturbed_target_symbol"] == p_gene
        ]
        if effects:
            data.append(
                {
                    "perturbation_gene_name": p_gene,
                    "n_total": summary["n_total"],
                    "n_up": summary["n_up"],
                    "n_down": summary["n_down"],
                    "effects": effects,
                }
            )
    return dataset_id, data


async def group_by_effect(
    pg_pool, es_dataset_ids, perturbation_gene_name, effect_gene_name
):
    # Get top 3 effects per dataset
    summary_query = """
        SELECT dataset_id, gene
        FROM (
            SELECT *, ROW_NUMBER() OVER(PARTITION BY dataset_id ORDER BY n_total DESC) as rn
            FROM perturb_seq_summary_effect
            WHERE 1=1
    """
    params = []
    if es_dataset_ids:
        params.append(tuple(es_dataset_ids))
        summary_query += f" AND dataset_id = ANY(${len(params)})"
    if effect_gene_name:
        params.append(effect_gene_name)
        summary_query += f" AND gene = ${len(params)}"

    summary_query += ") as ranked WHERE rn <= 3"

    async with pg_pool.acquire() as connection:
        top_effects = await connection.fetch(summary_query, *params)
    if not top_effects:
        return []

    dataset_to_effects = {}
    for row in top_effects:
        if row["dataset_id"] not in dataset_to_effects:
            dataset_to_effects[row["dataset_id"]] = []
        dataset_to_effects[row["dataset_id"]].append(row["gene"])

    # Get top 10 perturbations for each of those effects
    tasks = []
    for dataset_id, effects in dataset_to_effects.items():
        task = get_perturbations_for_effects(
            pg_pool, dataset_id, effects, perturbation_gene_name
        )
        tasks.append(task)

    results_from_tasks = await asyncio.gather(*tasks)

    all_data = {}
    for dataset_id, data in results_from_tasks:
        if data:
            all_data[dataset_id] = data

    dataset_ids = list(all_data.keys())
    dataset_details = await get_dataset_details(dataset_ids)

    output = []
    for ds_id, data in all_data.items():
        ds_details = dataset_details.get(ds_id)
        if ds_details:
            ds_details["data"] = data
            output.append(ds_details)

    return output


async def get_perturbations_for_effects(
    pg_pool, dataset_id, effects, perturbation_gene_name
):
    query = """
        WITH ranked_perturbations AS (
            SELECT
                gene,
                perturbed_target_symbol as perturbation_gene_name,
                log2foldchange,
                padj,
                ROW_NUMBER() OVER(PARTITION BY gene ORDER BY padj ASC) as rn
            FROM perturb_seq
            WHERE dataset_id = $1 AND gene = ANY($2)
    """
    params = [dataset_id, effects]
    if perturbation_gene_name:
        query += " AND perturbed_target_symbol = $3"
        params.append(perturbation_gene_name)

    query += """
        )
        SELECT * FROM ranked_perturbations WHERE rn <= 10
    """

    async with pg_pool.acquire() as connection:
        perturbation_rows = await connection.fetch(query, *params)

        # Need summary data too
        summary_query = "SELECT gene, n_total, n_up, n_down, base_mean FROM perturb_seq_summary_effect WHERE dataset_id = $1 AND gene = ANY($2)"
        summary_rows = await connection.fetch(summary_query, dataset_id, effects)

    summary_map = {r["gene"]: r for r in summary_rows}

    data = []
    for e_gene in effects:
        summary = summary_map.get(e_gene)
        if not summary:
            continue

        perturbations = [
            {
                "perturbation_gene_name": r["perturbation_gene_name"],
                "log2foldchange": r["log2foldchange"],
                "padj": r["padj"],
            }
            for r in perturbation_rows
            if r["gene"] == e_gene
        ]
        if perturbations:
            data.append(
                {
                    "effect_gene_name": e_gene,
                    "n_total": summary["n_total"],
                    "n_up": summary["n_up"],
                    "n_down": summary["n_down"],
                    "base_mean": summary["base_mean"],
                    "perturbations": perturbations,
                }
            )
    return dataset_id, data
