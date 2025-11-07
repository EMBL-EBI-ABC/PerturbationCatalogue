import asyncio
import os
from contextlib import asynccontextmanager

import asyncpg
from elasticsearch import AsyncElasticsearch
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import JSONResponse

# --- Environment Variables ---
ES_URL = os.getenv("ES_URL")
ES_USERNAME = os.getenv("ES_USERNAME")
ES_PASSWORD = os.getenv("ES_PASSWORD")
PG_HOST = os.getenv("PG_HOST_EXTERNAL")
PG_PORT = os.getenv("PG_PORT")
PG_USER = os.getenv("PG_USER")
PG_PASSWORD = os.getenv("PG_PASSWORD")
PG_DB = os.getenv("PG_DB")

# --- Global clients ---
es_client: AsyncElasticsearch | None = None
pg_pool: asyncpg.Pool | None = None


# --- FastAPI Lifespan ---
@asynccontextmanager
async def lifespan(app: FastAPI):
    global es_client, pg_pool
    es_client = AsyncElasticsearch(
        [ES_URL],
        basic_auth=(ES_USERNAME, ES_PASSWORD),
        verify_certs=True,
        request_timeout=60,
        max_retries=3,
        retry_on_timeout=True,
    )
    pg_pool = await asyncpg.create_pool(
        user=PG_USER,
        password=PG_PASSWORD,
        database=PG_DB,
        host=PG_HOST,
        port=PG_PORT,
    )
    yield
    await es_client.close()
    if pg_pool:
        await pg_pool.close()


app = FastAPI(lifespan=lifespan)

# --- Allowed parameters ---
ALLOWED_PARAMS = {
    "dataset_metadata",
    "perturbation_gene_name",
    "change_direction",
    "phenotype_gene_name",
    "modalities",
    "group_by",
    "max_datasets_per_modality",
    "max_top_level",
    "max_rows",
}


# --- Helper Functions ---
def validate_query_params(query_params):
    unrecognized_params = set(query_params.keys()) - ALLOWED_PARAMS
    if unrecognized_params:
        raise HTTPException(
            status_code=400,
            detail=f"Unrecognized query parameters: {', '.join(unrecognized_params)}",
        )


async def get_es_datasets(
    dataset_metadata: str | None,
    max_datasets_per_modality: int,
    dataset_ids: list[str] | None = None,
):
    es_query = {
        "size": 0,
        "query": {"bool": {"must": []}},
        "aggs": {
            "modalities": {
                "terms": {"field": "data_modalities", "size": 10},
                "aggs": {
                    "datasets": {
                        "top_hits": {
                            "size": max_datasets_per_modality,
                            "_source": [
                                "dataset_id",
                                "tissue_labels",
                                "cell_type_labels",
                                "cell_line_labels",
                                "sex_labels",
                                "developmental_stage_labels",
                                "disease_labels",
                                "library_perturbation_type_labels",
                            ],
                        }
                    }
                },
            },
            "tissue": {"terms": {"field": "tissue_labels", "size": 100}},
            "cell_type": {"terms": {"field": "cell_type_labels", "size": 100}},
            "cell_line": {"terms": {"field": "cell_line_labels", "size": 100}},
            "sex": {"terms": {"field": "sex_labels", "size": 100}},
            "developmental_stage": {
                "terms": {"field": "developmental_stage_labels", "size": 100}
            },
            "disease": {"terms": {"field": "disease_labels", "size": 100}},
            "library_perturbation_type": {
                "terms": {"field": "library_perturbation_type_labels", "size": 100}
            },
        },
    }

    if dataset_metadata:
        es_query["query"]["bool"]["must"].append(
            {"query_string": {"query": dataset_metadata}}
        )

    if dataset_ids:
        es_query["query"]["bool"]["must"].append({"terms": {"dataset_id": dataset_ids}})

    es_result = await es_client.search(index="dataset-summary-v2", body=es_query)

    datasets_by_modality = {}
    for bucket in es_result["aggregations"]["modalities"]["buckets"]:
        modality = bucket["key"].lower().replace(" ", "-")
        datasets = []
        for hit in bucket["datasets"]["hits"]["hits"]:
            source = hit["_source"]
            dataset = {"dataset_id": source["dataset_id"]}
            for key, value in source.items():
                if key.endswith("_labels"):
                    new_key = key.replace("_labels", "")
                    assert isinstance(value, list) and len(value) <= 1
                    dataset[new_key] = value[0] if value else None
            datasets.append(dataset)
        datasets_by_modality[modality] = datasets

    facet_counts = {}
    facet_fields = [
        "tissue",
        "cell_type",
        "cell_line",
        "sex",
        "developmental_stage",
        "disease",
        "library_perturbation_type",
    ]
    for facet in facet_fields:
        if facet in es_result["aggregations"]:
            facet_counts[facet] = [
                {"value": b["key"], "count": b["doc_count"]}
                for b in es_result["aggregations"][facet]["buckets"]
            ]

    total_counts = {
        b["key"].lower().replace(" ", "-"): b["doc_count"]
        for b in es_result["aggregations"]["modalities"]["buckets"]
    }

    return datasets_by_modality, facet_counts, total_counts


async def get_perturb_seq_dataset_ids(
    perturbation_gene_name: str | None,
    phenotype_gene_name: str | None,
    change_direction: str | None,
):
    async with pg_pool.acquire() as conn:
        query = "SELECT DISTINCT dataset_id FROM perturb_seq_2 WHERE 1=1"
        params = []
        param_idx = 1
        if perturbation_gene_name:
            query += f" AND perturbed_target_symbol = ${param_idx}"
            params.append(perturbation_gene_name)
            param_idx += 1
        if phenotype_gene_name:
            query += f" AND gene = ${param_idx}"
            params.append(phenotype_gene_name)
            param_idx += 1
        if change_direction:
            if change_direction == "increased":
                query += " AND log2foldchange > 0"
            elif change_direction == "decreased":
                query += " AND log2foldchange < 0"
        rows = await conn.fetch(query, *params)
        return {row["dataset_id"] for row in rows}


async def get_crispr_dataset_ids(perturbation_gene_name: str | None):
    async with pg_pool.acquire() as conn:
        query = "SELECT DISTINCT dataset_id FROM crispr_data WHERE 1=1"
        params = []
        if perturbation_gene_name:
            query += " AND perturbed_target_symbol = $1"
            params.append(perturbation_gene_name)
        rows = await conn.fetch(query, *params)
        return {row["dataset_id"] for row in rows}


async def get_mave_dataset_ids(perturbation_gene_name: str | None):
    async with pg_pool.acquire() as conn:
        query = "SELECT DISTINCT dataset_id FROM mave_data WHERE 1=1"
        params = []
        if perturbation_gene_name:
            query += " AND perturbed_target_symbol = $1"
            params.append(perturbation_gene_name)
        rows = await conn.fetch(query, *params)
        return {row["dataset_id"] for row in rows}


async def get_perturb_seq_data(
    datasets: list,
    group_by: str,
    perturbation_gene_name: str | None,
    phenotype_gene_name: str | None,
    change_direction: str | None,
    max_top_level: int,
    max_rows: int,
):
    if not datasets:
        return []

    results = []

    async with pg_pool.acquire() as conn:
        for dataset in datasets:
            dataset_id = dataset["dataset_id"]

            if group_by == "perturbation_gene_name":
                summary_view = "perturb_seq_summary_perturbation"
                group_col = "perturbed_target_symbol"
                detail_group_col = "perturbed_target_symbol"
                result_key = "by_perturbation"
                detail_key = "change_phenotype"

            else:  # phenotype_gene_name
                summary_view = "perturb_seq_summary_effect"
                group_col = "gene"
                detail_group_col = "gene"
                result_key = "by_phenotype"
                detail_key = "perturbation_change"

            summary_query = f"""
                SELECT {group_col}, n_total
                FROM {summary_view}
                WHERE dataset_id = $1
            """
            params = [dataset_id]
            param_idx = 2

            if perturbation_gene_name:
                summary_query += f" AND perturbed_target_symbol = ${param_idx}"
                params.append(perturbation_gene_name)
                param_idx += 1

            if phenotype_gene_name and group_by == "phenotype_gene_name":
                summary_query += f" AND gene = ${param_idx}"
                params.append(phenotype_gene_name)
                param_idx += 1

            summary_query += f" ORDER BY n_total DESC LIMIT ${param_idx}"
            params.append(max_top_level)

            top_level_entities = await conn.fetch(summary_query, *params)

            if not top_level_entities:
                continue

            dataset_result = {"dataset": dataset, result_key: []}

            for entity in top_level_entities:
                group_name = entity[group_col]

                detail_query = f"""
                    SELECT gene, log2foldchange, padj, basemean, perturbed_target_symbol
                    FROM perturb_seq_2
                    WHERE dataset_id = $1 AND {detail_group_col} = $2
                """
                detail_params = [dataset_id, group_name]
                param_idx_detail = 3

                if phenotype_gene_name and group_by == "perturbation_gene_name":
                    detail_query += f" AND gene = ${param_idx_detail}"
                    detail_params.append(phenotype_gene_name)
                    param_idx_detail += 1

                if change_direction:
                    if change_direction == "increased":
                        detail_query += " AND log2foldchange > 0"
                    elif change_direction == "decreased":
                        detail_query += " AND log2foldchange < 0"

                detail_query += f" ORDER BY padj ASC LIMIT ${param_idx_detail}"
                detail_params.append(max_rows)

                details = await conn.fetch(detail_query, *detail_params)

                group_items = []
                for row in details:
                    item = {
                        "change": {
                            "direction": (
                                "increased"
                                if row["log2foldchange"] > 0
                                else "decreased"
                            ),
                            "log2fc": row["log2foldchange"],
                            "padj": row["padj"],
                        }
                    }
                    if group_by == "perturbation_gene_name":
                        item["phenotype"] = {
                            "phenotype_gene_name": row["gene"],
                            "base_mean": row["basemean"],
                        }
                    else:  # phenotype_gene_name
                        item["perturbation"] = {
                            "perturbation_type": "gene_knockout",
                            "perturbation_gene_name": row["perturbed_target_symbol"],
                        }
                    group_items.append(item)

                if group_items:
                    top_level_item = {}
                    if group_by == "perturbation_gene_name":
                        top_level_item["perturbation"] = {
                            "perturbation_type": "gene_knockout",
                            "perturbation_gene_name": group_name,
                        }
                    else:  # phenotype_gene_name
                        top_level_item["phenotype"] = {
                            "phenotype_gene_name": group_name,
                        }
                    top_level_item[detail_key] = group_items
                    dataset_result[result_key].append(top_level_item)

            if dataset_result[result_key]:
                results.append(dataset_result)

    return results


async def get_crispr_data(
    datasets: list,
    perturbation_gene_name: str | None,
    max_rows: int,
):
    if not datasets:
        return []

    results = []
    async with pg_pool.acquire() as conn:
        for dataset in datasets:
            query = "SELECT perturbed_target_symbol, score_name, score_value FROM crispr_data WHERE dataset_id = $1"
            params = [dataset["dataset_id"]]
            if perturbation_gene_name:
                query += " AND perturbed_target_symbol = $2"
                params.append(perturbation_gene_name)

            query += f" LIMIT ${len(params) + 1}"
            params.append(max_rows)

            rows = await conn.fetch(query, *params)
            if rows:
                data = [
                    {
                        "perturbation": {
                            "perturbation_type": "gene_knockout",
                            "perturbation_gene_name": row["perturbed_target_symbol"],
                        },
                        "change": {"score_value": row["score_value"]},
                        "phenotype": {"score_name": row["score_name"]},
                    }
                    for row in rows
                ]
                results.append({"dataset": dataset, "data": data})
    return results


async def get_mave_data(
    datasets: list,
    perturbation_gene_name: str | None,
    max_rows: int,
):
    if not datasets:
        return []

    results = []
    async with pg_pool.acquire() as conn:
        for dataset in datasets:
            query = "SELECT perturbed_target_symbol, score_name, score_value, perturbation_name FROM mave_data WHERE dataset_id = $1"
            params = [dataset["dataset_id"]]
            if perturbation_gene_name:
                query += " AND perturbed_target_symbol = $2"
                params.append(perturbation_gene_name)

            query += f" LIMIT ${len(params) + 1}"
            params.append(max_rows)

            rows = await conn.fetch(query, *params)
            if rows:
                data = [
                    {
                        "perturbation": {
                            "perturbation_type": "mave",
                            "perturbation_gene_name": row["perturbed_target_symbol"],
                            "perturbation_name": row["perturbation_name"],
                        },
                        "change": {"score_value": row["score_value"]},
                        "phenotype": {"score_name": row["score_name"]},
                    }
                    for row in rows
                ]
                results.append({"dataset": dataset, "data": data})
    return results


# --- API Endpoint ---
@app.get("/v1/search")
async def search(
    request: Request,
    dataset_metadata: str = Query(None),
    perturbation_gene_name: str = Query(None),
    change_direction: str = Query(None, pattern="^(increased|decreased)$"),
    phenotype_gene_name: str = Query(None),
    modalities: str = Query("perturb-seq,crispr-screen,mave"),
    group_by: str = Query(
        ..., pattern="^(perturbation_gene_name|phenotype_gene_name)$"
    ),
    max_datasets_per_modality: int = Query(10, ge=1),
    max_top_level: int = Query(10, ge=1),
    max_rows: int = Query(10, ge=1),
):
    validate_query_params(request.query_params)

    target_modalities = {m.strip() for m in modalities.split(",")}

    # 1. Get dataset IDs from Postgres
    pg_tasks = []
    if "perturb-seq" in target_modalities:
        pg_tasks.append(
            get_perturb_seq_dataset_ids(
                perturbation_gene_name, phenotype_gene_name, change_direction
            )
        )
    if "crispr-screen" in target_modalities:
        pg_tasks.append(get_crispr_dataset_ids(perturbation_gene_name))
    if "mave" in target_modalities:
        pg_tasks.append(get_mave_dataset_ids(perturbation_gene_name))

    pg_results = await asyncio.gather(*pg_tasks)
    all_dataset_ids = list(set().union(*pg_results))

    # 2. Get datasets from ElasticSearch
    es_datasets, facet_counts, total_counts = await get_es_datasets(
        dataset_metadata, max_datasets_per_modality, all_dataset_ids
    )

    # 3. Get data from Postgres for each modality in parallel
    tasks = []
    if "perturb-seq" in target_modalities:
        tasks.append(
            get_perturb_seq_data(
                es_datasets.get("perturb-seq", []),
                group_by,
                perturbation_gene_name,
                phenotype_gene_name,
                change_direction,
                max_top_level,
                max_rows,
            )
        )
    if "crispr-screen" in target_modalities:
        tasks.append(
            get_crispr_data(
                es_datasets.get("crispr-screen", []),
                perturbation_gene_name,
                max_rows,
            )
        )
    if "mave" in target_modalities:
        tasks.append(
            get_mave_data(
                es_datasets.get("mave", []),
                perturbation_gene_name,
                max_rows,
            )
        )

    results = await asyncio.gather(*tasks)

    # 4. Assemble the final response
    response_modalities = []
    result_idx = 0

    all_modalities = ["perturb-seq", "crispr-screen", "mave"]

    for modality_name in all_modalities:
        if modality_name in target_modalities:
            response_modalities.append(
                {
                    "modality": modality_name,
                    "total_datasets_count": total_counts.get(modality_name, 0),
                    "datasets": (
                        results[result_idx] if result_idx < len(results) else []
                    ),
                }
            )
            result_idx += 1
        else:
            response_modalities.append(
                {
                    "modality": modality_name,
                    "total_datasets_count": total_counts.get(modality_name, 0),
                    "datasets": [],
                }
            )

    return JSONResponse(
        content={"modalities": response_modalities, "facet_counts": facet_counts}
    )
