import asyncio
from typing import Optional

from elasticsearch import AsyncElasticsearch
from fastapi import APIRouter, Depends, Request

router = APIRouter()


def get_es_client(request: Request) -> AsyncElasticsearch:
    return request.app.state.es_client


def get_pg_pool(request: Request):
    return request.app.state.pg_pool


@router.get("/v1/search")
async def search(
    dataset_metadata: Optional[str] = None,
    perturbation_gene_name: Optional[str] = None,
    effect_gene_name: Optional[str] = None,
    es_client: AsyncElasticsearch = Depends(get_es_client),
    pg_pool=Depends(get_pg_pool),
):
    es_task = asyncio.create_task(search_elastic(es_client, dataset_metadata))
    pg_task = asyncio.create_task(
        search_postgres(pg_pool, perturbation_gene_name, effect_gene_name)
    )

    es_results, pg_results = await asyncio.gather(es_task, pg_task)

    # Stitch results
    response = []
    for es_res in es_results:
        dataset_id = es_res["dataset_id"]
        perturbations = pg_results.get(dataset_id, [])
        es_res["perturbations"] = perturbations
        response.append(es_res)

    return response


async def search_elastic(
    es_client: AsyncElasticsearch, dataset_metadata: Optional[str]
):
    query = (
        {
            "bool": {
                "must": {
                    "multi_match": {
                        "query": dataset_metadata,
                        "fields": ["*"],
                    }
                }
            }
        }
        if dataset_metadata
        else {"match_all": {}}
    )

    es_response = await es_client.search(
        index="dataset-summary-v1",
        query=query,
        _source=[
            "dataset_id",
            "tissue_labels",
            "cell_type_labels.text",
            "cell_line_labels",
            "sex_labels",
            "developmental_stage_labels.text",
            "disease_labels",
            "library_perturbation_type_labels",
        ],
        size=1000,  # Assuming no more than 1000 datasets
        sort=[{"dataset_id": "asc"}],
    )
    results = []
    for hit in es_response["hits"]["hits"]:
        source = hit["_source"]
        result = {
            "dataset_id": source.get("dataset_id"),
            "tissue_label": (
                source.get("tissue_labels")[0] if source.get("tissue_labels") else None
            ),
            "cell_type_label": (
                source.get("cell_type_labels.text")[0]
                if source.get("cell_type_labels.text")
                else None
            ),
            "cell_line_label": (
                source.get("cell_line_labels")[0]
                if source.get("cell_line_labels")
                else None
            ),
            "sex_label": (
                source.get("sex_labels")[0] if source.get("sex_labels") else None
            ),
            "developmental_stage_label": (
                source.get("developmental_stage_labels.text")[0]
                if source.get("developmental_stage_labels.text")
                else None
            ),
            "disease_label": (
                source.get("disease_labels")[0]
                if source.get("disease_labels")
                else None
            ),
            "library_perturbation_type_label": (
                source.get("library_perturbation_type_labels")[0]
                if source.get("library_perturbation_type_labels")
                else None
            ),
        }
        results.append(result)
    return results


async def search_postgres(
    pg_pool, perturbation_gene_name: Optional[str], effect_gene_name: Optional[str]
):
    queries = []
    params = []
    if perturbation_gene_name:
        queries.append("perturbed_target_symbol = $1")
        params.append(perturbation_gene_name)
    if effect_gene_name:
        queries.append(f"gene = ${len(params) + 1}")
        params.append(effect_gene_name)

    where_clause = "WHERE " + " AND ".join(queries) if queries else ""

    query = f"""
        WITH ranked_perturbations AS (
            SELECT
                dataset_id,
                perturbed_target_symbol,
                gene,
                log2foldchange,
                padj,
                basemean,
                ROW_NUMBER() OVER(PARTITION BY dataset_id ORDER BY padj ASC) as rn
            FROM perturb_seq
            {where_clause}
        )
        SELECT
            dataset_id,
            perturbed_target_symbol,
            gene,
            log2foldchange,
            padj,
            basemean
        FROM ranked_perturbations
        WHERE rn <= 20
    """

    async with pg_pool.acquire() as connection:
        rows = await connection.fetch(query, *params)

    results = {}
    for row in rows:
        dataset_id = row["dataset_id"]
        if dataset_id not in results:
            results[dataset_id] = []
        results[dataset_id].append(
            {
                "perturbation_gene_name": row["perturbed_target_symbol"],
                "effect_gene_name": row["gene"],
                "log2foldchange": row["log2foldchange"],
                "padj": row["padj"],
                "base_mean": row["basemean"],
            }
        )
    return results
