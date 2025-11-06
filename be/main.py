import asyncio
import os
from contextlib import asynccontextmanager
from typing import List, Dict, Any, Literal, Optional

import asyncpg
from elasticsearch import AsyncElasticsearch
from fastapi import FastAPI, HTTPException, Request
from pydantic import BaseModel, Field
from pydantic_settings import BaseSettings


# --- Settings and Configuration ---


class Settings(BaseSettings):
    es_url: str = "http://localhost:9200"
    es_username: str = "elastic"
    es_password: str = ""
    pg_host: str = "localhost"
    pg_port: int = 5432
    pg_user: str = "postgres"
    pg_password: str = ""
    pg_db: str = "postgres"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


settings = Settings()

# --- Global Clients ---

clients = {}


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup: Initialize clients
    clients["es"] = AsyncElasticsearch(
        hosts=[settings.es_url], basic_auth=(settings.es_username, settings.es_password)
    )
    try:
        clients["pg"] = await asyncpg.create_pool(
            host=settings.pg_host,
            port=settings.pg_port,
            user=settings.pg_user,
            password=settings.pg_password,
            database=settings.pg_db,
        )
        yield
    finally:
        # Shutdown: Close clients
        if "pg" in clients and clients["pg"]:
            await clients["pg"].close()
        if "es" in clients and clients["es"]:
            await clients["es"].close()


app = FastAPI(lifespan=lifespan)

# --- Pydantic Models for Response Structure ---


# Dataset
class Dataset(BaseModel):
    dataset_id: str
    tissue: Optional[str] = None
    cell_type: Optional[str] = None
    cell_line: Optional[str] = None
    sex: Optional[str] = None
    developmental_stage: Optional[str] = None
    disease: Optional[str] = None
    library_perturbation_type: Optional[str] = None


# Perturb-seq Models
class PerturbSeqChange(BaseModel):
    direction: str
    log2fc: float = Field(..., alias="log2foldchange")
    padj: float


class PerturbSeqPerturbation(BaseModel):
    perturbation_type: str = "gene_knockout"
    perturbation_gene_name: str = Field(..., alias="perturbed_target_symbol")
    n_total: int
    n_up: int
    n_down: int


class PerturbSeqPhenotype(BaseModel):
    phenotype_gene_name: str = Field(..., alias="gene")
    base_mean: Optional[float] = None
    n_total: Optional[int] = None
    n_down: Optional[int] = None
    n_up: Optional[int] = None


class ChangePhenotype(BaseModel):
    change: PerturbSeqChange
    phenotype: PerturbSeqPhenotype


class PerturbationChange(BaseModel):
    perturbation: PerturbSeqPerturbation
    change: PerturbSeqChange


class ByPerturbation(BaseModel):
    perturbation: PerturbSeqPerturbation
    change_phenotype: List[ChangePhenotype]


class ByPhenotype(BaseModel):
    phenotype: PerturbSeqPhenotype
    perturbation_change: List[PerturbationChange]


class PerturbSeqDatasetData(BaseModel):
    dataset: Dataset
    by_perturbation: Optional[List[ByPerturbation]] = None
    by_phenotype: Optional[List[ByPhenotype]] = None


# CRISPR/MAVE Models
class GenericPerturbation(BaseModel):
    perturbation_type: str
    perturbation_gene_name: str = Field(..., alias="perturbed_target_symbol")
    perturbation_name: Optional[str] = None


class GenericChange(BaseModel):
    score_value: float


class GenericPhenotype(BaseModel):
    score_name: str


class GenericDataRow(BaseModel):
    perturbation: GenericPerturbation
    change: GenericChange
    phenotype: GenericPhenotype


class GenericDatasetData(BaseModel):
    dataset: Dataset
    data: List[GenericDataRow]


# Top-level Response Models
class ModalityData(BaseModel):
    modality: str
    datasets: List[Any]  # List[PerturbSeqDatasetData] or List[GenericDatasetData]


# --- Helper Functions ---


def validate_query_params(request: Request):
    allowed_params = {
        "dataset_metadata",
        "perturbation_gene_name",
        "change_direction",
        "phenotype_gene_name",
        "group_by",
    }
    unrecognized_params = [
        key for key in request.query_params if key not in allowed_params
    ]
    if unrecognized_params:
        raise HTTPException(
            status_code=400,
            detail=f"Unrecognized query parameters: {', '.join(unrecognized_params)}",
        )


async def get_datasets_from_es(
    es_client: AsyncElasticsearch,
    dataset_metadata: Optional[str] = None,
    dataset_ids: Optional[set] = None,
) -> Dict[str, List[Dict]]:
    query_fields = [
        "dataset_id",
        "tissue_labels",
        "cell_type_labels",
        "cell_line_labels",
        "sex_labels",
        "developmental_stage_labels",
        "disease_labels",
        "library_perturbation_type_labels",
    ]

    query = {"bool": {"must": []}}

    if dataset_metadata:
        query["bool"]["must"].append(
            {"multi_match": {"query": dataset_metadata, "fields": query_fields}}
        )

    if dataset_ids:
        query["bool"]["must"].append({"terms": {"dataset_id": list(dataset_ids)}})

    if not query["bool"]["must"]:
        query = {"match_all": {}}

    response = await es_client.search(
        index="dataset-summary-v2",
        query=query,
        source=query_fields + ["data_modalities"],
        size=1000,  # Fetch a large number to be filtered and limited in code
    )

    datasets_by_modality = {}
    for hit in response["hits"]["hits"]:
        source = hit["_source"]
        modality_val = source.get("data_modalities")
        if not modality_val:
            continue

        # Per spec, data_modalities is a single value, but ES might return a list.
        modality = modality_val[0] if isinstance(modality_val, list) else modality_val

        # Limit to 10 datasets per modality
        if len(datasets_by_modality.get(modality, [])) >= 10:
            continue

        dataset_info = {
            "dataset_id": source.get("dataset_id"),
            "data_modalities": modality,
            "tissue": (
                (v[0] if v else None)
                if isinstance((v := source.get("tissue_labels")), list)
                else v
            ),
            "cell_type": (
                (v[0] if v else None)
                if isinstance((v := source.get("cell_type_labels")), list)
                else v
            ),
            "cell_line": (
                (v[0] if v else None)
                if isinstance((v := source.get("cell_line_labels")), list)
                else v
            ),
            "sex": (
                (v[0] if v else None)
                if isinstance((v := source.get("sex_labels")), list)
                else v
            ),
            "developmental_stage": (
                (v[0] if v else None)
                if isinstance((v := source.get("developmental_stage_labels")), list)
                else v
            ),
            "disease": (
                (v[0] if v else None)
                if isinstance((v := source.get("disease_labels")), list)
                else v
            ),
            "library_perturbation_type": (
                (v[0] if v else None)
                if isinstance(
                    (v := source.get("library_perturbation_type_labels")), list
                )
                else v
            ),
        }

        if modality not in datasets_by_modality:
            datasets_by_modality[modality] = []
        datasets_by_modality[modality].append(dataset_info)

    return datasets_by_modality


async def get_perturb_seq_data(
    pool: asyncpg.Pool,
    datasets: List[Dict],
    group_by: str,
    perturbation_gene_name: Optional[str],
    phenotype_gene_name: Optional[str],
    change_direction: Optional[str],
) -> List[PerturbSeqDatasetData]:

    dataset_ids = [d["dataset_id"] for d in datasets]
    if not dataset_ids:
        return []

    # Determine summary view and grouping keys
    if group_by == "perturbation_gene_name":
        summary_view = "perturb_seq_summary_perturbation"
        group_key = "perturbed_target_symbol"
        detail_key = "gene"
    else:  # phenotype_gene_name
        summary_view = "perturb_seq_summary_effect"
        group_key = "gene"
        detail_key = "perturbed_target_symbol"

    # 1. Get top 3 entities from the summary view
    params = [dataset_ids]
    param_idx = 2

    filter_clauses = ["dataset_id = ANY($1)"]
    if perturbation_gene_name:
        # In the new logic, datasets are already pre-filtered, but we still need to filter the summary view.
        if group_by == "perturbation_gene_name":
            filter_clauses.append(f"{group_key} = ${param_idx}")
            params.append(perturbation_gene_name)
            param_idx += 1
    if phenotype_gene_name:
        if group_by == "phenotype_gene_name":
            filter_clauses.append(f"{group_key} = ${param_idx}")
            params.append(phenotype_gene_name)
            param_idx += 1

    summary_query = f"""
        WITH ranked_entities AS (
            SELECT *, ROW_NUMBER() OVER(PARTITION BY dataset_id ORDER BY n_total DESC) as rn
            FROM {summary_view}
            WHERE {" AND ".join(filter_clauses)}
        )
        SELECT * FROM ranked_entities WHERE rn <= 3;
    """
    top_entities = await pool.fetch(summary_query, *params)

    if not top_entities:
        return []

    # 2. Get up to 10 detail rows for each top entity
    detail_queries = []
    for entity in top_entities:
        change_filter = ""
        if change_direction == "increased":
            change_filter = "AND log2foldchange > 0"
        elif change_direction == "decreased":
            change_filter = "AND log2foldchange < 0"

        # Detail query needs to filter by both grouping key and detail key if they are provided
        detail_filters = ["dataset_id = $1", f"{group_key} = $2"]
        detail_params = [entity["dataset_id"], entity[group_key]]
        detail_param_idx = 3

        if perturbation_gene_name:
            detail_filters.append(f"perturbed_target_symbol = ${detail_param_idx}")
            detail_params.append(perturbation_gene_name)
            detail_param_idx += 1
        if phenotype_gene_name:
            detail_filters.append(f"gene = ${detail_param_idx}")
            detail_params.append(phenotype_gene_name)
            detail_param_idx += 1

        detail_query_sql = f"""
            SELECT * FROM perturb_seq_2
            WHERE {" AND ".join(detail_filters)} {change_filter}
            ORDER BY padj ASC
            LIMIT 10;
        """
        detail_queries.append(pool.fetch(detail_query_sql, *detail_params))

    detail_results_list = await asyncio.gather(*detail_queries)

    # 3. Assemble the response structure
    datasets_map = {d["dataset_id"]: d for d in datasets}
    results_by_dataset = {}

    for i, entity in enumerate(top_entities):
        dataset_id = entity["dataset_id"]
        if dataset_id not in results_by_dataset:
            results_by_dataset[dataset_id] = {
                "dataset": Dataset(**datasets_map[dataset_id]),
                "by_perturbation": [] if group_by == "perturbation_gene_name" else None,
                "by_phenotype": [] if group_by == "phenotype_gene_name" else None,
            }

        detail_results = detail_results_list[i]
        if not detail_results:
            continue

        if group_by == "perturbation_gene_name":
            change_phenotypes = [
                ChangePhenotype(
                    change=PerturbSeqChange(
                        direction=(
                            "increased" if row["log2foldchange"] > 0 else "decreased"
                        ),
                        **row,
                    ),
                    phenotype=PerturbSeqPhenotype(**row),
                )
                for row in detail_results
            ]
            if change_phenotypes:  # Only add if there are matching changes
                results_by_dataset[dataset_id]["by_perturbation"].append(
                    ByPerturbation(
                        perturbation=PerturbSeqPerturbation(**entity),
                        change_phenotype=change_phenotypes,
                    )
                )
        else:  # by phenotype
            perturbation_changes = [
                PerturbationChange(
                    perturbation=PerturbSeqPerturbation(**row),
                    change=PerturbSeqChange(
                        direction=(
                            "increased" if row["log2foldchange"] > 0 else "decreased"
                        ),
                        **row,
                    ),
                )
                for row in detail_results
            ]
            if perturbation_changes:
                results_by_dataset[dataset_id]["by_phenotype"].append(
                    ByPhenotype(
                        phenotype=PerturbSeqPhenotype(**entity),
                        perturbation_change=perturbation_changes,
                    )
                )

    return [
        PerturbSeqDatasetData(**data)
        for data in results_by_dataset.values()
        if (data.get("by_perturbation") or data.get("by_phenotype"))
    ]


async def get_generic_data(
    pool: asyncpg.Pool,
    datasets: List[Dict],
    table_name: str,
    modality_type: str,
    perturbation_gene_name: Optional[str],
) -> List[GenericDatasetData]:

    dataset_ids = [d["dataset_id"] for d in datasets]
    if not dataset_ids:
        return []

    params = [dataset_ids]
    param_idx = 2

    filter_clauses = ["dataset_id = ANY($1)"]
    if perturbation_gene_name:
        filter_clauses.append(f"perturbed_target_symbol = ${param_idx}")
        params.append(perturbation_gene_name)

    query = f"""
        WITH ranked_rows AS (
            SELECT *, ROW_NUMBER() OVER(PARTITION BY dataset_id ORDER BY perturbed_target_symbol) as rn
            FROM {table_name}
            WHERE {" AND ".join(filter_clauses)}
        )
        SELECT * FROM ranked_rows WHERE rn <= 10;
    """

    rows = await pool.fetch(query, *params)

    datasets_map = {d["dataset_id"]: d for d in datasets}
    results_by_dataset = {}

    for row in rows:
        dataset_id = row["dataset_id"]
        if dataset_id not in results_by_dataset:
            results_by_dataset[dataset_id] = {
                "dataset": Dataset(**datasets_map[dataset_id]),
                "data": [],
            }

        perturbation_type = (
            "gene_knockout" if modality_type == "CRISPR screen" else "mave"
        )

        results_by_dataset[dataset_id]["data"].append(
            GenericDataRow(
                perturbation=GenericPerturbation(
                    perturbation_type=perturbation_type, **row
                ),
                change=GenericChange(**row),
                phenotype=GenericPhenotype(**row),
            )
        )

    return [GenericDatasetData(**data) for data in results_by_dataset.values()]


# --- API Endpoint ---


@app.get("/v1/search", response_model=List[ModalityData])
async def search(
    request: Request,
    group_by: Literal["perturbation_gene_name", "phenotype_gene_name"],
    dataset_metadata: Optional[str] = None,
    perturbation_gene_name: Optional[str] = None,
    change_direction: Optional[Literal["increased", "decreased"]] = None,
    phenotype_gene_name: Optional[str] = None,
):
    validate_query_params(request)

    es_client = clients["es"]
    pg_pool = clients["pg"]

    # If any gene/change filters are applied, pre-filter dataset IDs from Postgres
    pg_filtered_ids = None
    if perturbation_gene_name or phenotype_gene_name or change_direction:
        tasks = []
        # Perturb-seq
        ps_q = "SELECT DISTINCT dataset_id FROM perturb_seq_2 WHERE 1=1"
        ps_p = []
        if perturbation_gene_name:
            ps_p.append(perturbation_gene_name)
            ps_q += f" AND perturbed_target_symbol = ${len(ps_p)}"
        if phenotype_gene_name:
            ps_p.append(phenotype_gene_name)
            ps_q += f" AND gene = ${len(ps_p)}"
        if change_direction == "increased":
            ps_q += " AND log2foldchange > 0"
        elif change_direction == "decreased":
            ps_q += " AND log2foldchange < 0"
        tasks.append(pg_pool.fetch(ps_q, *ps_p))

        # CRISPR & MAVE (only filter by perturbation_gene_name)
        if perturbation_gene_name and not phenotype_gene_name:
            tasks.append(
                pg_pool.fetch(
                    "SELECT DISTINCT dataset_id FROM crispr_data WHERE perturbed_target_symbol = $1",
                    perturbation_gene_name,
                )
            )
            tasks.append(
                pg_pool.fetch(
                    "SELECT DISTINCT dataset_id FROM mave_data WHERE perturbed_target_symbol = $1",
                    perturbation_gene_name,
                )
            )

        list_of_results = await asyncio.gather(*tasks)
        pg_filtered_ids = {
            row["dataset_id"] for result_list in list_of_results for row in result_list
        }

        if not pg_filtered_ids:
            return []

    datasets_by_modality = await get_datasets_from_es(
        es_client, dataset_metadata, pg_filtered_ids
    )

    tasks = []

    # Perturb-seq
    if "Perturb-seq" in datasets_by_modality:
        tasks.append(
            get_perturb_seq_data(
                pg_pool,
                datasets_by_modality["Perturb-seq"],
                group_by,
                perturbation_gene_name,
                phenotype_gene_name,
                change_direction,
            )
        )

    # CRISPR screen
    if "CRISPR screen" in datasets_by_modality and not phenotype_gene_name:
        tasks.append(
            get_generic_data(
                pg_pool,
                datasets_by_modality["CRISPR screen"],
                "crispr_data",
                "CRISPR screen",
                perturbation_gene_name,
            )
        )

    # MAVE
    if "MAVE" in datasets_by_modality and not phenotype_gene_name:
        tasks.append(
            get_generic_data(
                pg_pool,
                datasets_by_modality["MAVE"],
                "mave_data",
                "MAVE",
                perturbation_gene_name,
            )
        )

    results_from_tasks = await asyncio.gather(*tasks)

    final_response = []
    task_idx = 0

    if "Perturb-seq" in datasets_by_modality:
        if (
            results_from_tasks
            and task_idx < len(results_from_tasks)
            and results_from_tasks[task_idx]
        ):
            final_response.append(
                ModalityData(
                    modality="perturb-seq", datasets=results_from_tasks[task_idx]
                )
            )
        task_idx += 1

    if "CRISPR screen" in datasets_by_modality and not phenotype_gene_name:
        if (
            results_from_tasks
            and task_idx < len(results_from_tasks)
            and results_from_tasks[task_idx]
        ):
            final_response.append(
                ModalityData(
                    modality="crispr-screen", datasets=results_from_tasks[task_idx]
                )
            )
        task_idx += 1

    if "MAVE" in datasets_by_modality and not phenotype_gene_name:
        if (
            results_from_tasks
            and task_idx < len(results_from_tasks)
            and results_from_tasks[task_idx]
        ):
            final_response.append(
                ModalityData(modality="mave", datasets=results_from_tasks[task_idx])
            )
        task_idx += 1

    return final_response
