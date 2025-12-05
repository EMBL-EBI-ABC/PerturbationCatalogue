import os
import re
from typing import Any, Dict, List, Literal, Optional, Tuple

import asyncpg
from fastapi import APIRouter, Depends, HTTPException, Query, Request
from pydantic import BaseModel, Field


# --- Database Connection Management ---
db_pools: Dict[str, Any] = {}


router = APIRouter()

# --- Constants and Mappings ---
MODALITIES = Literal["perturb-seq", "crispr-screen", "mave"]

PG_TABLES = {
    "perturb-seq": "perturb_seq_data",
    "crispr-screen": "crispr_data",
    "mave": "mave_data",
}

# Field mappings from API to database
PERTURB_SEQ_PG_MAPPING = {
    "perturbation_gene_name": "perturbed_target_symbol",
    "effect_gene_name": "gene",
    "effect_log2fc": "log2foldchange",
    "effect_padj": "padj",
    "effect_base_mean": "basemean",
}
CRISPR_PG_MAPPING = {
    "perturbation_gene_name": "perturbed_target_symbol",
    "effect_score_name": "score_name",
    "effect_score_value": "score_value",
    "effect_significant": "significant",
    "effect_significance_criteria": "significance_criteria",
}
MAVE_PG_MAPPING = {
    "perturbation_gene_name": "perturbed_target_symbol",
    "perturbation_name": "perturbation_name",
    "effect_score_name": "score_name",
    "effect_score_value": "score_value",
}
PG_MAPPINGS = {
    "perturb-seq": PERTURB_SEQ_PG_MAPPING,
    "crispr-screen": CRISPR_PG_MAPPING,
    "mave": MAVE_PG_MAPPING,
}

ELASTIC_FIELD_MAPPING = {
    "dataset_id": "dataset_id",
    "dataset_tissue": "tissue_labels",
    "dataset_cell_type": "cell_type_labels",
    "dataset_cell_line": "cell_line_labels",
    "dataset_sex": "sex_labels",
    "dataset_developmental_stage": "developmental_stage_labels",
    "dataset_disease": "disease_labels",
    "dataset_library_perturbation_type": "library_perturbation_type_labels",
    "dataset_license_id": "license_ids",
    "dataset_license_label": "license_labels",
    "dataset_score_interpretation": "score_interpretation",
}

# --- Pydantic Models for API Response ---


# Perturbation Models
class PerturbationBase(BaseModel):
    gene_name: str = Field(..., alias="perturbation_gene_name")


class MavePerturbation(PerturbationBase):
    name: Optional[str] = Field(None, alias="perturbation_name")


class PerturbSeqPerturbation(PerturbationBase):
    n_total: int = Field(..., alias="perturbation_n_total")
    n_up: int = Field(..., alias="perturbation_n_up")
    n_down: int = Field(..., alias="perturbation_n_down")


# Effect Models
class EffectBase(BaseModel):
    pass


class PerturbSeqEffect(EffectBase):
    gene_name: str = Field(..., alias="effect_gene_name")
    direction: str = Field(..., alias="effect_direction")
    log2fc: float = Field(..., alias="effect_log2fc")
    padj: float = Field(..., alias="effect_padj")
    base_mean: float = Field(..., alias="effect_base_mean")
    n_total: int = Field(..., alias="effect_n_total")
    n_up: int = Field(..., alias="effect_n_up")
    n_down: int = Field(..., alias="effect_n_down")


class ScoreEffect(EffectBase):
    score_name: str = Field(..., alias="effect_score_name")
    score_value: float = Field(..., alias="effect_score_value")


# Result Models
class Result(BaseModel):
    perturbation: Dict
    effect: Dict


# Dataset Models
class DatasetMetadata(BaseModel):
    id: str = Field(..., alias="dataset_id")
    tissue: Optional[str] = Field(None, alias="dataset_tissue")
    cell_type: Optional[str] = Field(None, alias="dataset_cell_type")
    cell_line: Optional[str] = Field(None, alias="dataset_cell_line")
    sex: Optional[str] = Field(None, alias="dataset_sex")
    developmental_stage: Optional[str] = Field(
        None, alias="dataset_developmental_stage"
    )
    disease: Optional[str] = Field(None, alias="dataset_disease")
    library_perturbation_type: Optional[str] = Field(
        None, alias="dataset_library_perturbation_type"
    )
    license_id: Optional[str] = Field(None, alias="dataset_license_id")
    license_label: Optional[str] = Field(None, alias="dataset_license_label")
    score_interpretation: Optional[str] = Field(
        None, alias="dataset_score_interpretation"
    )


class DatasetResult(BaseModel):
    dataset: DatasetMetadata
    results: List[Result]


# Facet Models
class FacetValue(BaseModel):
    value: Any
    count: int


# Top-level Response Models
class ModalitySearchResponse(BaseModel):
    total_datasets_count: int
    facet_counts: Dict[str, List[FacetValue]]
    datasets: List[DatasetResult]


class DatasetSearchResponse(BaseModel):
    total_rows_count: int
    offset: int
    limit: int
    results: List[Result]


# --- Dependency Classes for Query Parameters ---


class CommonModalitySearchParams:
    def __init__(
        self,
        dataset_metadata: Optional[str] = Query(
            None, description="Search term for dataset metadata"
        ),
        dataset_tissue: Optional[str] = Query(None, description="Filter by tissue"),
        dataset_cell_type: Optional[str] = Query(
            None, description="Filter by cell type"
        ),
        dataset_cell_line: Optional[str] = Query(
            None, description="Filter by cell line"
        ),
        dataset_sex: Optional[str] = Query(None, description="Filter by sex"),
        dataset_developmental_stage: Optional[str] = Query(
            None, description="Filter by developmental stage"
        ),
        dataset_disease: Optional[str] = Query(None, description="Filter by disease"),
        dataset_library_perturbation_type: Optional[str] = Query(
            None, description="Filter by library perturbation type"
        ),
        dataset_license_id: Optional[str] = Query(
            None, description="Filter by license ID"
        ),
        dataset_license_label: Optional[str] = Query(
            None, description="Filter by license label"
        ),
        dataset_score_interpretation: Optional[str] = Query(
            None, description="Filter by score interpretation"
        ),
        dataset_limit: int = Query(10, description="Number of datasets to return"),
        dataset_offset: int = Query(0, description="Offset for datasets"),
        rows_per_dataset_limit: int = Query(
            10, description="Number of rows per dataset to return"
        ),
        sort: Optional[str] = Query(
            None, description="Sort order (e.g., 'field:asc,other:desc')"
        ),
    ):
        self.dataset_metadata = dataset_metadata
        self.dataset_tissue = dataset_tissue
        self.dataset_cell_type = dataset_cell_type
        self.dataset_cell_line = dataset_cell_line
        self.dataset_sex = dataset_sex
        self.dataset_developmental_stage = dataset_developmental_stage
        self.dataset_disease = dataset_disease
        self.dataset_library_perturbation_type = dataset_library_perturbation_type
        self.dataset_license_id = dataset_license_id
        self.dataset_license_label = dataset_license_label
        self.dataset_score_interpretation = dataset_score_interpretation
        self.dataset_limit = dataset_limit
        self.dataset_offset = dataset_offset
        self.rows_per_dataset_limit = rows_per_dataset_limit
        self.sort = sort

    def dict(self):
        return {k: v for k, v in self.__dict__.items() if v is not None}


class CommonDatasetSearchParams:
    def __init__(
        self,
        limit: int = Query(50, description="Number of rows to return"),
        offset: int = Query(0, description="Offset for rows"),
        sort: Optional[str] = Query(
            None, description="Sort order (e.g., 'field:asc,other:desc')"
        ),
    ):
        self.limit = limit
        self.offset = offset
        self.sort = sort

    def dict(self):
        return {k: v for k, v in self.__dict__.items() if v is not None}


class MaveParams:
    def __init__(
        self,
        perturbation_gene_name: Optional[str] = Query(
            None, description="Filter by perturbation gene name"
        ),
        perturbation_name: Optional[str] = Query(
            None, description="Filter by perturbation name"
        ),
        effect_score_name: Optional[str] = Query(
            None, description="Filter by effect score name"
        ),
        effect_score_value: Optional[str] = Query(
            None,
            description="Filter by effect score value (supports ranges e.g., '1_10', '1_', '_10')",
        ),
    ):
        self.perturbation_gene_name = perturbation_gene_name
        self.perturbation_name = perturbation_name
        self.effect_score_name = effect_score_name
        self.effect_score_value = effect_score_value

    def dict(self):
        return {k: v for k, v in self.__dict__.items() if v is not None}


class CrisprScreenParams:
    def __init__(
        self,
        perturbation_gene_name: Optional[str] = Query(
            None, description="Filter by perturbation gene name"
        ),
        effect_score_name: Optional[str] = Query(
            None, description="Filter by effect score name"
        ),
        effect_score_value: Optional[str] = Query(
            None, description="Filter by effect score value (supports ranges)"
        ),
        effect_significant: Optional[str] = Query(
            None, description="Filter by effect significant (true/false)"
        ),
        effect_significance_criteria: Optional[str] = Query(
            None, description="Filter by effect significance criteria"
        ),
    ):
        self.perturbation_gene_name = perturbation_gene_name
        self.effect_score_name = effect_score_name
        self.effect_score_value = effect_score_value
        self.effect_significant = effect_significant
        self.effect_significance_criteria = effect_significance_criteria

    def dict(self):
        return {k: v for k, v in self.__dict__.items() if v is not None}


class PerturbSeqParams:
    def __init__(
        self,
        perturbation_gene_name: Optional[str] = Query(
            None, description="Filter by perturbation gene name"
        ),
        effect_gene_name: Optional[str] = Query(
            None, description="Filter by effect gene name"
        ),
        effect_log2fc: Optional[str] = Query(
            None, description="Filter by effect log2fc (supports ranges)"
        ),
        effect_padj: Optional[str] = Query(
            None, description="Filter by effect padj (supports ranges)"
        ),
        effect_base_mean: Optional[str] = Query(
            None, description="Filter by effect base mean (supports ranges)"
        ),
    ):
        self.perturbation_gene_name = perturbation_gene_name
        self.effect_gene_name = effect_gene_name
        self.effect_log2fc = effect_log2fc
        self.effect_padj = effect_padj
        self.effect_base_mean = effect_base_mean

    def dict(self):
        return {k: v for k, v in self.__dict__.items() if v is not None}


# --- Helper Functions ---


def parse_numeric_filter(param_name: str, value: str) -> Tuple[str, List[Any]]:
    """Parses numeric filter syntax (e.g., 1_10, 1_, _10) into SQL."""
    if "_" in value:
        min_val, max_val = value.split("_", 1)
        conditions = []
        params = []
        if min_val:
            conditions.append(f"{param_name} >= $... ")
            params.append(float(min_val))
        if max_val:
            conditions.append(f"{param_name} <= $... ")
            params.append(float(max_val))
        return " AND ".join(conditions), params
    else:
        return f"{param_name} = $... ", [float(value)]


def get_api_to_db_mapping(modality: MODALITIES) -> Dict[str, str]:
    """Returns the combined API to DB field mapping for a modality."""
    return PG_MAPPINGS.get(modality, {})


def validate_query_params(
    query_params: Dict[str, Any], modality: MODALITIES, dataset_id: Optional[str] = None
):
    """Validates that all query params are known for the endpoint."""
    valid_params = {
        "dataset_metadata",
        "dataset_tissue",
        "dataset_cell_type",
        "dataset_cell_line",
        "dataset_sex",
        "dataset_developmental_stage",
        "dataset_disease",
        "dataset_library_perturbation_type",
        "dataset_license_id",
        "dataset_license_label",
        "dataset_score_interpretation",
        "sort",
        "dataset_limit",
        "dataset_offset",
        "rows_per_dataset_limit",
        "limit",
        "offset",
    }
    if dataset_id:
        valid_params = {"sort", "limit", "offset"}

    # Add all filterable perturbation and effect fields to valid_params
    pg_mapping = get_api_to_db_mapping(modality)
    valid_params.update(pg_mapping.keys())

    for param in query_params:
        if param not in valid_params:
            raise HTTPException(
                status_code=400, detail=f"Invalid query parameter: {param}"
            )


async def enrich_perturb_seq_rows(
    conn: asyncpg.Connection, dataset_id: str, rows: List[Dict]
) -> List[Dict]:
    """Enriches perturb-seq rows with data from summary views."""
    if not rows:
        return []

    # Fetch perturbation summaries
    pert_keys = list(
        set((row["dataset_id"], row["perturbed_target_symbol"]) for row in rows)
    )
    pert_summary_map = {}
    if pert_keys:
        pert_dataset_ids = [k[0] for k in pert_keys]
        pert_symbols = [k[1] for k in pert_keys]
        pert_summary_rows = await conn.fetch(
            """
            SELECT t.dataset_id, t.perturbed_target_symbol, t.n_total, t.n_up, t.n_down
            FROM perturb_seq_summary_perturbation AS t
            JOIN unnest($1::text[], $2::text[]) AS keys(did, pts)
            ON t.dataset_id = keys.did AND t.perturbed_target_symbol = keys.pts
            """,
            pert_dataset_ids,
            pert_symbols,
        )
        for r in pert_summary_rows:
            pert_summary_map[(r["dataset_id"], r["perturbed_target_symbol"])] = r

    # Fetch effect summaries
    effect_keys = list(set((row["dataset_id"], row["gene"]) for row in rows))
    effect_summary_map = {}
    if effect_keys:
        effect_dataset_ids = [k[0] for k in effect_keys]
        effect_genes = [k[1] for k in effect_keys]
        effect_summary_rows = await conn.fetch(
            """
            SELECT t.dataset_id, t.gene, t.n_total, t.n_up, t.n_down
            FROM perturb_seq_summary_effect AS t
            JOIN unnest($1::text[], $2::text[]) AS keys(did, g)
            ON t.dataset_id = keys.did AND t.gene = keys.g
            """,
            effect_dataset_ids,
            effect_genes,
        )
        for r in effect_summary_rows:
            effect_summary_map[(r["dataset_id"], r["gene"])] = r

    # Enrich rows
    for row in rows:
        pert_summary = pert_summary_map.get(
            (row["dataset_id"], row["perturbed_target_symbol"]), {}
        )
        row["perturbation_n_total"] = pert_summary.get("n_total")
        row["perturbation_n_up"] = pert_summary.get("n_up")
        row["perturbation_n_down"] = pert_summary.get("n_down")

        effect_summary = effect_summary_map.get((row["dataset_id"], row["gene"]), {})
        row["effect_n_total"] = effect_summary.get("n_total")
        row["effect_n_up"] = effect_summary.get("n_up")
        row["effect_n_down"] = effect_summary.get("n_down")

        log2fc = row.get("log2foldchange")
        if log2fc is None:
            row["effect_direction"] = "not available"
        elif log2fc > 0:
            row["effect_direction"] = "increased"
        elif log2fc < 0:
            row["effect_direction"] = "decreased"
        else:
            row["effect_direction"] = "no change"

    return rows


# --- Shared Implementation Functions ---


async def _search_modality_impl(
    modality: MODALITIES,
    query_params: Dict[str, Any],
):
    """Search across all datasets within a modality (Shared Implementation)."""
    dataset_limit = query_params.get("dataset_limit", 10)
    dataset_offset = query_params.get("dataset_offset", 0)
    rows_per_dataset_limit = query_params.get("rows_per_dataset_limit", 10)
    sort = query_params.get("sort")

    validate_query_params(query_params, modality)

    pg_conn = db_pools["pg"]
    es_client = db_pools["es"]

    pg_table = PG_TABLES[modality]
    api_to_db = get_api_to_db_mapping(modality)

    # 1. Pre-filter Datasets (Postgres)
    pg_filters = []
    pg_params: List[Any] = []

    for key, value in query_params.items():
        if key in api_to_db:
            db_field = api_to_db[key]
            # Simple string filter
            if isinstance(value, str) and "_" not in value:
                pg_filters.append(f"{db_field} = ${len(pg_params) + 1}")
                pg_params.append(value)
            # Numeric range filter
            elif isinstance(value, str):
                condition, params = parse_numeric_filter(db_field, value)
                # This is a bit tricky because parse_numeric_filter doesn't know the param index
                condition = condition.replace("$...", f"${len(pg_params) + 1}")
                if " AND " in condition:
                    condition = condition.replace("$...", f"${len(pg_params) + 2}")
                pg_filters.append(condition)
                pg_params.extend(params)

    where_clause = f"WHERE {' AND '.join(pg_filters)}" if pg_filters else ""

    prefilter_query = f"SELECT DISTINCT dataset_id FROM {pg_table} {where_clause}"

    try:
        prefiltered_dataset_ids = [
            row["dataset_id"]
            for row in await pg_conn.fetch(prefilter_query, *pg_params)
        ]
    except asyncpg.exceptions.UndefinedColumnError as e:
        raise HTTPException(status_code=400, detail=f"Invalid filter field: {e}")

    if not prefiltered_dataset_ids:
        return {"total_datasets_count": 0, "facet_counts": {}, "datasets": []}

    # 2. Filter Datasets & Get Facets (Elastic)
    MODALITY_CASE_MAPPING = {
        "perturb-seq": "Perturb-seq",
        "crispr-screen": "CRISPR screen",
        "mave": "MAVE",
    }
    es_modality = MODALITY_CASE_MAPPING.get(modality, modality)

    es_query_body: Dict[str, Any] = {
        "query": {
            "bool": {
                "filter": [
                    {"term": {"data_modalities": es_modality}},
                    {"terms": {"dataset_id": prefiltered_dataset_ids}},
                ]
            }
        },
        "aggs": {
            field: {"terms": {"field": es_field, "size": 100}}
            for field, es_field in ELASTIC_FIELD_MAPPING.items()
            if field not in ("dataset_id", "dataset_score_interpretation")
        },
    }

    # Add dataset_* filters to ES query
    for key, value in query_params.items():
        if key in ELASTIC_FIELD_MAPPING:
            es_field = ELASTIC_FIELD_MAPPING[key]
            es_query_body["query"]["bool"]["filter"].append({"term": {es_field: value}})
        elif key == "dataset_metadata" and value:
            es_query_body["query"]["bool"]["must"] = [
                {"query_string": {"query": value}}
            ]

    es_result = await es_client.search(
        index="dataset-summary-v3",
        body=es_query_body,
        size=10000,  # Get all matching datasets to apply pagination later
    )

    total_datasets_count = es_result["hits"]["total"]["value"]
    es_datasets = [hit["_source"] for hit in es_result["hits"]["hits"]]

    facet_counts = {
        api_field: [
            {"value": bucket["key"], "count": bucket["doc_count"]}
            for bucket in es_result["aggregations"][api_field]["buckets"]
        ]
        for api_field in es_result.get("aggregations", {})
    }

    # 3. Paginate Datasets
    paginated_datasets = es_datasets[dataset_offset : dataset_offset + dataset_limit]

    # 4. Fetch Data (Postgres)
    final_datasets = []
    for es_dataset in paginated_datasets:
        dataset_id = es_dataset["dataset_id"]

        # Re-apply filters for this specific dataset
        current_pg_filters = [f"dataset_id = ${len(pg_params) + 1}"] + pg_filters
        current_pg_params = pg_params + [dataset_id]

        where_clause = f"WHERE {' AND '.join(current_pg_filters)}"

        order_by_clause = ""
        if sort:
            sort_clauses = []
            for sort_param in sort.split(","):
                field, __, direction = sort_param.partition(":")
                direction = "DESC" if direction == "desc" else "ASC"
                if field in api_to_db:
                    sort_clauses.append(f"{api_to_db[field]} {direction}")
            if sort_clauses:
                order_by_clause = f"ORDER BY {', '.join(sort_clauses)}"

        data_query = f"SELECT * FROM {pg_table} {where_clause} {order_by_clause} LIMIT {rows_per_dataset_limit}"

        pg_rows = await pg_conn.fetch(data_query, *current_pg_params)
        pg_rows_dict = [dict(row) for row in pg_rows]

        if modality == "perturb-seq":
            pg_rows_dict = await enrich_perturb_seq_rows(
                pg_conn, dataset_id, pg_rows_dict
            )

        # 5. Assemble Response
        results = []
        for row in pg_rows_dict:
            perturbation = {
                k.replace("perturbation_", ""): row.get(v)
                for k, v in api_to_db.items()
                if k.startswith("perturbation_")
            }
            effect = {
                k.replace("effect_", ""): row.get(v)
                for k, v in api_to_db.items()
                if k.startswith("effect_")
            }

            # Manual additions/transformations
            if modality == "perturb-seq":
                perturbation.update(
                    {
                        "n_total": row.get("perturbation_n_total"),
                        "n_up": row.get("perturbation_n_up"),
                        "n_down": row.get("perturbation_n_down"),
                    }
                )
                effect.update(
                    {
                        "direction": row.get("effect_direction"),
                        "n_total": row.get("effect_n_total"),
                        "n_up": row.get("effect_n_up"),
                        "n_down": row.get("effect_n_down"),
                    }
                )

            results.append({"perturbation": perturbation, "effect": effect})

        # Map ES fields to final dataset metadata
        def get_first_or_none(data: Optional[list]):
            return data[0] if data else None

        dataset_meta = {
            "dataset_id": es_dataset.get("dataset_id"),
            "dataset_tissue": get_first_or_none(es_dataset.get("tissue_labels")),
            "dataset_cell_type": get_first_or_none(es_dataset.get("cell_type_labels")),
            "dataset_cell_line": get_first_or_none(es_dataset.get("cell_line_labels")),
            "dataset_sex": get_first_or_none(es_dataset.get("sex_labels")),
            "dataset_developmental_stage": get_first_or_none(
                es_dataset.get("developmental_stage_labels")
            ),
            "dataset_disease": get_first_or_none(es_dataset.get("disease_labels")),
            "dataset_library_perturbation_type": get_first_or_none(
                es_dataset.get("library_perturbation_type_labels")
            ),
            "dataset_license_id": get_first_or_none(es_dataset.get("license_ids")),
            "dataset_license_label": get_first_or_none(
                es_dataset.get("license_labels")
            ),
            "dataset_score_interpretation": es_dataset.get("score_interpretation"),
        }

        final_datasets.append({"dataset": dataset_meta, "results": results})

    return {
        "total_datasets_count": total_datasets_count,
        "facet_counts": facet_counts,
        "datasets": final_datasets,
    }


async def _search_dataset_impl(
    modality: MODALITIES,
    dataset_id: str,
    query_params: Dict[str, Any],
):
    """Search within a specific dataset in a modality (Shared Implementation)."""
    limit = query_params.get("limit", 50)
    offset = query_params.get("offset", 0)
    sort = query_params.get("sort")

    validate_query_params(query_params, modality, dataset_id)

    pg_conn = db_pools["pg"]
    pg_table = PG_TABLES[modality]
    api_to_db = get_api_to_db_mapping(modality)

    # Build filters
    pg_filters = [f"dataset_id = ${1}"]
    pg_params: List[Any] = [dataset_id]

    for key, value in query_params.items():
        if key in api_to_db:
            db_field = api_to_db[key]
            if isinstance(value, str) and "_" not in value:
                pg_filters.append(f"{db_field} = ${len(pg_params) + 1}")
                pg_params.append(value)
            elif isinstance(value, str):
                condition, params = parse_numeric_filter(db_field, value)
                condition = condition.replace("$...", f"${len(pg_params) + 1}")
                if " AND " in condition:
                    condition = condition.replace("$...", f"${len(pg_params) + 2}")
                pg_filters.append(condition)
                pg_params.extend(params)

    where_clause = f"WHERE {' AND '.join(pg_filters)}"

    # 1. Count Rows
    count_query = f"SELECT COUNT(*) FROM {pg_table} {where_clause}"
    try:
        total_rows_count = await pg_conn.fetchval(count_query, *pg_params)
    except asyncpg.exceptions.UndefinedColumnError as e:
        raise HTTPException(status_code=400, detail=f"Invalid filter field: {e}")

    # 2. Fetch Rows
    order_by_clause = ""
    if sort:
        sort_clauses = []
        for sort_param in sort.split(","):
            field, __, direction = sort_param.partition(":")
            direction = "DESC" if direction == "desc" else "ASC"
            if field in api_to_db:
                sort_clauses.append(f"{api_to_db[field]} {direction}")
        if sort_clauses:
            order_by_clause = f"ORDER BY {', '.join(sort_clauses)}"

    data_query = f"SELECT * FROM {pg_table} {where_clause} {order_by_clause} LIMIT {limit} OFFSET {offset}"

    pg_rows = await pg_conn.fetch(data_query, *pg_params)
    pg_rows_dict = [dict(row) for row in pg_rows]

    if modality == "perturb-seq":
        pg_rows_dict = await enrich_perturb_seq_rows(pg_conn, dataset_id, pg_rows_dict)

    # 3. Assemble Response
    results = []
    for row in pg_rows_dict:
        perturbation = {
            k.replace("perturbation_", ""): row.get(v)
            for k, v in api_to_db.items()
            if k.startswith("perturbation_")
        }
        effect = {
            k.replace("effect_", ""): row.get(v)
            for k, v in api_to_db.items()
            if k.startswith("effect_")
        }

        if modality == "perturb-seq":
            perturbation.update(
                {
                    "n_total": row.get("perturbation_n_total"),
                    "n_up": row.get("perturbation_n_up"),
                    "n_down": row.get("perturbation_n_down"),
                }
            )
            effect.update(
                {
                    "direction": row.get("effect_direction"),
                    "n_total": row.get("effect_n_total"),
                    "n_up": row.get("effect_n_up"),
                    "n_down": row.get("effect_n_down"),
                }
            )
        results.append({"perturbation": perturbation, "effect": effect})

    return {
        "total_rows_count": total_rows_count,
        "offset": offset,
        "limit": limit,
        "results": results,
    }


# --- API Endpoints ---


@router.get(
    "/v1/mave/search",
    response_model=ModalitySearchResponse,
    response_model_by_alias=False,
)
async def search_mave(
    common: CommonModalitySearchParams = Depends(),
    modality_params: MaveParams = Depends(),
):
    """Search across all MAVE datasets."""
    params = {**common.dict(), **modality_params.dict()}
    return await _search_modality_impl("mave", params)


@router.get(
    "/v1/crispr-screen/search",
    response_model=ModalitySearchResponse,
    response_model_by_alias=False,
)
async def search_crispr_screen(
    common: CommonModalitySearchParams = Depends(),
    modality_params: CrisprScreenParams = Depends(),
):
    """Search across all CRISPR Screen datasets."""
    params = {**common.dict(), **modality_params.dict()}
    return await _search_modality_impl("crispr-screen", params)


@router.get(
    "/v1/perturb-seq/search",
    response_model=ModalitySearchResponse,
    response_model_by_alias=False,
)
async def search_perturb_seq(
    common: CommonModalitySearchParams = Depends(),
    modality_params: PerturbSeqParams = Depends(),
):
    """Search across all Perturb-seq datasets."""
    params = {**common.dict(), **modality_params.dict()}
    return await _search_modality_impl("perturb-seq", params)


@router.get(
    "/v1/mave/{dataset_id}/search",
    response_model=DatasetSearchResponse,
    response_model_by_alias=False,
)
async def search_mave_dataset(
    dataset_id: str,
    common: CommonDatasetSearchParams = Depends(),
    modality_params: MaveParams = Depends(),
):
    """Search within a specific MAVE dataset."""
    params = {**common.dict(), **modality_params.dict()}
    return await _search_dataset_impl("mave", dataset_id, params)


@router.get(
    "/v1/crispr-screen/{dataset_id}/search",
    response_model=DatasetSearchResponse,
    response_model_by_alias=False,
)
async def search_crispr_screen_dataset(
    dataset_id: str,
    common: CommonDatasetSearchParams = Depends(),
    modality_params: CrisprScreenParams = Depends(),
):
    """Search within a specific CRISPR Screen dataset."""
    params = {**common.dict(), **modality_params.dict()}
    return await _search_dataset_impl("crispr-screen", dataset_id, params)


@router.get(
    "/v1/perturb-seq/{dataset_id}/search",
    response_model=DatasetSearchResponse,
    response_model_by_alias=False,
)
async def search_perturb_seq_dataset(
    dataset_id: str,
    common: CommonDatasetSearchParams = Depends(),
    modality_params: PerturbSeqParams = Depends(),
):
    """Search within a specific Perturb-seq dataset."""
    params = {**common.dict(), **modality_params.dict()}
    return await _search_dataset_impl("perturb-seq", dataset_id, params)
