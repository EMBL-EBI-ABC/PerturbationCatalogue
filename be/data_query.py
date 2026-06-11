import asyncio
import csv
import io
import os
import json
from collections import defaultdict
from typing import Any, Dict, List, Literal, Optional, Tuple

import asyncpg
from fastapi import APIRouter, Depends, HTTPException, Query, Request
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field, create_model

try:
    from .es_indexes import ES_DATASET_SUMMARY  # type: ignore
except ImportError:  # pragma: no cover - fallback for running as a script
    from es_indexes import ES_DATASET_SUMMARY  # type: ignore

# --- Database Connection Management ---
db_pools: Dict[str, Any] = {}


router = APIRouter()

# --- Constants and Mappings ---
MODALITIES = Literal["perturb-seq", "crispr-screen", "mave"]

PG_TABLE_SUFFIX = os.getenv("PG_TABLE_SUFFIX", "")
if PG_TABLE_SUFFIX and not PG_TABLE_SUFFIX.replace("_", "").isalnum():
    raise ValueError(
        "PG_TABLE_SUFFIX may only contain letters, numbers, and underscores"
    )

PG_TABLES = {
    "perturb-seq": f"perturb_seq_dea{PG_TABLE_SUFFIX}",
    "crispr-screen": f"crispr_data{PG_TABLE_SUFFIX}",
    "mave": f"mave_data{PG_TABLE_SUFFIX}",
}

# Field mappings from API to database
PERTURB_SEQ_PG_MAPPING = {
    "perturbation_gene_name": "perturbed_target_symbol",
    "effect_gene_name": "gene",
    "effect_log2fc": "log2foldchange",
    "effect_padj": "padj",
    "effect_score_name": "score_name",
    "effect_score_value": "score_value",
    "effect_cell_type": "cell_type",
}
PERTURB_SEQ_GSEA_PG_MAPPING = {
    "perturbation_gene_name": "perturbed_target_symbol",
    "gsea_term": "term",
    "gsea_sidak": "sidak",
    "effect_term": "term",
    "effect_es": "es",
    "effect_nes": "nes",
    "effect_pval": "pval",
    "effect_sidak": "sidak",
    "effect_fdr": "fdr",
    "effect_geneset_size": "geneset_size",
    "effect_leading_edge": "leading_edge",
    "effect_cell_type": "cell_type",
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
    "perturbation_position": "perturbation_position",
    "perturbation_aa_wt": "perturbation_aa_wt",
    "perturbation_aa_change": "perturbation_aa_change",
    "effect_score_name": "score_name",
    "effect_score_value": "score_value",
}
PG_MAPPINGS = {
    "perturb-seq": PERTURB_SEQ_PG_MAPPING,
    "crispr-screen": CRISPR_PG_MAPPING,
    "mave": MAVE_PG_MAPPING,
}

# Numeric field mappings: "int" for integer fields, "float" for float fields
NUMERIC_FIELDS = {
    "perturb-seq": {
        "effect_log2fc": "float",
        "effect_padj": "float",
        "effect_score_value": "float",
        "gsea_sidak": "float",
    },
    "crispr-screen": {
        "effect_score_value": "float",
    },
    "mave": {
        "perturbation_position": "int",
        "effect_score_value": "float",
    },
}

# Default sorts for different modalities
DEFAULT_SORTS = {
    "crispr-screen": "effect_significant:desc",
    "perturb-seq": "effect_padj:asc",
}

# Load dataset metadata configuration
METADATA_PATH = os.path.join(os.path.dirname(__file__), "dataset_metadata.json")
with open(METADATA_PATH) as f:
    METADATA_CONFIG = json.load(f)

DATASET_FIELDS = METADATA_CONFIG["fields"]

ELASTIC_FIELD_MAPPING = {f["api_name"]: f["es_field"] for f in DATASET_FIELDS}
ELASTIC_AGG_FIELDS = {
    f["api_name"]: f["es_field"]
    for f in DATASET_FIELDS
    if f.get("es_type") != "text"
    and f["api_name"] not in ("dataset_id", "dataset_score_interpretation")
}


# --- Pydantic Models for API Response ---


# Dynamic Dataset Metadata Model
def _build_dataset_metadata_model():
    fields = {}
    for f in DATASET_FIELDS:
        name = f["api_name"].replace("dataset_", "")
        if name == "id":
            fields[name] = (str, Field(..., alias=f["api_name"]))
        elif f.get("es_type") == "boolean":
            fields[name] = (Optional[bool], Field(None, alias=f["api_name"]))
        elif f.get("is_array"):
            fields[name] = (Optional[List[str]], Field(None, alias=f["api_name"]))
        else:
            fields[name] = (Optional[str], Field(None, alias=f["api_name"]))
    return create_model("DatasetMetadata", **fields)


DatasetMetadata = _build_dataset_metadata_model()

# --- Pydantic Models for API Response ---


# Perturbation Models
class PerturbationBase(BaseModel):
    gene_name: str = Field(..., alias="perturbation_gene_name")


class MavePerturbation(PerturbationBase):
    name: Optional[str] = Field(None, alias="perturbation_name")
    position: Optional[int] = Field(None, alias="perturbation_position")
    aa_wt: Optional[str] = Field(None, alias="perturbation_aa_wt")
    aa_change: Optional[str] = Field(None, alias="perturbation_aa_change")


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
    score_name: Optional[str] = Field(None, alias="effect_score_name")
    score_value: Optional[float] = Field(None, alias="effect_score_value")
    cell_type: Optional[str] = Field(None, alias="effect_cell_type")
    n_total: Optional[int] = Field(None, alias="effect_n_total")
    n_up: Optional[int] = Field(None, alias="effect_n_up")
    n_down: Optional[int] = Field(None, alias="effect_n_down")


class PerturbSeqGseaEffect(EffectBase):
    term: str = Field(..., alias="effect_term")
    es: float = Field(..., alias="effect_es")
    nes: float = Field(..., alias="effect_nes")
    pval: float = Field(..., alias="effect_pval")
    sidak: float = Field(..., alias="effect_sidak")
    fdr: float = Field(..., alias="effect_fdr")
    geneset_size: int = Field(..., alias="effect_geneset_size")
    leading_edge: Optional[str] = Field(None, alias="effect_leading_edge")
    cell_type: Optional[str] = Field(None, alias="effect_cell_type")


class ScoreEffect(EffectBase):
    score_name: str = Field(..., alias="effect_score_name")
    score_value: float = Field(..., alias="effect_score_value")


# Result Models
class Result(BaseModel):
    perturbation: Dict
    effect: Dict


class GseaResult(BaseModel):
    perturbation: Dict
    effects: List[Dict]


# Dataset Models


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


# Dynamic Search Params
def _build_search_params_class():
    fields = {
        "dataset_metadata": (
            Optional[str],
            Query(None, description="Search term for dataset metadata"),
        )
    }
    for f in DATASET_FIELDS:
        fields[f["api_name"]] = (
            Optional[str],
            Query(None, description=f.get("description", f"Filter by {f['api_name']}")),
        )

    # Standard paging and sort
    fields.update(
        {
            "dataset_limit": (
                int,
                Query(10, description="Number of datasets to return"),
            ),
            "dataset_offset": (int, Query(0, description="Offset for datasets")),
            "rows_per_dataset_limit": (
                int,
                Query(20, description="Number of rows per dataset to return"),
            ),
            "sort": (
                Optional[str],
                Query(None, description="Sort order (e.g., 'field:asc,other:desc')"),
            ),
        }
    )

    # Create a Pydantic model and then wrap it or use it as-is in Depends
    Model = create_model("CommonModalitySearchParams", **fields)

    # FastAPI's Depends() works with Pydantic models. We'll use this model.
    return Model


CommonModalitySearchParams = _build_search_params_class()


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
        perturbation_position: Optional[str] = Query(
            None, description="Filter by perturbation position (supports ranges)"
        ),
        perturbation_aa_wt: Optional[str] = Query(
            None, description="Filter by source amino acid"
        ),
        perturbation_aa_change: Optional[str] = Query(
            None, description="Filter by target amino acid"
        ),
    ):
        self.perturbation_gene_name = perturbation_gene_name
        self.perturbation_name = perturbation_name
        self.effect_score_name = effect_score_name
        self.effect_score_value = effect_score_value
        self.perturbation_position = perturbation_position
        self.perturbation_aa_wt = perturbation_aa_wt
        self.perturbation_aa_change = perturbation_aa_change

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
        effect_score_name: Optional[str] = Query(
            None, description="Filter by effect score name"
        ),
        effect_score_value: Optional[str] = Query(
            None, description="Filter by effect score value (supports ranges)"
        ),
        effect_cell_type: Optional[str] = Query(
            None, description="Filter by cell type"
        ),
        gsea_term: Optional[str] = Query(None, description="Filter GSEA by term"),
        gsea_sidak: Optional[str] = Query(
            None, description="Filter GSEA by sidak (supports ranges)"
        ),
    ):
        self.perturbation_gene_name = perturbation_gene_name
        self.effect_gene_name = effect_gene_name
        self.effect_log2fc = effect_log2fc
        self.effect_padj = effect_padj
        self.effect_score_name = effect_score_name
        self.effect_score_value = effect_score_value
        self.effect_cell_type = effect_cell_type
        self.gsea_term = gsea_term
        self.gsea_sidak = gsea_sidak

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
        "sort",
        "dataset_limit",
        "dataset_offset",
        "rows_per_dataset_limit",
        "limit",
        "offset",
    }
    # Add all dynamic dataset params
    valid_params.update(ELASTIC_FIELD_MAPPING.keys())

    if dataset_id:
        valid_params = {"sort", "limit", "offset"}

    # Add all filterable perturbation and effect fields to valid_params
    pg_mapping = get_api_to_db_mapping(modality)
    valid_params.update(pg_mapping.keys())
    if modality == "perturb-seq":
        valid_params.update(PERTURB_SEQ_GSEA_PG_MAPPING.keys())

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
    # Fetch effect summaries keys
    effect_keys = list(set((row["dataset_id"], row["gene"]) for row in rows))
    effect_summary_map = {}

    # Fetch summaries in parallel
    pert_task = None
    if pert_keys:
        pert_dataset_ids = [k[0] for k in pert_keys]
        pert_symbols = [k[1] for k in pert_keys]
        pert_task = conn.fetch(
            """
            SELECT t.dataset_id, t.perturbed_target_symbol, t.n_total, t.n_up, t.n_down
            FROM perturb_seq_summary_perturbation AS t
            JOIN unnest($1::text[], $2::text[]) AS keys(did, pts)
            ON t.dataset_id = keys.did AND t.perturbed_target_symbol = keys.pts
            """,
            pert_dataset_ids,
            pert_symbols,
        )

    effect_task = None
    if effect_keys:

        effect_dataset_ids = [k[0] for k in effect_keys]
        effect_genes = [k[1] for k in effect_keys]
        effect_task = conn.fetch(
            """
            SELECT t.dataset_id, t.gene, t.n_total, t.n_up, t.n_down
            FROM perturb_seq_summary_effect AS t
            JOIN unnest($1::text[], $2::text[]) AS keys(did, g)
            ON t.dataset_id = keys.did AND t.gene = keys.g
            """,
            effect_dataset_ids,
            effect_genes,
        )

    # Wait for both tasks if they were created
    tasks = []
    if pert_task:
        tasks.append(pert_task)
    if effect_task:
        tasks.append(effect_task)

    if tasks:
        results = await asyncio.gather(*tasks)

        # Unpack results
        res_idx = 0
        if pert_task:
            pert_summary_rows = results[res_idx]
            res_idx += 1
            for r in pert_summary_rows:
                pert_summary_map[(r["dataset_id"], r["perturbed_target_symbol"])] = r

        if effect_task:
            effect_summary_rows = results[res_idx]
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

    return rows


async def _fetch_perturb_seq_gsea(
    conn: asyncpg.Connection,
    dataset_id: str,
    query_params: Dict[str, Any],
) -> List[Dict]:
    """Fetches GSEA data for a perturb-seq dataset."""
    pg_filters = ["dataset_id = $1"]
    pg_params = [dataset_id]

    # Re-use perturbation_gene_name filter if present
    if "perturbation_gene_name" in query_params:
        pg_filters.append(f"perturbed_target_symbol = ${len(pg_params) + 1}")
        pg_params.append(query_params["perturbation_gene_name"])

    # GSEA specific filters
    if "gsea_term" in query_params:
        pg_filters.append(f"term = ${len(pg_params) + 1}")
        pg_params.append(query_params["gsea_term"])

    if "gsea_sidak" in query_params:
        condition, params = parse_numeric_filter("sidak", query_params["gsea_sidak"])
        condition = condition.replace("$...", f"${len(pg_params) + 1}", 1)
        if " AND " in condition:
            condition = condition.replace("$...", f"${len(pg_params) + 2}", 1)
        pg_filters.append(condition)
        pg_params.extend(params)

    where_clause = f"WHERE {' AND '.join(pg_filters)}"
    query = f"""
        SELECT *
        FROM perturb_seq_gsea
        {where_clause}
        ORDER BY sidak ASC
        LIMIT 50
    """
    rows = await conn.fetch(query, *pg_params)
    return [dict(row) for row in rows]


# --- Shared Implementation Functions ---


async def _search_modality_impl(
    modality: MODALITIES,
    query_params: Dict[str, Any],
):
    """Search across all datasets within a modality (Shared Implementation)."""
    dataset_limit = query_params.get("dataset_limit", 10)
    dataset_offset = query_params.get("dataset_offset", 0)
    rows_per_dataset_limit = query_params.get("rows_per_dataset_limit", 20)
    sort = query_params.get("sort") or DEFAULT_SORTS.get(modality)

    # Check if position range filter is specified - if so, return all matching rows
    has_position_range = (
        modality == "mave"
        and "perturbation_position" in query_params
        and isinstance(query_params["perturbation_position"], str)
        and "_" in query_params["perturbation_position"]
    )

    validate_query_params(query_params, modality)

    pg_conn = db_pools["pg"]
    es_client = db_pools["es"]

    pg_table = PG_TABLES[modality]
    api_to_db = get_api_to_db_mapping(modality)

    # 1. Pre-filter Datasets (Postgres)
    pg_filters = []

    # Exclude "null" rows
    essential_columns = {
        "perturb-seq": "gene",
        "crispr-screen": "score_name",
        "mave": "score_name",
    }
    if modality in essential_columns:
        pg_filters.append(f"{essential_columns[modality]} IS NOT NULL")

    pg_params: List[Any] = []
    numeric_fields = NUMERIC_FIELDS.get(modality, {})

    for key, value in query_params.items():
        if key in api_to_db:
            db_field = api_to_db[key]
            field_type = numeric_fields.get(key)

            # Numeric field handling
            if field_type and isinstance(value, str):
                # Numeric range filter (contains "_")
                if "_" in value:
                    condition, params = parse_numeric_filter(db_field, value)
                    # This is a bit tricky because parse_numeric_filter doesn't know the param index
                    condition = condition.replace("$...", f"${len(pg_params) + 1}", 1)
                    if " AND " in condition:
                        condition = condition.replace(
                            "$...", f"${len(pg_params) + 2}", 1
                        )
                    pg_filters.append(condition)
                    pg_params.extend(params)
                # Simple numeric filter (no "_")
                else:
                    pg_filters.append(f"{db_field} = ${len(pg_params) + 1}")
                    if field_type == "int":
                        pg_params.append(int(value))
                    else:  # float
                        pg_params.append(float(value))
            # Simple string filter
            elif isinstance(value, str) and "_" not in value:
                pg_filters.append(f"{db_field} = ${len(pg_params) + 1}")
                pg_params.append(value)
            # Numeric range filter (for fields not explicitly in NUMERIC_FIELDS but using range syntax)
            elif isinstance(value, str):
                condition, params = parse_numeric_filter(db_field, value)
                # This is a bit tricky because parse_numeric_filter doesn't know the param index
                condition = condition.replace("$...", f"${len(pg_params) + 1}", 1)
                if " AND " in condition:
                    condition = condition.replace("$...", f"${len(pg_params) + 2}", 1)
                pg_filters.append(condition)
                pg_params.extend(params)

    where_clause = f"WHERE {' AND '.join(pg_filters)}" if pg_filters else ""

    if modality == "crispr-screen":
        prefilter_query = f"""
            SELECT dataset_id
            FROM {pg_table}
            {where_clause}
            GROUP BY dataset_id
            ORDER BY MAX(CASE WHEN significant = 'True' THEN 1 ELSE 0 END) DESC
        """
    elif modality == "perturb-seq":
        prefilter_query = f"""
            SELECT dataset_id
            FROM {pg_table}
            {where_clause}
            GROUP BY dataset_id
            ORDER BY COUNT(*) FILTER (WHERE padj < 0.05) DESC
        """
    else:
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

    es_query_body: Dict[str, Any] = {
        "query": {
            "bool": {
                "filter": [
                    {
                        "term": {
                            "data_modalities": MODALITY_CASE_MAPPING.get(
                                modality, modality
                            )
                        }
                    },
                    {"terms": {"dataset_id": prefiltered_dataset_ids}},
                ]
            }
        },
        "aggs": {
            field: {"terms": {"field": es_field, "size": 100}}
            for field, es_field in ELASTIC_AGG_FIELDS.items()
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
        index=ES_DATASET_SUMMARY,
        body=es_query_body,
        size=10000,  # Get all matching datasets to apply pagination later
    )

    total_datasets_count = es_result["hits"]["total"]["value"]
    es_datasets = [hit["_source"] for hit in es_result["hits"]["hits"]]

    # Re-order es_datasets to match prefiltered_dataset_ids order (Postgres significance sort)
    dataset_id_to_order = {did: i for i, did in enumerate(prefiltered_dataset_ids)}
    es_datasets.sort(key=lambda x: dataset_id_to_order.get(x["dataset_id"], 999999))

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

        # Only apply LIMIT if not using position range filter (which should return all matching rows)
        limit_clause = (
            f"LIMIT {rows_per_dataset_limit}" if not has_position_range else ""
        )
        data_query = f"SELECT * FROM {pg_table} {where_clause} {order_by_clause} {limit_clause}".strip()

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
                log2fc = row.get("log2foldchange")
                if log2fc is None:
                    effect["direction"] = "not available"
                elif log2fc > 0:
                    effect["direction"] = "increased"
                elif log2fc < 0:
                    effect["direction"] = "decreased"
                else:
                    effect["direction"] = "no change"
                effect.update(
                    {
                        "n_total": row.get("effect_n_total"),
                        "n_up": row.get("effect_n_up"),
                        "n_down": row.get("effect_n_down"),
                    }
                )

            results.append({"perturbation": perturbation, "effect": effect})

        # Map ES fields to final dataset metadata
        dataset_meta = {}
        for f in DATASET_FIELDS:
            val = es_dataset.get(f["es_field"])
            dataset_meta[f["api_name"]] = val

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
    sort = query_params.get("sort") or DEFAULT_SORTS.get(modality)

    # Check if position range filter is specified - if so, return all matching rows
    has_position_range = (
        modality == "mave"
        and "perturbation_position" in query_params
        and isinstance(query_params["perturbation_position"], str)
        and "_" in query_params["perturbation_position"]
    )

    validate_query_params(query_params, modality, dataset_id)

    pg_conn = db_pools["pg"]
    pg_table = PG_TABLES[modality]
    api_to_db = get_api_to_db_mapping(modality)

    # Build filters
    pg_filters = [f"dataset_id = ${1}"]

    # Exclude "null" rows
    essential_columns = {
        "perturb-seq": "gene",
        "crispr-screen": "score_name",
        "mave": "score_name",
    }
    if modality in essential_columns:
        pg_filters.append(f"{essential_columns[modality]} IS NOT NULL")

    pg_params: List[Any] = [dataset_id]
    numeric_fields = NUMERIC_FIELDS.get(modality, {})

    for key, value in query_params.items():
        if key in api_to_db:
            db_field = api_to_db[key]
            field_type = numeric_fields.get(key)

            # Numeric field handling
            if field_type and isinstance(value, str):
                # Numeric range filter (contains "_")
                if "_" in value:
                    condition, params = parse_numeric_filter(db_field, value)
                    condition = condition.replace("$...", f"${len(pg_params) + 1}", 1)
                    if " AND " in condition:
                        condition = condition.replace(
                            "$...", f"${len(pg_params) + 2}", 1
                        )
                    pg_filters.append(condition)
                    pg_params.extend(params)
                # Simple numeric filter (no "_")
                else:
                    pg_filters.append(f"{db_field} = ${len(pg_params) + 1}")
                    if field_type == "int":
                        pg_params.append(int(value))
                    else:  # float
                        pg_params.append(float(value))
            # Simple string filter
            elif isinstance(value, str) and "_" not in value:
                pg_filters.append(f"{db_field} = ${len(pg_params) + 1}")
                pg_params.append(value)
            # Numeric range filter (for fields not explicitly in NUMERIC_FIELDS but using range syntax)
            elif isinstance(value, str):
                condition, params = parse_numeric_filter(db_field, value)
                condition = condition.replace("$...", f"${len(pg_params) + 1}", 1)
                if " AND " in condition:
                    condition = condition.replace("$...", f"${len(pg_params) + 2}", 1)
                pg_filters.append(condition)
                pg_params.extend(params)

    where_clause = f"WHERE {' AND '.join(pg_filters)}"

    # 1. Count Rows
    no_user_filters = not any(k in api_to_db for k in query_params)
    if modality == "perturb-seq" and no_user_filters:
        count_query = (
            "SELECT n_total FROM perturb_seq_summary_dataset WHERE dataset_id = $1"
        )
        count_params = [dataset_id]
    else:
        count_query = f"SELECT COUNT(*) FROM {pg_table} {where_clause}"
        count_params = pg_params

    try:
        total_rows_count = await pg_conn.fetchval(count_query, *count_params) or 0
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

    # Only apply LIMIT/OFFSET if not using position range filter (which should return all matching rows)
    if has_position_range:
        pagination_clause = ""
    else:
        pagination_clause = f"LIMIT {limit} OFFSET {offset}"
    data_query = f"SELECT * FROM {pg_table} {where_clause} {order_by_clause} {pagination_clause}".strip()

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
            log2fc = row.get("log2foldchange")
            if log2fc is None:
                effect["direction"] = "not available"
            elif log2fc > 0:
                effect["direction"] = "increased"
            elif log2fc < 0:
                effect["direction"] = "decreased"
            else:
                effect["direction"] = "no change"
            perturbation.update(
                {
                    "n_total": row.get("perturbation_n_total"),
                    "n_up": row.get("perturbation_n_up"),
                    "n_down": row.get("perturbation_n_down"),
                }
            )
            effect.update(
                {
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
    response_model_by_alias=True,
)
async def search_mave(
    common: CommonModalitySearchParams = Depends(),
    modality_params: MaveParams = Depends(),
):
    """Search across all MAVE datasets."""
    params = {**common.model_dump(exclude_none=True), **modality_params.dict()}
    return await _search_modality_impl("mave", params)


@router.get(
    "/v1/crispr-screen/search",
    response_model=ModalitySearchResponse,
    response_model_by_alias=True,
)
async def search_crispr_screen(
    common: CommonModalitySearchParams = Depends(),
    modality_params: CrisprScreenParams = Depends(),
):
    """Search across all CRISPR Screen datasets."""
    params = {**common.model_dump(exclude_none=True), **modality_params.dict()}
    return await _search_modality_impl("crispr-screen", params)


@router.get(
    "/v1/perturb-seq/search",
    response_model=ModalitySearchResponse,
    response_model_by_alias=True,
)
async def search_perturb_seq(
    common: CommonModalitySearchParams = Depends(),
    modality_params: PerturbSeqParams = Depends(),
):
    """Search across all Perturb-seq datasets."""
    params = {**common.model_dump(exclude_none=True), **modality_params.dict()}
    return await _search_modality_impl("perturb-seq", params)


@router.get(
    "/v1/mave/{dataset_id}/search",
    response_model=DatasetSearchResponse,
    response_model_by_alias=True,
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
    response_model_by_alias=True,
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
    response_model_by_alias=True,
)
async def search_perturb_seq_dataset(
    dataset_id: str,
    common: CommonDatasetSearchParams = Depends(),
    modality_params: PerturbSeqParams = Depends(),
):
    """Search within a specific Perturb-seq dataset."""
    params = {**common.dict(), **modality_params.dict()}
    return await _search_dataset_impl("perturb-seq", dataset_id, params)


@router.get(
    "/v1/perturb-seq-gsea",
    response_model=List[GseaResult],
    response_model_by_alias=False,
)
async def get_perturb_seq_gsea(
    dataset_id: str = Query(..., description="Mandatory dataset ID"),
    perturbed_gene_name: str = Query(
        ..., description="Mandatory perturbed gene symbol"
    ),
):
    """Retrieve GSEA results for a specific gene in a dataset."""
    pg_pool = db_pools.get("pg")
    if not pg_pool:
        raise HTTPException(status_code=500, detail="Database pool not initialized")
    async with pg_pool.acquire() as conn:
        # Fetch rows for this gene (no default filtering, return all)
        rows = await conn.fetch(
            """
            SELECT * FROM perturb_seq_gsea 
            WHERE dataset_id = $1 AND perturbed_target_symbol = $2
            ORDER BY sidak ASC
            """,
            dataset_id,
            perturbed_gene_name,
        )
        if not rows:
            return []

        # Group by perturbation
        gsea_by_pert = defaultdict(list)
        for r in rows:
            effect = {
                k.replace("effect_", ""): r.get(v)
                for k, v in PERTURB_SEQ_GSEA_PG_MAPPING.items()
                if k.startswith("effect_")
            }
            gsea_by_pert[r["perturbed_target_symbol"]].append(effect)

        # Enrich perturbation
        pert_summary = await conn.fetchrow(
            """
            SELECT n_total, n_up, n_down
            FROM perturb_seq_summary_perturbation
            WHERE dataset_id = $1 AND perturbed_target_symbol = $2
            """,
            dataset_id,
            perturbed_gene_name,
        )
        pert_summary = dict(pert_summary) if pert_summary else {}

        results = []
        for pert_symbol, effects in gsea_by_pert.items():
            results.append(
                {
                    "perturbation": {
                        "gene_name": pert_symbol,
                        "n_total": pert_summary.get("n_total"),
                        "n_up": pert_summary.get("n_up"),
                        "n_down": pert_summary.get("n_down"),
                    },
                    "effects": effects,
                }
            )
        return results


# CSV column definitions for each modality
CSV_COLUMNS = {
    "perturb-seq": [
        ("perturbation_gene_name", "Perturbation Gene"),
        ("effect_gene_name", "Effect Gene"),
        ("effect_log2fc", "Log2FC"),
        ("effect_padj", "Padj"),
        ("effect_score_name", "Score Name"),
        ("effect_score_value", "Score Value"),
        ("effect_cell_type", "Cell Type"),
    ],
    "crispr-screen": [
        ("perturbation_gene_name", "Perturbation Gene"),
        ("effect_score_name", "Score Name"),
        ("effect_score_value", "Score Value"),
        ("effect_significant", "Significant"),
        ("effect_significance_criteria", "Significance Criteria"),
    ],
    "mave": [
        ("perturbation_gene_name", "Perturbation Gene"),
        ("perturbation_name", "Perturbation Name"),
        ("perturbation_position", "Position"),
        ("perturbation_aa_wt", "AA WT"),
        ("perturbation_aa_change", "AA Change"),
        ("effect_score_name", "Score Name"),
        ("effect_score_value", "Score Value"),
    ],
}


def _results_to_csv(results: List[Dict], modality: MODALITIES) -> str:
    """Convert results to CSV format."""
    output = io.StringIO()
    columns = CSV_COLUMNS.get(modality, [])

    writer = csv.writer(output)
    # Write header
    writer.writerow([col[1] for col in columns])

    # Write data rows
    for result in results:
        perturbation = result.get("perturbation", {})
        effect = result.get("effect", {})

        row = []
        for field_key, _ in columns:
            if field_key.startswith("perturbation_"):
                key = field_key.replace("perturbation_", "")
                value = perturbation.get(key, "")
            elif field_key.startswith("effect_"):
                key = field_key.replace("effect_", "")
                value = effect.get(key, "")
            else:
                value = ""
            row.append(value if value is not None else "")
        writer.writerow(row)

    return output.getvalue()


@router.get("/v1/{modality}/download")
async def download_modality_data(
    modality: MODALITIES,
    common: CommonModalitySearchParams = Depends(),
    perturbation_gene_name: Optional[str] = Query(None),
    effect_gene_name: Optional[str] = Query(None),
    effect_log2fc: Optional[str] = Query(None),
    effect_padj: Optional[str] = Query(None),
    effect_score_name: Optional[str] = Query(None),
    effect_score_value: Optional[str] = Query(None),
    effect_cell_type: Optional[str] = Query(None),
    effect_significant: Optional[str] = Query(None),
    effect_significance_criteria: Optional[str] = Query(None),
    perturbation_name: Optional[str] = Query(None),
    perturbation_position: Optional[str] = Query(None),
    perturbation_aa_wt: Optional[str] = Query(None),
    perturbation_aa_change: Optional[str] = Query(None),
):
    """Download data for a modality as CSV."""
    # Build params dict from all provided parameters
    params = common.model_dump(exclude_none=True)

    # Add modality-specific params
    modality_params = {
        "perturbation_gene_name": perturbation_gene_name,
        "effect_gene_name": effect_gene_name,
        "effect_log2fc": effect_log2fc,
        "effect_padj": effect_padj,
        "effect_score_name": effect_score_name,
        "effect_score_value": effect_score_value,
        "effect_cell_type": effect_cell_type,
        "effect_significant": effect_significant,
        "effect_significance_criteria": effect_significance_criteria,
        "perturbation_name": perturbation_name,
        "perturbation_position": perturbation_position,
        "perturbation_aa_wt": perturbation_aa_wt,
        "perturbation_aa_change": perturbation_aa_change,
    }
    params.update({k: v for k, v in modality_params.items() if v is not None})

    # Override limits for download - get more data
    params["dataset_limit"] = 1000
    params["rows_per_dataset_limit"] = 10000

    result = await _search_modality_impl(modality, params)

    # Flatten all results from all datasets
    all_results = []
    for dataset in result.get("datasets", []):
        all_results.extend(dataset.get("results", []))

    csv_content = _results_to_csv(all_results, modality)

    # Generate filename
    gene_name = perturbation_gene_name or effect_gene_name or "all"
    filename = f"{modality}_{gene_name}_data.csv"

    return StreamingResponse(
        iter([csv_content]),
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )


# GSEA CSV columns (excluding leading_edge)
GSEA_CSV_COLUMNS = [
    ("term", "Term"),
    ("es", "ES"),
    ("nes", "NES"),
    ("pval", "P-value"),
    ("sidak", "Sidak"),
    ("fdr", "FDR"),
    ("geneset_size", "Geneset Size"),
    ("cell_type", "Cell Type"),
]


def _gsea_results_to_csv(gsea_results: List[Dict]) -> str:
    """Convert GSEA results to CSV format."""
    output = io.StringIO()
    writer = csv.writer(output)

    # Write header
    writer.writerow([col[1] for col in GSEA_CSV_COLUMNS])

    # Flatten and write data rows
    for result in gsea_results:
        effects = result.get("effects") or []
        for effect in effects:
            row = []
            for field_key, _ in GSEA_CSV_COLUMNS:
                value = effect.get(field_key, "")
                row.append(value if value is not None else "")
            writer.writerow(row)

    return output.getvalue()


@router.get("/v1/{modality}/{dataset_id}/download")
async def download_dataset_data(
    modality: MODALITIES,
    dataset_id: str,
    # Common params
    limit: int = Query(100000, description="Maximum rows to download"),
    offset: int = Query(0, description="Offset for rows"),
    sort: Optional[str] = Query(None, description="Sort order"),
    # Perturb-seq params
    perturbation_gene_name: Optional[str] = Query(None),
    effect_gene_name: Optional[str] = Query(None),
    effect_log2fc: Optional[str] = Query(None),
    effect_padj: Optional[str] = Query(None),
    effect_score_name: Optional[str] = Query(None),
    effect_score_value: Optional[str] = Query(None),
    effect_cell_type: Optional[str] = Query(None),
    # CRISPR params
    effect_significant: Optional[str] = Query(None),
    effect_significance_criteria: Optional[str] = Query(None),
    # MAVE params
    perturbation_name: Optional[str] = Query(None),
    perturbation_position: Optional[str] = Query(None),
    perturbation_aa_wt: Optional[str] = Query(None),
    perturbation_aa_change: Optional[str] = Query(None),
):
    """Download data for a specific dataset as CSV."""
    # Build params dict
    params = {
        "limit": limit,
        "offset": offset,
    }
    if sort:
        params["sort"] = sort

    # Add modality-specific params
    modality_params = {
        "perturbation_gene_name": perturbation_gene_name,
        "effect_gene_name": effect_gene_name,
        "effect_log2fc": effect_log2fc,
        "effect_padj": effect_padj,
        "effect_score_name": effect_score_name,
        "effect_score_value": effect_score_value,
        "effect_cell_type": effect_cell_type,
        "effect_significant": effect_significant,
        "effect_significance_criteria": effect_significance_criteria,
        "perturbation_name": perturbation_name,
        "perturbation_position": perturbation_position,
        "perturbation_aa_wt": perturbation_aa_wt,
        "perturbation_aa_change": perturbation_aa_change,
    }
    params.update({k: v for k, v in modality_params.items() if v is not None})

    result = await _search_dataset_impl(modality, dataset_id, params)

    csv_content = _results_to_csv(result.get("results", []), modality)

    filename = f"{modality}_{dataset_id}_data.csv"

    return StreamingResponse(
        iter([csv_content]),
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )


@router.get("/v1/perturb-seq-gsea/download")
async def download_perturb_seq_gsea(
    dataset_id: str = Query(..., description="Mandatory dataset ID"),
    perturbed_gene_name: str = Query(
        ..., description="Mandatory perturbed gene symbol"
    ),
):
    """Download GSEA data for a specific gene in a dataset as CSV."""
    # Reuse the existing GSEA endpoint logic
    pg_pool = db_pools.get("pg")
    if not pg_pool:
        raise HTTPException(status_code=500, detail="Database pool not initialized")

    async with pg_pool.acquire() as conn:
        rows = await conn.fetch(
            """
            SELECT * FROM perturb_seq_gsea
            WHERE dataset_id = $1 AND perturbed_target_symbol = $2
            ORDER BY sidak ASC
            """,
            dataset_id,
            perturbed_gene_name,
        )

        if not rows:
            csv_content = _gsea_results_to_csv([])
        else:
            # Build results in the same format as the main GSEA endpoint
            gsea_by_pert = defaultdict(list)
            for r in rows:
                effect = {
                    k.replace("effect_", ""): r.get(v)
                    for k, v in PERTURB_SEQ_GSEA_PG_MAPPING.items()
                    if k.startswith("effect_")
                }
                gsea_by_pert[r["perturbed_target_symbol"]].append(effect)

            results = [
                {"perturbation": {"gene_name": pert_symbol}, "effects": effects}
                for pert_symbol, effects in gsea_by_pert.items()
            ]
            csv_content = _gsea_results_to_csv(results)

    filename = f"gsea_{perturbed_gene_name}_{dataset_id}_data.csv"

    return StreamingResponse(
        iter([csv_content]),
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )
