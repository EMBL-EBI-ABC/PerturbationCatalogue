import json
import logging
import os
import uuid
from typing import Any, Dict, Optional

from fastapi import APIRouter
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

from google import genai
from google.genai import types

from data_query import db_pools

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/v1/chat", tags=["AI Chat"])

# --- Configuration ---

_config: Dict[str, Any] = {}

# In-memory session store: session_id -> list of content dicts
_sessions: Dict[str, list] = {}


def configure(
    google_cloud_project: str = "",
    gemini_model: str = "gemini-2.5-flash",
    gemini_api_key: str = "",
):
    _config["google_cloud_project"] = google_cloud_project
    _config["gemini_model"] = gemini_model
    _config["gemini_api_key"] = gemini_api_key


# --- Request model ---


class ChatRequest(BaseModel):
    message: str
    session_id: Optional[str] = None


# --- Tool definitions for Gemini function calling ---

SEARCH_DATASETS_DECLARATION = types.FunctionDeclaration(
    name="search_datasets",
    description=(
        "Search for datasets in the Perturbation Catalogue by metadata. "
        "Returns dataset metadata including titles, authors, modalities, tissues, cell lines, diseases."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "query": types.Schema(
                type="STRING",
                description="Free-text search query for dataset metadata (titles, authors, summaries)",
            ),
            "modality": types.Schema(
                type="STRING",
                description="Data modality filter: 'Perturb-seq', 'CRISPR screen', or 'MAVE'",
            ),
            "tissue": types.Schema(
                type="STRING",
                description="Tissue filter (e.g., 'brain', 'lung', 'blood')",
            ),
            "cell_line": types.Schema(
                type="STRING",
                description="Cell line filter (e.g., 'K562', 'A549')",
            ),
            "disease": types.Schema(
                type="STRING",
                description="Disease filter (e.g., 'cancer', 'leukemia')",
            ),
            "max_results": types.Schema(
                type="INTEGER",
                description="Maximum number of datasets to return (default 10)",
            ),
        },
    ),
)

SEARCH_TARGET_SUMMARY_DECLARATION = types.FunctionDeclaration(
    name="search_target_summary",
    description=(
        "Look up gene target information aggregated across datasets. "
        "Shows what modalities, tissues, cell types, and diseases have been studied for a gene."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "gene_name": types.Schema(
                type="STRING",
                description="Gene symbol to look up (e.g., 'BRCA2', 'TP53', 'KRAS')",
            ),
            "query": types.Schema(
                type="STRING",
                description="Free-text search if gene_name is not known",
            ),
            "max_results": types.Schema(
                type="INTEGER",
                description="Maximum number of targets to return (default 5)",
            ),
        },
    ),
)

QUERY_PERTURBATION_DATA_DECLARATION = types.FunctionDeclaration(
    name="query_perturbation_data",
    description=(
        "Query row-level perturbation data from PostgreSQL. "
        "Returns actual experimental results: CRISPR scores, Perturb-seq differential expression, MAVE variant effects."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "modality": types.Schema(
                type="STRING",
                description="Data modality: 'perturb-seq', 'crispr-screen', or 'mave'",
            ),
            "gene_name": types.Schema(
                type="STRING",
                description="Perturbation gene name to filter by",
            ),
            "dataset_id": types.Schema(
                type="STRING",
                description="Specific dataset ID to query within",
            ),
            "limit": types.Schema(
                type="INTEGER",
                description="Maximum rows to return (default 20, max 50)",
            ),
        },
        required=["modality"],
    ),
)

GET_CATALOGUE_SUMMARY_DECLARATION = types.FunctionDeclaration(
    name="get_catalogue_summary",
    description=(
        "Get high-level summary statistics about the entire Perturbation Catalogue: "
        "total datasets, targets, modality counts, etc."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={},
    ),
)

CREATE_VISUALIZATION_DECLARATION = types.FunctionDeclaration(
    name="create_visualization",
    description=(
        "Create a visualization (table, pie chart, or bar chart) to display in the data portal below the chat. "
        "Always create visualizations when you have data to show the user. "
        "Use pie charts for 2-6 categories, bar charts for comparisons or >6 categories, tables for detailed data."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "viz_type": types.Schema(
                type="STRING",
                description="Visualization type: 'table', 'pie_chart', or 'bar_chart'",
            ),
            "title": types.Schema(
                type="STRING",
                description="Title for the visualization",
            ),
            "data": types.Schema(
                type="OBJECT",
                description=(
                    "Data for the visualization. "
                    "For table: {headers: [...], rows: [[...], ...]}. "
                    "For pie_chart: {labels: [...], values: [...]}. "
                    "For bar_chart: {labels: [...], values: [...], xlabel: '...', ylabel: '...'}."
                ),
                properties={
                    "headers": types.Schema(
                        type="ARRAY",
                        items=types.Schema(type="STRING"),
                        description="Column headers for table",
                    ),
                    "rows": types.Schema(
                        type="ARRAY",
                        items=types.Schema(
                            type="ARRAY",
                            items=types.Schema(type="STRING"),
                        ),
                        description="Row data for table",
                    ),
                    "labels": types.Schema(
                        type="ARRAY",
                        items=types.Schema(type="STRING"),
                        description="Labels for chart",
                    ),
                    "values": types.Schema(
                        type="ARRAY",
                        items=types.Schema(type="NUMBER"),
                        description="Values for chart",
                    ),
                    "xlabel": types.Schema(
                        type="STRING", description="X-axis label for bar chart"
                    ),
                    "ylabel": types.Schema(
                        type="STRING", description="Y-axis label for bar chart"
                    ),
                },
            ),
        },
        required=["viz_type", "title", "data"],
    ),
)

ALL_TOOL_DECLARATIONS = [
    SEARCH_DATASETS_DECLARATION,
    SEARCH_TARGET_SUMMARY_DECLARATION,
    QUERY_PERTURBATION_DATA_DECLARATION,
    GET_CATALOGUE_SUMMARY_DECLARATION,
    CREATE_VISUALIZATION_DECLARATION,
]

# --- Tool implementations ---

ES_DATASET_SUMMARY = "dataset-summary"
ES_TARGET_SUMMARY = "target-summary"
ES_LANDING_PAGE_SUMMARY = "landing-page-summary"

PG_TABLES = {
    "perturb-seq": "perturb_seq_dea",
    "crispr-screen": "crispr_data",
    "mave": "mave_data",
}


async def _tool_search_datasets(args: dict) -> dict:
    es = db_pools["es"]
    query_text = args.get("query", "")
    modality = args.get("modality")
    tissue = args.get("tissue")
    cell_line = args.get("cell_line")
    disease = args.get("disease")
    max_results = min(args.get("max_results", 10), 20)

    filters = []
    if modality:
        filters.append({"term": {"data_modalities": modality}})
    if tissue:
        filters.append(
            {
                "wildcard": {
                    "tissue_labels": {
                        "value": f"*{tissue.lower()}*",
                        "case_insensitive": True,
                    }
                }
            }
        )
    if cell_line:
        filters.append(
            {
                "wildcard": {
                    "cell_line_labels": {
                        "value": f"*{cell_line.lower()}*",
                        "case_insensitive": True,
                    }
                }
            }
        )
    if disease:
        filters.append(
            {
                "wildcard": {
                    "disease_labels": {
                        "value": f"*{disease.lower()}*",
                        "case_insensitive": True,
                    }
                }
            }
        )

    bool_query: Dict[str, Any] = {}
    if filters:
        bool_query["filter"] = filters
    if query_text:
        bool_query["must"] = [
            {
                "multi_match": {
                    "query": query_text,
                    "fields": [
                        "dataset_id",
                        "study_title^2",
                        "experiment_title^2",
                        "experiment_summary",
                        "first_author",
                        "last_author",
                        "tissue_labels",
                        "cell_type_labels",
                        "cell_line_labels",
                        "disease_labels",
                    ],
                    "type": "best_fields",
                    "fuzziness": "AUTO",
                }
            }
        ]

    es_query = {"bool": bool_query} if bool_query else {"match_all": {}}

    result = await es.search(
        index=ES_DATASET_SUMMARY,
        query=es_query,
        size=max_results,
        source=[
            "dataset_id",
            "study_title",
            "experiment_title",
            "experiment_summary",
            "data_modalities",
            "first_author",
            "last_author",
            "tissue_labels",
            "cell_type_labels",
            "cell_line_labels",
            "disease_labels",
            "license_labels",
        ],
    )

    hits = result.get("hits", {})
    total = hits.get("total", {}).get("value", 0)
    datasets = [hit["_source"] for hit in hits.get("hits", [])]

    return {"total": total, "datasets": datasets}


async def _tool_search_target_summary(args: dict) -> dict:
    es = db_pools["es"]
    gene_name = args.get("gene_name", "")
    query_text = args.get("query", "")
    max_results = min(args.get("max_results", 5), 20)

    if gene_name:
        es_query = {
            "bool": {
                "should": [
                    {"term": {"perturbed_target_symbol": gene_name.upper()}},
                    {
                        "wildcard": {
                            "perturbed_target_symbol": {
                                "value": f"*{gene_name}*",
                                "case_insensitive": True,
                            }
                        }
                    },
                ],
                "minimum_should_match": 1,
            }
        }
    elif query_text:
        es_query = {
            "multi_match": {
                "query": query_text,
                "fields": [
                    "perturbed_target_symbol^2",
                    "tissues_tested",
                    "cell_types_tested",
                    "diseases_tested",
                    "data_modalities",
                ],
                "fuzziness": "AUTO",
            }
        }
    else:
        es_query = {"match_all": {}}

    result = await es.search(
        index=ES_TARGET_SUMMARY,
        query=es_query,
        size=max_results,
    )

    hits = result.get("hits", {})
    total = hits.get("total", {}).get("value", 0)
    targets = [hit["_source"] for hit in hits.get("hits", [])]

    return {"total": total, "targets": targets}


async def _tool_query_perturbation_data(args: dict) -> dict:
    modality = args.get("modality", "crispr-screen")
    gene_name = args.get("gene_name")
    dataset_id = args.get("dataset_id")
    limit = min(args.get("limit", 20), 50)

    pg_table = PG_TABLES.get(modality)
    if not pg_table:
        return {"error": f"Unknown modality: {modality}"}

    pg_pool = db_pools["pg"]

    pg_filters = []
    pg_params = []

    # Exclude null essential columns
    essential_columns = {
        "perturb-seq": "gene",
        "crispr-screen": "score_name",
        "mave": "score_name",
    }
    if modality in essential_columns:
        pg_filters.append(f"{essential_columns[modality]} IS NOT NULL")

    if gene_name:
        pg_filters.append(f"perturbed_target_symbol = ${len(pg_params) + 1}")
        pg_params.append(gene_name.upper())

    if dataset_id:
        pg_filters.append(f"dataset_id = ${len(pg_params) + 1}")
        pg_params.append(dataset_id)

    where_clause = f"WHERE {' AND '.join(pg_filters)}" if pg_filters else ""

    query = f"SELECT * FROM {pg_table} {where_clause} LIMIT {limit}"
    rows = await pg_pool.fetch(query, *pg_params)
    rows_dict = [dict(r) for r in rows]

    # Convert to serializable format
    clean_rows = []
    for row in rows_dict:
        clean = {}
        for k, v in row.items():
            if v is None:
                clean[k] = None
            elif isinstance(v, (int, float, str, bool)):
                clean[k] = v
            else:
                clean[k] = str(v)
        clean_rows.append(clean)

    return {"total_returned": len(clean_rows), "rows": clean_rows}


async def _tool_get_catalogue_summary() -> dict:
    es = db_pools["es"]
    try:
        response = await es.get(index=ES_LANDING_PAGE_SUMMARY, id="summary")
        source = response.get("_source", {})
        return source
    except Exception as exc:
        return {"error": str(exc)}


TOOL_HANDLERS = {
    "search_datasets": _tool_search_datasets,
    "search_target_summary": _tool_search_target_summary,
    "query_perturbation_data": _tool_query_perturbation_data,
    "get_catalogue_summary": lambda args: _tool_get_catalogue_summary(),
}


# --- System instruction ---

SYSTEM_INSTRUCTION = """You are the AI Explorer for the Perturbation Catalogue, a comprehensive resource for human gene perturbation experiments (CRISPR screens, Perturb-seq, MAVE).

Your role:
- Help researchers explore perturbation data through natural language
- Search datasets, look up gene targets, and query experimental results
- Be concise and scientific, lead with key findings

CRITICAL RULES FOR RESPONSES:
- NEVER include markdown tables, data listings, or raw data rows in your text responses.
- Your text response should ONLY contain a brief summary and interpretation of the findings.
- ALL data (tables, charts) MUST be sent via the create_visualization tool, which displays them in a separate data portal below the chat.
- After calling create_visualization, write a short summary (2-3 sentences) interpreting the results. Do NOT repeat the data in text.

Visualization guidelines:
- Use pie charts for distributions with 2-6 categories
- Use bar charts for comparisons or >6 categories
- Use tables for detailed data rows (limit to 20 rows for readability)
- When showing datasets, include dataset_id, title, modality, and key metadata
- When showing perturbation data, highlight significant results

Available data modalities:
- Perturb-seq: Single-cell transcriptomic readout of gene perturbations. Key fields: perturbation gene, effect gene, log2FC, padj
- CRISPR screen: Fitness/viability screens. Key fields: perturbation gene, score name, score value, significant
- MAVE: Multiplexed Assay of Variant Effect. Key fields: perturbation gene, variant, position, score

When users ask about a gene, first search for target summary to understand what's available, then query specific data if needed."""

# --- SSE streaming endpoint ---


def _sse_event(event: str, data: dict) -> str:
    return f"event: {event}\ndata: {json.dumps(data)}\n\n"


@router.post("/stream")
async def chat_stream(request: ChatRequest):
    session_id = request.session_id or str(uuid.uuid4())
    message = request.message

    async def generate():
        yield _sse_event("thinking", {"status": "Understanding your question..."})

        try:
            model_name = _config.get("gemini_model", "gemini-2.5-flash")
            project = _config.get("google_cloud_project", "")
            api_key = _config.get("gemini_api_key", "")

            if project:
                # Vertex AI path (Cloud Run / GCP)
                client = genai.Client(
                    vertexai=True,
                    project=project,
                    location=os.getenv("GOOGLE_CLOUD_LOCATION", "us-central1"),
                )
            elif api_key:
                # API key path (local dev)
                client = genai.Client(api_key=api_key)
            else:
                raise ValueError(
                    "Set GOOGLE_CLOUD_PROJECT (for Vertex AI) or GEMINI_API_KEY (for API key auth)"
                )

            # Build conversation history
            history = _sessions.get(session_id, [])

            contents = []
            for entry in history:
                contents.append(types.Content(**entry))
            contents.append(
                types.Content(role="user", parts=[types.Part.from_text(text=message)])
            )

            tools = [types.Tool(function_declarations=ALL_TOOL_DECLARATIONS)]
            config = types.GenerateContentConfig(
                system_instruction=SYSTEM_INSTRUCTION,
                tools=tools,
                temperature=0.3,
            )

            # Iterative tool-calling loop
            max_iterations = 8
            for _ in range(max_iterations):
                response = await client.aio.models.generate_content(
                    model=model_name,
                    contents=contents,
                    config=config,
                )

                candidate = response.candidates[0] if response.candidates else None
                if (
                    not candidate
                    or not candidate.content
                    or not candidate.content.parts
                ):
                    yield _sse_event(
                        "text", {"content": "I wasn't able to generate a response."}
                    )
                    break

                parts = candidate.content.parts

                # Check if any part is a function call
                function_calls = [p for p in parts if p.function_call]
                if not function_calls:
                    # No function calls — emit text parts and finish
                    for part in parts:
                        if part.text:
                            yield _sse_event("text", {"content": part.text})
                    # Store assistant turn in history
                    contents.append(candidate.content)
                    break

                # Handle function calls
                contents.append(candidate.content)
                function_response_parts = []

                for part in function_calls:
                    fc = part.function_call
                    tool_name = fc.name
                    tool_args = dict(fc.args) if fc.args else {}

                    # create_visualization is special — we intercept it
                    if tool_name == "create_visualization":
                        viz_data = {
                            "type": tool_args.get("viz_type", "table"),
                            "title": tool_args.get("title", ""),
                            "data": tool_args.get("data", {}),
                        }
                        yield _sse_event("visualization", viz_data)
                        function_response_parts.append(
                            types.Part.from_function_response(
                                name=tool_name,
                                response={
                                    "status": "success",
                                    "message": "Visualization displayed to user",
                                },
                            )
                        )
                    else:
                        yield _sse_event(
                            "tool_call",
                            {
                                "tool": tool_name,
                                "description": f"Calling {tool_name}...",
                            },
                        )
                        yield _sse_event(
                            "thinking",
                            {"status": f"Running {tool_name}..."},
                        )

                        handler = TOOL_HANDLERS.get(tool_name)
                        if handler:
                            try:
                                result = await handler(tool_args)
                            except Exception as exc:
                                logger.exception("Tool %s failed", tool_name)
                                result = {"error": str(exc)}
                        else:
                            result = {"error": f"Unknown tool: {tool_name}"}

                        function_response_parts.append(
                            types.Part.from_function_response(
                                name=tool_name,
                                response=result,
                            )
                        )

                contents.append(
                    types.Content(role="user", parts=function_response_parts)
                )
            else:
                yield _sse_event(
                    "text",
                    {
                        "content": "I reached the maximum number of tool calls. Here's what I found so far."
                    },
                )

            # Save session history (only user + model text turns)
            _sessions[session_id] = [
                {"role": c.role, "parts": [{"text": p.text} for p in c.parts if p.text]}
                for c in contents
                if any(p.text for p in c.parts)
            ]

        except Exception as exc:
            logger.exception("Chat stream error")
            yield _sse_event("error", {"message": str(exc)})

        yield _sse_event("done", {"session_id": session_id})

    return StreamingResponse(generate(), media_type="text/event-stream")
