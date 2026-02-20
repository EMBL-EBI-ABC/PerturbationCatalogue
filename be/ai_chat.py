import json
import logging
import os
import uuid
from typing import Any, Dict, List, Optional

import aiohttp

from fastapi import APIRouter
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

from google import genai
from google.genai import types

try:
    from mcp.client.streamable_http import streamablehttp_client
    from mcp import ClientSession

    MCP_AVAILABLE = True
except ImportError:
    MCP_AVAILABLE = False

from data_query import db_pools

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/v1/chat", tags=["AI Chat"])

# --- Configuration ---

_config: Dict[str, Any] = {}

# In-memory session store: session_id -> list of content dicts
_sessions: Dict[str, list] = {}

# Open Targets MCP state (populated at startup)
OT_MCP_DEFAULT_URL = "https://mcp.platform.opentargets.org/mcp"
_ot_tool_declarations: List[types.FunctionDeclaration] = []
_ot_tool_names: set = set()


def configure(
    google_cloud_project: str = "",
    gemini_model: str = "gemini-2.5-flash",
    gemini_api_key: str = "",
    ot_mcp_url: str = "",
):
    _config["google_cloud_project"] = google_cloud_project
    _config["gemini_model"] = gemini_model
    _config["gemini_api_key"] = gemini_api_key
    _config["ot_mcp_url"] = ot_mcp_url or OT_MCP_DEFAULT_URL


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
        "Create a visualization (table, pie chart, bar chart, or gene card) to display in the data portal below the chat. "
        "Always create visualizations when you have data to show the user. "
        "Use pie charts for 2-6 categories, bar charts for comparisons or >6 categories, tables for detailed data. "
        "Use gene_card to display a summary card for a gene with protein info, diseases, domains, and links."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "viz_type": types.Schema(
                type="STRING",
                description="Visualization type: 'table', 'pie_chart', 'bar_chart', or 'gene_card'",
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
                    "For bar_chart: {labels: [...], values: [...], xlabel: '...', ylabel: '...'}. "
                    "For gene_card: {gene_name, protein_name, function, uniprot_id, alphafold_id, diseases: [...], domains: [...], go_terms: [...], subcellular_location}."
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
                    "gene_name": types.Schema(
                        type="STRING", description="Gene symbol for gene_card"
                    ),
                    "protein_name": types.Schema(
                        type="STRING", description="Protein name for gene_card"
                    ),
                    "function": types.Schema(
                        type="STRING",
                        description="Protein function description for gene_card",
                    ),
                    "uniprot_id": types.Schema(
                        type="STRING", description="UniProt accession for gene_card"
                    ),
                    "alphafold_id": types.Schema(
                        type="STRING", description="AlphaFold ID for gene_card"
                    ),
                    "diseases": types.Schema(
                        type="ARRAY",
                        items=types.Schema(type="STRING"),
                        description="Disease associations for gene_card",
                    ),
                    "domains": types.Schema(
                        type="ARRAY",
                        items=types.Schema(type="STRING"),
                        description="Protein domains for gene_card",
                    ),
                    "go_terms": types.Schema(
                        type="ARRAY",
                        items=types.Schema(type="STRING"),
                        description="GO terms for gene_card",
                    ),
                    "subcellular_location": types.Schema(
                        type="STRING",
                        description="Subcellular location for gene_card",
                    ),
                },
            ),
        },
        required=["viz_type", "title", "data"],
    ),
)

LOOKUP_PROTEIN_DECLARATION = types.FunctionDeclaration(
    name="lookup_protein",
    description=(
        "Look up protein information from UniProt for a gene symbol. "
        "Returns function description, subcellular location, disease associations, "
        "protein domains, GO terms, and UniProt accession."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "gene_name": types.Schema(
                type="STRING",
                description="Gene symbol like 'TP53', 'BRCA2', 'KRAS'",
            ),
        },
        required=["gene_name"],
    ),
)

GET_PROTEIN_STRUCTURE_DECLARATION = types.FunctionDeclaration(
    name="get_protein_structure",
    description=(
        "Get AlphaFold predicted protein structure info for a UniProt accession. "
        "Returns structure URLs, confidence scores, and metadata for rendering "
        "an interactive 3D viewer."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "uniprot_id": types.Schema(
                type="STRING",
                description="UniProt accession like 'P04637', 'P51587'",
            ),
        },
        required=["uniprot_id"],
    ),
)

MAP_IDENTIFIERS_DECLARATION = types.FunctionDeclaration(
    name="map_identifiers",
    description=(
        "Map gene symbols, UniProt accessions, or Ensembl gene IDs between identifier systems. "
        "Bridges gene names to UniProt/Ensembl IDs and vice versa. Human proteins only."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "ids": types.Schema(
                type="ARRAY",
                items=types.Schema(type="STRING"),
                description="List of identifiers to map (e.g., ['TP53', 'BRCA2'] or ['P04637'])",
            ),
            "from_db": types.Schema(
                type="STRING",
                description="Source database: 'gene_name', 'uniprot', or 'ensembl'",
            ),
            "to_db": types.Schema(
                type="STRING",
                description="Target database: 'gene_name', 'uniprot', or 'ensembl'",
            ),
        },
        required=["ids", "from_db", "to_db"],
    ),
)

GET_PROTEIN_VARIANTS_DECLARATION = types.FunctionDeclaration(
    name="get_protein_variants",
    description=(
        "Get known protein variants and mutagenesis data from UniProt for a given accession. "
        "Returns natural variants (e.g., disease-associated SNPs) and experimental mutagenesis data. "
        "Particularly useful for interpreting MAVE variant effect scores."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "uniprot_id": types.Schema(
                type="STRING",
                description="UniProt accession like 'P04637', 'P51587'",
            ),
        },
        required=["uniprot_id"],
    ),
)

SEARCH_LITERATURE_DECLARATION = types.FunctionDeclaration(
    name="search_literature",
    description=(
        "Search biomedical literature via Europe PMC. Finds papers about genes, diseases, "
        "perturbation experiments, CRISPR screens, etc. Returns titles, authors, journals, "
        "citation counts, and links. Useful for finding relevant publications about a gene or topic."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "query": types.Schema(
                type="STRING",
                description="Search query (e.g., 'BRCA2 CRISPR screen', 'TP53 perturbation')",
            ),
            "max_results": types.Schema(
                type="INTEGER",
                description="Maximum number of papers to return (default 10, max 20)",
            ),
            "sort": types.Schema(
                type="STRING",
                description="Sort order: 'relevance' (default) or 'date'",
            ),
        },
        required=["query"],
    ),
)

GET_DRUGGABILITY_DECLARATION = types.FunctionDeclaration(
    name="get_druggability",
    description=(
        "Check if a gene/protein target is druggable using the Pharos database (NIH). "
        "Returns the Target Development Level (Tclin=approved drug, Tchem=active compound, "
        "Tbio=biological evidence, Tdark=understudied), protein family, description, "
        "and known drugs/ligands. Use after identifying interesting perturbation targets."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "gene_name": types.Schema(
                type="STRING",
                description="Gene symbol like 'TP53', 'BRCA2', 'KRAS'",
            ),
        },
        required=["gene_name"],
    ),
)

ANNOTATE_VARIANT_DECLARATION = types.FunctionDeclaration(
    name="annotate_variant",
    description=(
        "Annotate a human missense variant with molecular consequence data from ProtVar (EBI). "
        "Returns protein mapping, pathogenicity predictions (AlphaMissense, EVE, ESM-1b, Conservation), "
        "protein stability change (FoldX ddG), CADD score, gnomAD allele frequency, and affected "
        "gene/isoform details. Accepts variants in many formats: dbSNP IDs (rs1042779), "
        "gnomAD (19-1010539-G-C), VCF-like (19 1010539 G C), HGVS genomic (NC_000019.10:g.1010539G>C), "
        "or UniProt+change (P80404 Gln56Arg). Human missense variants only."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "variant": types.Schema(
                type="STRING",
                description=(
                    "Variant identifier in any supported format: "
                    "dbSNP ID (e.g., 'rs1042779'), "
                    "gnomAD format (e.g., '19-1010539-G-C'), "
                    "VCF-like (e.g., '19 1010539 G C'), "
                    "HGVS genomic (e.g., 'NC_000019.10:g.1010539G>C'), "
                    "or UniProt+change (e.g., 'P80404 Gln56Arg')"
                ),
            ),
        },
        required=["variant"],
    ),
)

GET_VARIANT_STRUCTURAL_CONTEXT_DECLARATION = types.FunctionDeclaration(
    name="get_variant_structural_context",
    description=(
        "Get structural context for a specific residue position in a protein from ProtVar (EBI). "
        "Returns PDB structure mappings, predicted binding pockets, protein-protein interaction "
        "interfaces, FoldX stability predictions for all possible substitutions, and functional "
        "annotations (domains, active sites, PTMs) at that position. "
        "Use after identifying a variant of interest to understand its structural impact."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "uniprot_id": types.Schema(
                type="STRING",
                description="UniProt accession like 'P04637', 'Q9NUW8'",
            ),
            "position": types.Schema(
                type="INTEGER",
                description="Amino acid residue position (1-based)",
            ),
            "variant_aa": types.Schema(
                type="STRING",
                description="Optional: variant amino acid (1- or 3-letter code, e.g., 'R' or 'Arg') for specific substitution predictions",
            ),
        },
        required=["uniprot_id", "position"],
    ),
)

INTERNAL_TOOL_DECLARATIONS = [
    SEARCH_DATASETS_DECLARATION,
    SEARCH_TARGET_SUMMARY_DECLARATION,
    QUERY_PERTURBATION_DATA_DECLARATION,
    GET_CATALOGUE_SUMMARY_DECLARATION,
    CREATE_VISUALIZATION_DECLARATION,
    LOOKUP_PROTEIN_DECLARATION,
    GET_PROTEIN_STRUCTURE_DECLARATION,
    MAP_IDENTIFIERS_DECLARATION,
    GET_PROTEIN_VARIANTS_DECLARATION,
    SEARCH_LITERATURE_DECLARATION,
    GET_DRUGGABILITY_DECLARATION,
    ANNOTATE_VARIANT_DECLARATION,
    GET_VARIANT_STRUCTURAL_CONTEXT_DECLARATION,
]


def _get_all_tool_declarations():
    return INTERNAL_TOOL_DECLARATIONS + _ot_tool_declarations

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


async def _tool_lookup_protein(args: dict) -> dict:
    gene_name = args.get("gene_name", "")
    if not gene_name:
        return {"error": "gene_name is required"}

    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene_exact:{gene_name} AND organism_id:9606 AND reviewed:true",
        "fields": (
            "accession,gene_primary,protein_name,cc_function,"
            "cc_subcellular_location,cc_disease,ft_domain,"
            "go_p,go_f,go_c,xref_pdb,xref_alphafolddb"
        ),
        "format": "json",
        "size": "1",
    }
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=15)) as resp:
                if resp.status != 200:
                    return {"error": f"UniProt API returned status {resp.status}"}
                data = await resp.json()
    except Exception as exc:
        return {"error": f"UniProt API error: {str(exc)}"}

    results = data.get("results", [])
    if not results:
        return {"error": f"No UniProt entry found for gene '{gene_name}' in human"}

    entry = results[0]

    # Extract accession
    uniprot_id = entry.get("primaryAccession", "")

    # Extract gene name
    gene = ""
    genes = entry.get("genes", [])
    if genes:
        gene = genes[0].get("geneName", {}).get("value", "")

    # Extract protein name
    protein_name = ""
    pn = entry.get("proteinDescription", {})
    rec_name = pn.get("recommendedName")
    if rec_name:
        protein_name = rec_name.get("fullName", {}).get("value", "")
    elif pn.get("submissionNames"):
        protein_name = pn["submissionNames"][0].get("fullName", {}).get("value", "")

    # Extract function
    function_text = ""
    comments = entry.get("comments", [])
    for c in comments:
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts", [])
            if texts:
                function_text = texts[0].get("value", "")
            break

    # Extract subcellular location
    subcellular_location = ""
    for c in comments:
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            locs = c.get("subcellularLocations", [])
            loc_names = []
            for loc in locs:
                loc_val = loc.get("location", {}).get("value", "")
                if loc_val:
                    loc_names.append(loc_val)
            subcellular_location = "; ".join(loc_names)
            break

    # Extract diseases
    diseases = []
    for c in comments:
        if c.get("commentType") == "DISEASE":
            disease = c.get("disease", {})
            disease_name = disease.get("diseaseId", "")
            if disease_name:
                diseases.append(disease_name)

    # Extract domains
    domains = []
    features = entry.get("features", [])
    for f in features:
        if f.get("type") == "Domain":
            desc = f.get("description", "")
            if desc and desc not in domains:
                domains.append(desc)

    # Extract GO terms
    go_terms = []
    xrefs = entry.get("uniProtKBCrossReferences", [])
    for xref in xrefs:
        if xref.get("database") == "GO":
            props = xref.get("properties", [])
            for prop in props:
                if prop.get("key") == "GoTerm":
                    term = prop.get("value", "")
                    # GO terms look like "P:apoptotic process" — strip prefix
                    if ":" in term:
                        term = term.split(":", 1)[1]
                    go_terms.append(term)

    # Extract PDB IDs
    pdb_ids = []
    for xref in xrefs:
        if xref.get("database") == "PDB":
            pdb_ids.append(xref.get("id", ""))

    # Extract AlphaFold ID
    alphafold_id = ""
    for xref in xrefs:
        if xref.get("database") == "AlphaFoldDB":
            alphafold_id = xref.get("id", "")
            break

    return {
        "uniprot_id": uniprot_id,
        "gene_name": gene or gene_name,
        "protein_name": protein_name,
        "function": function_text,
        "subcellular_location": subcellular_location,
        "diseases": diseases[:10],
        "domains": domains[:10],
        "go_terms": go_terms[:15],
        "pdb_ids": pdb_ids[:5],
        "alphafold_id": alphafold_id,
    }


async def _tool_get_protein_structure(args: dict) -> dict:
    uniprot_id = args.get("uniprot_id", "")
    if not uniprot_id:
        return {"error": "uniprot_id is required"}

    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=15)) as resp:
                if resp.status == 404:
                    return {"error": f"No AlphaFold structure found for {uniprot_id}"}
                if resp.status != 200:
                    return {"error": f"AlphaFold API returned status {resp.status}"}
                data = await resp.json()
    except Exception as exc:
        return {"error": f"AlphaFold API error: {str(exc)}"}

    # API returns a list; take the first entry
    entry = data[0] if isinstance(data, list) and data else data

    return {
        "entry_id": entry.get("entryId", ""),
        "gene": entry.get("gene", ""),
        "organism": entry.get("organismScientificName", ""),
        "uniprot_id": entry.get("uniprotAccession", uniprot_id),
        "pdb_url": entry.get("pdbUrl", ""),
        "cif_url": entry.get("cifUrl", ""),
        "pae_image_url": entry.get("paeImageUrl", ""),
        "plddt_url": entry.get("confidenceUrl", ""),
    }


async def _tool_map_identifiers(args: dict) -> dict:
    ids = args.get("ids", [])
    from_db = args.get("from_db", "gene_name")
    to_db = args.get("to_db", "uniprot")

    if not ids:
        return {"error": "ids list is required"}

    ids = ids[:20]
    results = []
    fields = "accession,gene_primary,xref_ensembl"

    async with aiohttp.ClientSession() as http:
        for id_val in ids:
            if from_db == "gene_name":
                query = f"gene_exact:{id_val} AND organism_id:9606 AND reviewed:true"
            elif from_db == "uniprot":
                query = f"accession:{id_val}"
            elif from_db == "ensembl":
                query = f"xref:ensembl-{id_val} AND organism_id:9606 AND reviewed:true"
            else:
                return {
                    "error": f"Unknown from_db: {from_db}. Use 'gene_name', 'uniprot', or 'ensembl'"
                }

            url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": query,
                "fields": fields,
                "format": "json",
                "size": "1",
            }

            try:
                async with http.get(
                    url, params=params, timeout=aiohttp.ClientTimeout(total=15)
                ) as resp:
                    if resp.status != 200:
                        results.append({"input": id_val, "error": f"HTTP {resp.status}"})
                        continue
                    data = await resp.json()
            except Exception as exc:
                results.append({"input": id_val, "error": str(exc)})
                continue

            entries = data.get("results", [])
            if not entries:
                results.append({"input": id_val, "error": "Not found"})
                continue

            entry = entries[0]
            mapping = {"input": id_val}
            mapping["uniprot_id"] = entry.get("primaryAccession", "")

            genes = entry.get("genes", [])
            mapping["gene_name"] = (
                genes[0].get("geneName", {}).get("value", "") if genes else ""
            )

            ensembl_ids = []
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "Ensembl":
                    for prop in xref.get("properties", []):
                        if prop.get("key") == "GeneId":
                            eid = prop.get("value", "")
                            if eid and eid not in ensembl_ids:
                                ensembl_ids.append(eid)
            mapping["ensembl_ids"] = ensembl_ids

            results.append(mapping)

    return {"mappings": results}


async def _tool_get_protein_variants(args: dict) -> dict:
    uniprot_id = args.get("uniprot_id", "")
    if not uniprot_id:
        return {"error": "uniprot_id is required"}

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    params = {
        "fields": "ft_variant,ft_mutagen,gene_primary,protein_name",
        "format": "json",
    }

    try:
        async with aiohttp.ClientSession() as http:
            async with http.get(
                url, params=params, timeout=aiohttp.ClientTimeout(total=15)
            ) as resp:
                if resp.status == 404:
                    return {"error": f"No UniProt entry found for {uniprot_id}"}
                if resp.status != 200:
                    return {"error": f"UniProt API returned status {resp.status}"}
                entry = await resp.json()
    except Exception as exc:
        return {"error": f"UniProt API error: {str(exc)}"}

    gene_name = ""
    genes = entry.get("genes", [])
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")

    variants = []
    mutagenesis = []

    for feat in entry.get("features", []):
        if feat.get("type") == "Natural variant":
            v = {
                "position": feat.get("location", {}).get("start", {}).get("value"),
                "original": feat.get("alternativeSequence", {}).get(
                    "originalSequence", ""
                ),
                "variant": ", ".join(
                    feat.get("alternativeSequence", {}).get(
                        "alternativeSequences", []
                    )
                ),
                "description": feat.get("description", ""),
            }
            if feat.get("featureId"):
                v["feature_id"] = feat["featureId"]
            variants.append(v)
        elif feat.get("type") == "Mutagenesis":
            m = {
                "position": feat.get("location", {}).get("start", {}).get("value"),
                "original": feat.get("alternativeSequence", {}).get(
                    "originalSequence", ""
                ),
                "variant": ", ".join(
                    feat.get("alternativeSequence", {}).get(
                        "alternativeSequences", []
                    )
                ),
                "description": feat.get("description", ""),
            }
            mutagenesis.append(m)

    return {
        "uniprot_id": uniprot_id,
        "gene_name": gene_name,
        "natural_variants": variants[:50],
        "mutagenesis": mutagenesis[:50],
        "total_variants": len(variants),
        "total_mutagenesis": len(mutagenesis),
    }


async def _tool_search_literature(args: dict) -> dict:
    query = args.get("query", "")
    if not query:
        return {"error": "query is required"}

    max_results = min(args.get("max_results", 10), 20)
    sort = args.get("sort", "relevance")

    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": query,
        "format": "json",
        "pageSize": str(max_results),
        "resultType": "lite",
    }
    if sort == "date":
        params["sort"] = "P_PDATE_D desc"

    try:
        async with aiohttp.ClientSession() as http:
            async with http.get(
                url, params=params, timeout=aiohttp.ClientTimeout(total=15)
            ) as resp:
                if resp.status != 200:
                    return {"error": f"Europe PMC API returned status {resp.status}"}
                data = await resp.json()
    except Exception as exc:
        return {"error": f"Europe PMC API error: {str(exc)}"}

    result_list = data.get("resultList", {}).get("result", [])
    hit_count = data.get("hitCount", 0)

    articles = []
    for article in result_list:
        a = {
            "title": article.get("title", ""),
            "authors": article.get("authorString", ""),
            "journal": article.get("journalTitle", ""),
            "year": article.get("pubYear", ""),
            "cited_by": article.get("citedByCount", 0),
            "pmid": article.get("pmid", ""),
            "doi": article.get("doi", ""),
        }
        if a["pmid"]:
            a["url"] = f"https://europepmc.org/article/MED/{a['pmid']}"
        elif a["doi"]:
            a["url"] = f"https://doi.org/{a['doi']}"

        abstract = article.get("abstractText", "")
        if abstract:
            a["abstract"] = (
                abstract[:300] + ("..." if len(abstract) > 300 else "")
            )

        articles.append(a)

    return {
        "total_hits": hit_count,
        "articles": articles,
    }


async def _tool_get_druggability(args: dict) -> dict:
    gene_name = args.get("gene_name", "")
    if not gene_name:
        return {"error": "gene_name is required"}

    url = "https://pharos-api.ncats.io/graphql"
    payload = {
        "query": """
        query TargetDruggability($term: String!) {
          targets(filter: {term: $term}, top: 10) {
            targets {
              name
              sym
              tdl
              fam
              description
              novelty
            }
          }
        }
        """,
        "variables": {"term": gene_name.upper()},
    }

    try:
        async with aiohttp.ClientSession() as http:
            async with http.post(
                url,
                json=payload,
                timeout=aiohttp.ClientTimeout(total=30),
                headers={"Content-Type": "application/json"},
            ) as resp:
                if resp.status != 200:
                    body = await resp.text()
                    return {"error": f"Pharos API returned status {resp.status}: {body[:200]}"}
                data = await resp.json()
    except Exception as exc:
        return {"error": f"Pharos API error: {str(exc)}"}

    if "errors" in data:
        return {
            "error": f"Pharos GraphQL error: {data['errors'][0].get('message', '')}"
        }

    targets = data.get("data", {}).get("targets", {}).get("targets", [])
    if not targets:
        return {"error": f"No Pharos entry found for gene '{gene_name}'"}

    # Pharos term search is fuzzy — find exact gene symbol match
    exact = [t for t in targets if t.get("sym", "").upper() == gene_name.upper()]
    if exact:
        targets = exact

    target = targets[0]
    tdl = target.get("tdl", "")
    tdl_desc = {
        "Tclin": "Approved drug target — has at least one approved drug",
        "Tchem": "Chemical tool target — has active compounds but no approved drug",
        "Tbio": "Biological target — has biological evidence but no active compounds",
        "Tdark": "Understudied target — little known about this protein",
    }

    return {
        "gene_name": target.get("sym", gene_name),
        "protein_name": target.get("name", ""),
        "tdl": tdl,
        "tdl_description": tdl_desc.get(tdl, "Unknown"),
        "protein_family": target.get("fam", ""),
        "description": target.get("description", ""),
        "novelty_score": target.get("novelty"),
        "pharos_url": f"https://pharos.nih.gov/targets/{target.get('sym', gene_name)}",
    }


PROTVAR_BASE = "https://www.ebi.ac.uk/ProtVar/api"


async def _tool_annotate_variant(args: dict) -> dict:
    """Annotate a human missense variant using ProtVar mapping + scores + FoldX."""
    variant = args.get("variant", "").strip()
    if not variant:
        return {"error": "variant is required"}

    # Step 1: Map the variant via ProtVar
    url = f"{PROTVAR_BASE}/mapping"
    params = {"input": variant}

    try:
        async with aiohttp.ClientSession() as http:
            async with http.get(
                url, params=params, timeout=aiohttp.ClientTimeout(total=20)
            ) as resp:
                if resp.status != 200:
                    body = await resp.text()
                    return {
                        "error": f"ProtVar mapping API returned status {resp.status}: {body[:300]}"
                    }
                data = await resp.json()
    except Exception as exc:
        return {"error": f"ProtVar API error: {str(exc)}"}

    inputs = data.get("inputs", [])
    if not inputs:
        return {"error": f"No ProtVar results for variant '{variant}'"}

    user_input = inputs[0]
    input_type = user_input.get("type", "")

    # Handle ID-based inputs (dbSNP, ClinVar, COSMIC) which have derivedGenomicInputs
    genomic_inputs = user_input.get("derivedGenomicInputs", [])
    if genomic_inputs:
        # Use the first derived genomic input for mapping
        user_input = genomic_inputs[0]

    mappings = user_input.get("mappings", [])
    if not mappings:
        messages = data.get("messages", [])
        msg_texts = [m.get("text", "") for m in messages if m.get("text")]
        return {
            "error": f"No gene mappings found for variant '{variant}'",
            "messages": msg_texts or None,
        }

    # Extract results from the first gene mapping
    results = []
    for mapping in mappings[:3]:  # Limit to top 3 gene mappings
        genes = mapping.get("genes", [])
        for gene in genes[:2]:
            gene_name = gene.get("geneName", "")
            ensg = gene.get("ensg", "")
            cadd_score = gene.get("caddScore")
            allele_freq = gene.get("alleleFreq")

            isoforms = gene.get("isoforms", [])
            for iso in isoforms[:2]:
                accession = iso.get("accession", "")
                canonical = iso.get("canonical", False)
                position = iso.get("isoformPosition")
                ref_aa = iso.get("refAA", "")
                var_aa = iso.get("variantAA", "")
                consequences = iso.get("consequences", [])
                codon_change = iso.get("codonChange", "")
                aa_change = iso.get("aminoAcidChange", "")

                # Inline scores from mapping response
                conserv_score = iso.get("conservScore")
                am_score = iso.get("amScore")
                esm_score = iso.get("esmScore")

                entry = {
                    "gene": gene_name,
                    "ensembl_gene": ensg,
                    "uniprot_accession": accession,
                    "canonical_isoform": canonical,
                    "position": position,
                    "ref_aa": ref_aa,
                    "variant_aa": var_aa,
                    "consequences": consequences,
                    "codon_change": codon_change,
                    "amino_acid_change": aa_change,
                    "cadd_score": cadd_score,
                    "gnomad_allele_frequency": allele_freq,
                }

                # Add inline scores
                predictions = {}
                if conserv_score is not None:
                    predictions["conservation"] = conserv_score
                if am_score is not None:
                    predictions["alphamissense"] = am_score
                if esm_score is not None:
                    predictions["esm1b"] = esm_score

                # Step 2: Fetch EVE score and full score details if we have accession + position + variant_aa
                if accession and position and var_aa:
                    try:
                        score_url = f"{PROTVAR_BASE}/score/{accession}/{position}"
                        score_params = {"mt": var_aa}
                        async with aiohttp.ClientSession() as http:
                            async with http.get(
                                score_url,
                                params=score_params,
                                timeout=aiohttp.ClientTimeout(total=10),
                            ) as score_resp:
                                if score_resp.status == 200:
                                    scores = await score_resp.json()
                                    for s in scores:
                                        name = s.get("name", "")
                                        if name == "EVE":
                                            predictions["eve_score"] = s.get("score")
                                            predictions["eve_class"] = s.get(
                                                "eveClass"
                                            )
                                        elif name == "AM":
                                            predictions["alphamissense"] = s.get(
                                                "amPathogenicity"
                                            )
                                            predictions["alphamissense_class"] = (
                                                s.get("amClass")
                                            )
                                        elif name == "CONSERV":
                                            predictions["conservation"] = s.get(
                                                "score"
                                            )
                                        elif name == "ESM":
                                            predictions["esm1b"] = s.get("score")
                    except Exception:
                        pass  # Scores are supplementary; don't fail the whole call

                    # Step 3: Fetch FoldX stability prediction
                    try:
                        foldx_url = (
                            f"{PROTVAR_BASE}/foldx/{accession}/{position}"
                        )
                        foldx_params = {"variantAA": var_aa}
                        async with aiohttp.ClientSession() as http:
                            async with http.get(
                                foldx_url,
                                params=foldx_params,
                                timeout=aiohttp.ClientTimeout(total=10),
                            ) as foldx_resp:
                                if foldx_resp.status == 200:
                                    foldx_data = await foldx_resp.json()
                                    if foldx_data:
                                        fx = foldx_data[0]
                                        predictions["foldx_ddg"] = fx.get(
                                            "foldxDdg"
                                        )
                                        predictions["alphafold_plddt"] = fx.get(
                                            "plddt"
                                        )
                                        # Interpret stability
                                        ddg = fx.get("foldxDdg")
                                        if ddg is not None:
                                            if ddg > 2:
                                                predictions[
                                                    "stability_effect"
                                                ] = "Destabilizing"
                                            elif ddg > 0.5:
                                                predictions[
                                                    "stability_effect"
                                                ] = "Mildly destabilizing"
                                            elif ddg < -2:
                                                predictions[
                                                    "stability_effect"
                                                ] = "Stabilizing"
                                            else:
                                                predictions[
                                                    "stability_effect"
                                                ] = "Neutral"
                    except Exception:
                        pass  # FoldX is supplementary

                if predictions:
                    entry["predictions"] = predictions

                results.append(entry)

    # Extract any messages from ProtVar
    messages = data.get("messages", [])
    msg_texts = [m.get("text", "") for m in messages if m.get("text")]

    output = {
        "input_variant": variant,
        "annotations": results,
        "protvar_url": f"https://www.ebi.ac.uk/ProtVar/query?search={variant}",
    }
    if msg_texts:
        output["messages"] = msg_texts

    return output


async def _tool_get_variant_structural_context(args: dict) -> dict:
    """Get structural context for a residue position from ProtVar."""
    uniprot_id = args.get("uniprot_id", "").strip()
    position = args.get("position")
    variant_aa = args.get("variant_aa", "")

    if not uniprot_id:
        return {"error": "uniprot_id is required"}
    if not position:
        return {"error": "position is required"}

    result = {
        "uniprot_id": uniprot_id,
        "position": position,
    }

    async with aiohttp.ClientSession() as http:
        # 1. Functional annotations (domains, active sites, PTMs)
        try:
            func_url = f"{PROTVAR_BASE}/function/{uniprot_id}/{position}"
            func_params = {}
            if variant_aa:
                func_params["variantAA"] = variant_aa
            async with http.get(
                func_url,
                params=func_params if func_params else None,
                timeout=aiohttp.ClientTimeout(total=15),
            ) as resp:
                if resp.status == 200:
                    func_data = await resp.json()
                    # Extract key functional features
                    features = []
                    for feat in func_data.get("features", []):
                        f = {
                            "type": feat.get("type", ""),
                            "description": feat.get("description", ""),
                        }
                        loc = feat.get("location", {})
                        if loc:
                            f["start"] = loc.get("start", {}).get("value")
                            f["end"] = loc.get("end", {}).get("value")
                        if f["type"]:
                            features.append(f)

                    if features:
                        result["functional_features"] = features[:20]

                    # Gene/protein info
                    gene_name = ""
                    genes = func_data.get("genes", [])
                    if genes:
                        gene_name = (
                            genes[0].get("geneName", {}).get("value", "")
                        )
                    if gene_name:
                        result["gene_name"] = gene_name

                    protein_name = (
                        func_data.get("proteinDescription", {})
                        .get("recommendedName", {})
                        .get("fullName", {})
                        .get("value", "")
                    )
                    if protein_name:
                        result["protein_name"] = protein_name

                    # Pockets from function endpoint
                    pockets = func_data.get("pockets", [])
                    if pockets:
                        pocket_info = []
                        for p in pockets[:5]:
                            pocket_info.append(
                                {
                                    "pocket_id": p.get("pocketId"),
                                    "score": p.get("score"),
                                    "energy_per_vol": p.get("energyPerVol"),
                                    "buriedness": p.get("buriedness"),
                                    "mean_plddt": p.get("meanPlddt"),
                                }
                            )
                        result["binding_pockets"] = pocket_info

                    # Interactions from function endpoint
                    interactions = func_data.get("interactions", [])
                    if interactions:
                        interaction_info = []
                        for inter in interactions[:5]:
                            interaction_info.append(
                                {
                                    "partner_a": inter.get("a", ""),
                                    "partner_b": inter.get("b", ""),
                                    "pdockq": inter.get("pdockq"),
                                }
                            )
                        result["protein_interactions"] = interaction_info

                    # FoldX data from function endpoint
                    foldx_list = func_data.get("foldxs", [])
                    if foldx_list:
                        foldx_info = []
                        for fx in foldx_list[:20]:
                            foldx_info.append(
                                {
                                    "wild_type": fx.get("wildType", ""),
                                    "mutated_type": fx.get("mutatedType", ""),
                                    "foldx_ddg": fx.get("foldxDdg"),
                                    "plddt": fx.get("plddt"),
                                }
                            )
                        result["foldx_predictions"] = foldx_info
                elif resp.status == 404:
                    return {
                        "error": f"No ProtVar data for {uniprot_id} position {position}"
                    }
                else:
                    result["function_error"] = (
                        f"Function endpoint returned status {resp.status}"
                    )
        except Exception as exc:
            result["function_error"] = str(exc)

        # 2. PDB structure mappings
        try:
            struct_url = f"{PROTVAR_BASE}/structure/{uniprot_id}/{position}"
            async with http.get(
                struct_url, timeout=aiohttp.ClientTimeout(total=10)
            ) as resp:
                if resp.status == 200:
                    struct_data = await resp.json()
                    if struct_data:
                        structures = []
                        for s in struct_data[:10]:
                            structures.append(
                                {
                                    "pdb_id": s.get("pdb_id", ""),
                                    "chain_id": s.get("chain_id", ""),
                                    "experimental_method": s.get(
                                        "experimental_method", ""
                                    ),
                                    "resolution": s.get("resolution"),
                                }
                            )
                        if structures:
                            result["pdb_structures"] = structures
        except Exception:
            pass  # Structural data is supplementary

        # 3. Co-located variants (population data)
        try:
            pop_url = f"{PROTVAR_BASE}/population/{uniprot_id}/{position}"
            async with http.get(
                pop_url, timeout=aiohttp.ClientTimeout(total=10)
            ) as resp:
                if resp.status == 200:
                    pop_data = await resp.json()
                    colocated = []
                    for var in pop_data.get(
                        "proteinColocatedVariant", []
                    )[:10]:
                        v = {
                            "wild_type": var.get("wildType", ""),
                            "variant": var.get("alternativeSequence", ""),
                        }
                        # Clinical significance
                        clin_sigs = var.get("clinicalSignificances", [])
                        if clin_sigs:
                            v["clinical_significance"] = clin_sigs[0]
                        # Disease associations
                        assocs = var.get("association", [])
                        if assocs:
                            v["disease_associations"] = [
                                a.get("name", "") for a in assocs[:3]
                            ]
                        # dbSNP xrefs
                        xrefs = var.get("xrefs", [])
                        for xref in xrefs:
                            if xref.get("source") == "dbSNP":
                                v["dbsnp_id"] = xref.get("id", "")
                                break
                        colocated.append(v)
                    if colocated:
                        result["colocated_variants"] = colocated
        except Exception:
            pass  # Population data is supplementary

        # 4. Pathogenicity scores (if variant_aa specified)
        if variant_aa:
            try:
                score_url = f"{PROTVAR_BASE}/score/{uniprot_id}/{position}"
                score_params = {"mt": variant_aa}
                async with http.get(
                    score_url,
                    params=score_params,
                    timeout=aiohttp.ClientTimeout(total=10),
                ) as resp:
                    if resp.status == 200:
                        scores = await resp.json()
                        predictions = {}
                        for s in scores:
                            name = s.get("name", "")
                            if name == "EVE":
                                predictions["eve_score"] = s.get("score")
                                predictions["eve_class"] = s.get("eveClass")
                            elif name == "AM":
                                predictions["alphamissense"] = s.get(
                                    "amPathogenicity"
                                )
                                predictions["alphamissense_class"] = s.get(
                                    "amClass"
                                )
                            elif name == "CONSERV":
                                predictions["conservation"] = s.get("score")
                            elif name == "ESM":
                                predictions["esm1b"] = s.get("score")
                        if predictions:
                            result["pathogenicity_scores"] = predictions
            except Exception:
                pass

    result["protvar_url"] = (
        f"https://www.ebi.ac.uk/ProtVar/query?search={uniprot_id}+{position}"
    )

    return result


TOOL_HANDLERS = {
    "search_datasets": _tool_search_datasets,
    "search_target_summary": _tool_search_target_summary,
    "query_perturbation_data": _tool_query_perturbation_data,
    "get_catalogue_summary": lambda args: _tool_get_catalogue_summary(),
    "lookup_protein": _tool_lookup_protein,
    "get_protein_structure": _tool_get_protein_structure,
    "map_identifiers": _tool_map_identifiers,
    "get_protein_variants": _tool_get_protein_variants,
    "search_literature": _tool_search_literature,
    "get_druggability": _tool_get_druggability,
    "annotate_variant": _tool_annotate_variant,
    "get_variant_structural_context": _tool_get_variant_structural_context,
}


# --- Open Targets MCP integration ---


def _json_schema_to_gemini(schema: dict) -> types.Schema:
    """Convert a JSON Schema dict to a Gemini types.Schema.

    Handles features Gemini doesn't support (anyOf, additionalProperties)
    by simplifying to the closest Gemini-compatible representation.
    """
    type_map = {
        "string": "STRING",
        "integer": "INTEGER",
        "number": "NUMBER",
        "boolean": "BOOLEAN",
        "array": "ARRAY",
        "object": "OBJECT",
    }

    # Handle anyOf / oneOf: pick the first non-null variant
    for union_key in ("anyOf", "oneOf"):
        if union_key in schema:
            variants = [v for v in schema[union_key] if v.get("type") != "null"]
            base = variants[0] if variants else {"type": "string"}
            # Carry over description and default from parent
            if "description" in schema:
                base.setdefault("description", schema["description"])
            return _json_schema_to_gemini(base)

    json_type = schema.get("type", "string")
    kwargs: Dict[str, Any] = {"type": type_map.get(json_type, "STRING")}

    if "description" in schema:
        kwargs["description"] = schema["description"]

    if "enum" in schema:
        kwargs["enum"] = schema["enum"]

    if json_type == "object" and "properties" in schema:
        kwargs["properties"] = {
            k: _json_schema_to_gemini(v) for k, v in schema["properties"].items()
        }
        if "required" in schema:
            kwargs["required"] = schema["required"]

    if json_type == "array" and "items" in schema:
        kwargs["items"] = _json_schema_to_gemini(schema["items"])

    return types.Schema(**kwargs)


async def init_open_targets_mcp():
    """Connect to OT MCP server, discover tools, and build Gemini declarations."""
    if not MCP_AVAILABLE:
        logger.warning("MCP SDK not installed — Open Targets tools will be unavailable")
        return

    url = _config.get("ot_mcp_url", OT_MCP_DEFAULT_URL)
    logger.info("Discovering Open Targets MCP tools from %s", url)

    try:
        async with streamablehttp_client(url) as (read_stream, write_stream, _):
            async with ClientSession(read_stream, write_stream) as session:
                await session.initialize()
                tools_result = await session.list_tools()

                # Skip schema tool — response is too large for LLM context;
                # example queries are provided in the system prompt instead.
                skip_tools = {"get_open_targets_graphql_schema"}

                for tool in tools_result.tools:
                    if tool.name in skip_tools:
                        logger.info("  Skipped OT tool: %s (too large)", tool.name)
                        continue
                    input_schema = tool.inputSchema or {"type": "object", "properties": {}}
                    gemini_params = _json_schema_to_gemini(input_schema)
                    decl = types.FunctionDeclaration(
                        name=tool.name,
                        description=tool.description or "",
                        parameters=gemini_params,
                    )
                    _ot_tool_declarations.append(decl)
                    _ot_tool_names.add(tool.name)
                    logger.info("  Registered OT tool: %s", tool.name)

        logger.info(
            "Open Targets MCP: %d tools available", len(_ot_tool_declarations)
        )
    except Exception:
        logger.exception("Failed to connect to Open Targets MCP — tools will be unavailable")


_OT_TOOL_DESCRIPTIONS = {
    "search_entities": "Searching Open Targets...",
    "query_open_targets_graphql": "Querying Open Targets Platform...",
    "get_open_targets_graphql_schema": "Fetching Open Targets schema...",
    "batch_query_open_targets_graphql": "Running batch query on Open Targets...",
    "lookup_protein": "Looking up protein information from UniProt...",
    "get_protein_structure": "Loading 3D protein structure...",
    "search_datasets": "Searching datasets...",
    "search_target_summary": "Looking up gene targets...",
    "query_perturbation_data": "Querying perturbation data...",
    "get_catalogue_summary": "Getting catalogue summary...",
    "map_identifiers": "Mapping identifiers across databases...",
    "get_protein_variants": "Fetching protein variants from UniProt...",
    "search_literature": "Searching biomedical literature...",
    "get_druggability": "Checking druggability on Pharos...",
    "annotate_variant": "Annotating variant with ProtVar...",
    "get_variant_structural_context": "Fetching structural context from ProtVar...",
}


async def _call_open_targets_tool(name: str, args: dict) -> dict:
    """Execute a tool call against the remote OT MCP server."""
    url = _config.get("ot_mcp_url", OT_MCP_DEFAULT_URL)
    try:
        async with streamablehttp_client(url) as (read_stream, write_stream, _):
            async with ClientSession(read_stream, write_stream) as session:
                await session.initialize()
                result = await session.call_tool(name, args)

                if result.isError:
                    error_text = " ".join(
                        getattr(c, "text", str(c)) for c in result.content
                    )
                    return {"error": f"Open Targets tool error: {error_text}"}

                # Combine text content from the result
                text_parts = []
                for item in result.content:
                    if hasattr(item, "text"):
                        text_parts.append(item.text)

                combined = "\n".join(text_parts)
                # Try to parse as JSON for structured data
                try:
                    return json.loads(combined)
                except (json.JSONDecodeError, TypeError):
                    # Truncate very large text to avoid blowing up context
                    if len(combined) > 30000:
                        combined = combined[:30000] + "\n... (truncated)"
                    return {"result": combined}
    except Exception as exc:
        logger.exception("Open Targets MCP call failed: %s", name)
        return {"error": f"Open Targets MCP error: {str(exc)}"}


# --- System instruction ---

SYSTEM_INSTRUCTION = """You are the AI Explorer for the Perturbation Catalogue, a comprehensive resource for human gene perturbation experiments (CRISPR screens, Perturb-seq, MAVE).

Your role:
- Help researchers explore perturbation data through natural language
- Search datasets, look up gene targets, and query experimental results
- Connect perturbation findings to broader biological context using Open Targets
- Be concise and scientific, lead with key findings

CRITICAL RULES FOR RESPONSES:
- NEVER include markdown tables, data listings, or raw data rows in your text responses.
- Your text response should ONLY contain a brief summary and interpretation of the findings.
- ALL data (tables, charts) MUST be sent via the create_visualization tool, which displays them in a separate data portal below the chat.
- After calling create_visualization, write a short summary (2-3 sentences) interpreting the results. Do NOT repeat the data in text.

CRITICAL TOOL ROUTING — follow these rules for tool selection:
- "Is X druggable?" / "drugs for X" / "target X" → call get_druggability FIRST (Pharos), NOT Open Targets
- "What papers..." / "literature on..." → call search_literature (Europe PMC)
- "What diseases are linked to X?" → use Open Targets
- "Show structure of X" → call get_protein_structure (AlphaFold)
- "What variants does X have?" → call get_protein_variants (UniProt)
- "What does variant rs123 do?" / "Is this variant pathogenic?" / "Effect of mutation X" → call annotate_variant (ProtVar)
- "Structural impact at position X" / "What's at residue 493?" → call get_variant_structural_context (ProtVar)

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

When users ask about a gene, first search for target summary to understand what's available, then query specific data if needed.

OPEN TARGETS INTEGRATION:
You also have access to the Open Targets Platform via its MCP tools. Open Targets aggregates data from 22+ sources (GWAS, ClinVar, ChEMBL, UniProt, Reactome, etc.) for target-disease associations, drug data, and genetic evidence.

Open Targets workflow:
1. Use search_entities to resolve gene/disease/drug names to standardized IDs (Ensembl for genes, EFO/MONDO for diseases, ChEMBL for drugs). The result contains objects with "id" and "entity" fields.
2. Use query_open_targets_graphql with the resolved IDs. Use the example queries below as templates — do NOT call get_open_targets_graphql_schema (it is too large).

IMPORTANT: Always use the GraphQL query templates below. Adapt them as needed but keep the structure.

Disease associations for a target (gene):
  query_string: "query { target(ensemblId: \\"ENSG00000139618\\") { approvedSymbol associatedDiseases(page: {index: 0, size: 25}) { count rows { disease { id name } score } } } }"

Target associations for a disease:
  query_string: "query { disease(efoId: \\"MONDO_0007254\\") { name associatedTargets(page: {index: 0, size: 25}) { count rows { target { id approvedSymbol } score } } } }"

Drugs for a target (ONLY use AFTER calling get_druggability first — never as the first tool for druggability questions):
  query_string: "query { target(ensemblId: \\"ENSG00000139618\\") { approvedSymbol knownDrugs(size: 25) { count rows { drug { id name } mechanismOfAction phase status } } } }"

Target details (pathways, GO terms):
  query_string: "query { target(ensemblId: \\"ENSG00000139618\\") { approvedSymbol biotype functionDescriptions pathways { pathway pathwayId } } }"

Use Open Targets when users ask about:
- Disease associations for a gene (e.g. "What diseases are linked to BRCA2?")
- Genetic evidence and GWAS associations
- Known pathways and biological functions from curated databases
- Detailed drug mechanisms of action and clinical trial phases (AFTER checking druggability with get_druggability first)

When combining data from both our Catalogue and Open Targets, clearly distinguish between the two sources in your response.

PROTEIN CONTEXT TOOLS:
You have access to UniProt and AlphaFold tools for protein-level information.

- lookup_protein: Look up protein function, domains, disease associations, and GO terms from UniProt for any gene. Use this when users ask "What does gene X do?" or when you want to provide biological context for a perturbation target.
- get_protein_structure: Fetch AlphaFold predicted 3D structure. This automatically renders an interactive 3D viewer in the data portal. Use when users ask to "show the structure" of a gene/protein.
- map_identifiers: Map between gene symbols, UniProt accessions, and Ensembl gene IDs. Use when you need to bridge identifiers across databases (e.g., gene symbol to Ensembl ID for Open Targets queries).
- get_protein_variants: Get known natural variants and mutagenesis data from UniProt for a given accession. Returns position, amino acid change, and clinical/functional description. Especially useful for interpreting MAVE variant effect scores.

After calling lookup_protein, use create_visualization with viz_type "gene_card" to display a summary card in the data portal. The gene card should include the gene name, protein function, UniProt ID, diseases, domains, and links.

LITERATURE SEARCH:
- search_literature: Search Europe PMC for biomedical papers. Use when users ask "What papers describe perturbation of gene X?" or want to find publications about a gene, disease, or experimental method. Returns titles, authors, journals, citation counts, and links. Display results as a table with create_visualization.

DRUGGABILITY (PHAROS) — ALWAYS USE FOR DRUGGABILITY QUESTIONS:
When users ask "Is gene X druggable?", "Can we target X?", "What drugs target X?", or anything about druggability, you MUST call get_druggability FIRST (before Open Targets). Pharos provides the authoritative Target Development Level classification that Open Targets does not have.

- get_druggability: Check if a gene/protein target is druggable using the Pharos database (NIH). Returns the Target Development Level:
  * Tclin = approved drug target — has at least one approved drug
  * Tchem = chemical tool target — has active compounds but no approved drug
  * Tbio = biological target — has biological evidence but no active compounds
  * Tdark = understudied target — little known about this protein
  Also returns protein family, description, disease count, and ligand count.
  After calling get_druggability, you may OPTIONALLY also query Open Targets knownDrugs for detailed drug names and mechanisms.

VARIANT MOLECULAR CONSEQUENCES (PROTVAR):
You have access to ProtVar (EBI) for deep molecular annotation of human missense variants.

- annotate_variant: Annotate a variant with full molecular consequences. Accepts many input formats: dbSNP IDs (rs1042779), gnomAD (19-1010539-G-C), VCF-like (19 1010539 G C), HGVS, or UniProt+change (P80404 Gln56Arg). Returns pathogenicity predictions (AlphaMissense, EVE, ESM-1b, Conservation), protein stability (FoldX ddG), CADD score, gnomAD allele frequency, and gene/isoform mapping. Use when users ask "What does this variant do?", "Is rs123456 pathogenic?", "What's the effect of this mutation?", or when interpreting MAVE variant results.

- get_variant_structural_context: Get detailed structural context at a specific protein residue position. Returns PDB structure mappings, binding pocket information, protein-protein interaction interfaces, FoldX stability predictions for all substitutions, functional features (domains, active sites), and co-located variants with clinical significance. Use when users ask "What's special about position 493 in Q9NUW8?" or to understand the structural impact of a variant identified by annotate_variant.

ProtVar workflow:
1. Use annotate_variant for initial variant annotation (pathogenicity + stability)
2. If the variant maps to a protein position, optionally use get_variant_structural_context for deeper structural analysis
3. Combine with get_protein_variants (UniProt) and get_protein_structure (AlphaFold) for full context

Workflow for gene queries:
1. First use our Catalogue tools (search_target_summary, query_perturbation_data) for perturbation data
2. Use lookup_protein to add protein context
3. Optionally use get_protein_structure for 3D visualization
4. Use get_druggability (Pharos) for druggability classification — ALWAYS before Open Targets for drug questions
5. Use Open Targets for disease associations and additional drug detail
6. Use search_literature for relevant publications"""

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

            tools = [types.Tool(function_declarations=_get_all_tool_declarations())]
            config = types.GenerateContentConfig(
                system_instruction=SYSTEM_INSTRUCTION,
                tools=tools,
                temperature=0.3,
            )

            # Iterative tool-calling loop
            max_iterations = 12
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
                    elif tool_name == "get_protein_structure":
                        # get_protein_structure emits a visualization AND returns data to LLM
                        yield _sse_event(
                            "tool_call",
                            {"tool": tool_name, "description": "Loading 3D protein structure..."},
                        )
                        yield _sse_event(
                            "thinking", {"status": "Loading 3D protein structure..."}
                        )
                        try:
                            result = await _tool_get_protein_structure(tool_args)
                        except Exception as exc:
                            logger.exception("Tool %s failed", tool_name)
                            result = {"error": str(exc)}

                        if "error" not in result:
                            viz_data = {
                                "type": "protein_structure",
                                "title": f"AlphaFold Structure: {result.get('gene', result.get('uniprot_id', ''))}",
                                "data": result,
                            }
                            yield _sse_event("visualization", viz_data)

                        function_response_parts.append(
                            types.Part.from_function_response(
                                name=tool_name,
                                response=result,
                            )
                        )
                    else:
                        # Choose user-friendly description for OT tools
                        desc = _OT_TOOL_DESCRIPTIONS.get(
                            tool_name, f"Calling {tool_name}..."
                        )
                        yield _sse_event(
                            "tool_call",
                            {"tool": tool_name, "description": desc},
                        )
                        yield _sse_event(
                            "thinking", {"status": desc}
                        )

                        if tool_name in _ot_tool_names:
                            try:
                                result = await _call_open_targets_tool(
                                    tool_name, tool_args
                                )
                            except Exception as exc:
                                logger.exception("OT tool %s failed", tool_name)
                                result = {"error": str(exc)}
                        elif tool_name in TOOL_HANDLERS:
                            try:
                                result = await TOOL_HANDLERS[tool_name](tool_args)
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
