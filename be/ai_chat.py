import json
import logging
import os
import uuid
from typing import Any, Dict, List, Optional
from urllib.parse import quote

import aiohttp

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import StreamingResponse
from starlette.responses import Response
from pydantic import BaseModel

from google import genai
from google.genai import types

try:
    from mcp.client.streamable_http import streamablehttp_client
    from mcp import ClientSession

    MCP_AVAILABLE = True
except ImportError:
    MCP_AVAILABLE = False

from auth import get_current_user
from data_query import db_pools

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/v1/chat", tags=["AI Chat"])

# --- Configuration ---

_config: Dict[str, Any] = {}

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


class SessionUpdateRequest(BaseModel):
    title: str


# --- Session CRUD endpoints ---


async def _verify_session_ownership(session_id: str, user_id: int):
    """Verify session exists and belongs to user. Returns row or raises 404."""
    row = await db_pools["pg"].fetchrow(
        "SELECT id, user_id, title, created_at, updated_at FROM chat_sessions WHERE id = $1",
        uuid.UUID(session_id),
    )
    if not row or row["user_id"] != user_id:
        raise HTTPException(status_code=404, detail="Session not found")
    return row


@router.get("/sessions")
async def list_sessions(user: dict = Depends(get_current_user)):
    rows = await db_pools["pg"].fetch(
        """SELECT s.id, s.title, s.created_at, s.updated_at,
                  COUNT(DISTINCT m.id) AS message_count,
                  COUNT(DISTINCT v.id) AS viz_count
           FROM chat_sessions s
           LEFT JOIN chat_messages m ON m.session_id = s.id
           LEFT JOIN chat_visualizations v ON v.session_id = s.id
           WHERE s.user_id = $1
           GROUP BY s.id
           ORDER BY s.updated_at DESC
           LIMIT 50""",
        user["id"],
    )
    return {
        "sessions": [
            {
                "id": str(r["id"]),
                "title": r["title"],
                "created_at": r["created_at"].isoformat(),
                "updated_at": r["updated_at"].isoformat(),
                "message_count": r["message_count"],
                "viz_count": r["viz_count"],
            }
            for r in rows
        ]
    }


@router.post("/sessions")
async def create_session(user: dict = Depends(get_current_user)):
    session_id = uuid.uuid4()
    await db_pools["pg"].execute(
        "INSERT INTO chat_sessions (id, user_id) VALUES ($1, $2)",
        session_id,
        user["id"],
    )
    row = await db_pools["pg"].fetchrow(
        "SELECT id, title, created_at, updated_at FROM chat_sessions WHERE id = $1",
        session_id,
    )
    return {
        "id": str(row["id"]),
        "title": row["title"],
        "created_at": row["created_at"].isoformat(),
        "updated_at": row["updated_at"].isoformat(),
    }


@router.get("/sessions/{session_id}")
async def get_session(session_id: str, user: dict = Depends(get_current_user)):
    row = await _verify_session_ownership(session_id, user["id"])
    messages = await db_pools["pg"].fetch(
        "SELECT role, content, created_at FROM chat_messages WHERE session_id = $1 ORDER BY id",
        uuid.UUID(session_id),
    )
    visualizations = await db_pools["pg"].fetch(
        "SELECT id, viz_type, title, data, created_at FROM chat_visualizations WHERE session_id = $1 ORDER BY id",
        uuid.UUID(session_id),
    )
    return {
        "id": str(row["id"]),
        "title": row["title"],
        "created_at": row["created_at"].isoformat(),
        "updated_at": row["updated_at"].isoformat(),
        "messages": [
            {
                "role": m["role"],
                "content": m["content"],
                "created_at": m["created_at"].isoformat(),
            }
            for m in messages
        ],
        "visualizations": [
            {
                "id": v["id"],
                "viz_type": v["viz_type"],
                "title": v["title"],
                "data": json.loads(v["data"]) if isinstance(v["data"], str) else v["data"],
            }
            for v in visualizations
        ],
    }


@router.patch("/sessions/{session_id}")
async def update_session(
    session_id: str,
    body: SessionUpdateRequest,
    user: dict = Depends(get_current_user),
):
    await _verify_session_ownership(session_id, user["id"])
    await db_pools["pg"].execute(
        "UPDATE chat_sessions SET title = $1, updated_at = NOW() WHERE id = $2",
        body.title[:200],
        uuid.UUID(session_id),
    )
    row = await db_pools["pg"].fetchrow(
        "SELECT id, title, created_at, updated_at FROM chat_sessions WHERE id = $1",
        uuid.UUID(session_id),
    )
    return {
        "id": str(row["id"]),
        "title": row["title"],
        "created_at": row["created_at"].isoformat(),
        "updated_at": row["updated_at"].isoformat(),
    }


@router.delete("/sessions/{session_id}")
async def delete_session(session_id: str, user: dict = Depends(get_current_user)):
    await _verify_session_ownership(session_id, user["id"])
    await db_pools["pg"].execute(
        "DELETE FROM chat_sessions WHERE id = $1", uuid.UUID(session_id)
    )
    return Response(status_code=204)


@router.delete("/sessions/{session_id}/visualizations/{viz_id}")
async def delete_visualization(
    session_id: str, viz_id: int, user: dict = Depends(get_current_user)
):
    await _verify_session_ownership(session_id, user["id"])
    await db_pools["pg"].execute(
        "DELETE FROM chat_visualizations WHERE id = $1 AND session_id = $2",
        viz_id,
        uuid.UUID(session_id),
    )
    return Response(status_code=204)


@router.delete("/sessions/{session_id}/visualizations")
async def delete_all_visualizations(
    session_id: str, user: dict = Depends(get_current_user)
):
    await _verify_session_ownership(session_id, user["id"])
    await db_pools["pg"].execute(
        "DELETE FROM chat_visualizations WHERE session_id = $1",
        uuid.UUID(session_id),
    )
    return Response(status_code=204)


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

FIND_DATASETS_FOR_TARGET_DECLARATION = types.FunctionDeclaration(
    name="find_datasets_for_target",
    description=(
        "Find all datasets that contain perturbation data for a specific gene target. "
        "Returns real dataset IDs with titles and metadata across all modalities "
        "(Perturb-seq, CRISPR screen, MAVE). "
        "ALWAYS use this tool when users ask what data is available for a gene."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "gene_name": types.Schema(
                type="STRING",
                description="Gene symbol to look up (e.g., 'TP53', 'BRCA2', 'KRAS')",
            ),
            "modality": types.Schema(
                type="STRING",
                description="Optional: filter to a specific modality ('perturb-seq', 'crispr-screen', or 'mave')",
            ),
        },
        required=["gene_name"],
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
                description="Maximum rows to return (default 20, max 500). Use higher limits (e.g. 500) for volcano plots.",
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
        "Use gene_card to display a summary card for a gene with protein info, diseases, domains, and links. "
        "For volcano plots, use the dedicated create_volcano_plot tool instead."
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

CREATE_VOLCANO_PLOT_DECLARATION = types.FunctionDeclaration(
    name="create_volcano_plot",
    description=(
        "Create a volcano plot for Perturb-seq differential expression data. "
        "The backend queries the database and renders the plot directly — much faster than using create_visualization. "
        "Use this whenever a user asks about differential expression effects of a perturbation, "
        "or explicitly requests a volcano plot."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "gene_name": types.Schema(
                type="STRING",
                description="Perturbed gene name to show DEA results for (e.g. 'BRCA2', 'TP53')",
            ),
            "dataset_id": types.Schema(
                type="STRING",
                description="Optional dataset ID to filter by. If omitted, uses data from all datasets for this gene.",
            ),
            "fc_threshold": types.Schema(
                type="NUMBER",
                description="log2 fold-change threshold for significance coloring (default 0, meaning padj-only). Set to e.g. 1.0 to require both statistical and fold-change significance.",
            ),
            "padj_threshold": types.Schema(
                type="NUMBER",
                description="Adjusted p-value threshold for significance coloring (default 0.05)",
            ),
        },
        required=["gene_name"],
    ),
)

CREATE_MAVE_HEATMAP_DECLARATION = types.FunctionDeclaration(
    name="create_mave_heatmap",
    description=(
        "Create an interactive heatmap of MAVE (Multiplexed Assay of Variant Effect) data. "
        "Shows functional scores as a position × amino acid matrix: columns are protein positions, "
        "rows are amino acid substitutions, and cell color indicates the effect score. "
        "The backend queries the database and renders the heatmap directly. "
        "Use this whenever a user asks about MAVE variant effects, functional scores, "
        "or deep mutational scanning results for a gene."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "gene_name": types.Schema(
                type="STRING",
                description="Gene name to show MAVE data for (e.g. 'BRCA1', 'TP53')",
            ),
            "dataset_id": types.Schema(
                type="STRING",
                description="Optional dataset ID. If omitted and multiple datasets exist, returns a list for user selection.",
            ),
            "score_name": types.Schema(
                type="STRING",
                description="Score metric name to display (default 'score'). Use query_perturbation_data with modality='mave' to discover available score names first if unsure.",
            ),
            "position_start": types.Schema(
                type="INTEGER",
                description="Start position for the heatmap window (default: first available position). Use to zoom into a region of interest.",
            ),
            "position_end": types.Schema(
                type="INTEGER",
                description="End position for the heatmap window (default: position_start + 29, i.e. 30 positions). Use to zoom into a region of interest.",
            ),
        },
        required=["gene_name"],
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

PREDICT_VARIANT_CONSEQUENCE_DECLARATION = types.FunctionDeclaration(
    name="predict_variant_consequence",
    description=(
        "Predict the functional consequence of a genetic variant using Ensembl VEP (Variant Effect Predictor). "
        "Returns consequence type (missense, frameshift, splice, regulatory, etc.), impact severity, "
        "affected gene/transcript, protein change, and in-silico predictions: SIFT, PolyPhen, CADD, "
        "SpliceAI (splicing impact), AlphaMissense (missense pathogenicity), LOFTEE (loss-of-function), "
        "and conservation scores. Also reports regulatory feature consequences (enhancer/promoter disruption) "
        "and colocated known variants with population frequencies. "
        "Handles ALL variant types (SNPs, indels, frameshifts) — not limited to missense like ProtVar. "
        "Use for: 'What is the consequence of variant X?', 'Predict effect of rs123', "
        "'Is this a splice variant?', 'Does this variant affect regulatory regions?', "
        "'VEP annotation for variant X'. "
        "For protein stability (FoldX) and structural context (pockets, interfaces), use ProtVar tools instead."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "variant": types.Schema(
                type="STRING",
                description=(
                    "Variant identifier in any format: "
                    "rsID (e.g. 'rs56116432'), "
                    "HGVS (e.g. 'ENST00000366667:c.803C>T', '9:g.22125504G>C'), "
                    "or VCF-style region (e.g. '9:22125504:G:C')"
                ),
            ),
        },
        required=["variant"],
    ),
)

BATCH_VARIANT_CONSEQUENCES_DECLARATION = types.FunctionDeclaration(
    name="batch_variant_consequences",
    description=(
        "Predict functional consequences for a batch of variants (up to 200) using Ensembl VEP. "
        "Returns a summary per variant: most severe consequence, impact, gene, protein change, "
        "and key scores (SIFT, PolyPhen, CADD, AlphaMissense). Efficient for annotating variant lists "
        "from MAVE, CRISPR, or other screening results. "
        "Use for: 'Annotate these variants', 'VEP for this list of variants', "
        "'What are the consequences of these SNPs?'"
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "variants": types.Schema(
                type="ARRAY",
                items=types.Schema(
                    type="STRING",
                    description="Variant identifier (rsID, HGVS, or VCF-style region)",
                ),
                description=(
                    "List of variant identifiers (max 200). "
                    "All formats accepted: rsIDs, HGVS, VCF-style regions. Can be mixed."
                ),
            ),
        },
        required=["variants"],
    ),
)

CREATE_GENE_NETWORK_DECLARATION = types.FunctionDeclaration(
    name="create_gene_interaction_network",
    description=(
        "Create an interactive gene interaction network from Perturb-seq differential expression data. "
        "Shows the perturbed gene in the center with differentially expressed genes radiating outward. "
        "Edge width reflects the magnitude of the effect (|log2FC|) and nodes are colored by direction "
        "(up/down-regulated). The backend queries the database and renders the network directly. "
        "Use this when a user asks for a network view of perturbation effects, gene interactions, "
        "or wants to visualize which genes are affected by a perturbation."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "gene_name": types.Schema(
                type="STRING",
                description="Perturbed gene name to show interaction network for (e.g. 'BRCA2', 'TP53')",
            ),
            "dataset_id": types.Schema(
                type="STRING",
                description="Optional dataset ID to filter by. If omitted and multiple datasets exist, returns a list for user selection.",
            ),
            "padj_threshold": types.Schema(
                type="NUMBER",
                description="Adjusted p-value threshold for including genes in the network (default 0.05). Only genes with padj below this are shown.",
            ),
            "fc_threshold": types.Schema(
                type="NUMBER",
                description="Minimum absolute log2 fold-change threshold for including genes (default 0.5). Filters out small effects.",
            ),
            "max_genes": types.Schema(
                type="INTEGER",
                description="Maximum number of differentially expressed genes to show (default 30, max 50). Top genes selected by significance.",
            ),
        },
        required=["gene_name"],
    ),
)

GET_STRING_INTERACTIONS_DECLARATION = types.FunctionDeclaration(
    name="get_string_interactions",
    description=(
        "Fetch protein-protein interaction partners for a gene from the STRING database. "
        "Returns interaction partners with combined confidence scores and evidence channel "
        "breakdown (coexpression, experimental, database, textmining). "
        "Renders an interactive network in the data portal: query gene at center, "
        "partners radiating outward, edge width encoding confidence. "
        "Use when users ask 'What proteins interact with X?', 'Show STRING network for X', "
        "'What are the known interaction partners of X?', or to provide protein interaction "
        "context alongside perturbation data."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "gene_name": types.Schema(
                type="STRING",
                description="Human gene symbol (e.g. 'TP53', 'BRCA2', 'KRAS')",
            ),
            "required_score": types.Schema(
                type="INTEGER",
                description=(
                    "Minimum combined confidence score (0-1000). "
                    "400=medium confidence (default), 700=high, 900=highest confidence"
                ),
            ),
            "limit": types.Schema(
                type="INTEGER",
                description="Maximum number of interaction partners to return (default 20, max 50)",
            ),
            "network_type": types.Schema(
                type="STRING",
                description=(
                    "Network type: 'functional' (default, all evidence channels) "
                    "or 'physical' (physical binding evidence only)"
                ),
            ),
        },
        required=["gene_name"],
    ),
)

GET_FUNCTIONAL_ENRICHMENT_DECLARATION = types.FunctionDeclaration(
    name="get_functional_enrichment",
    description=(
        "Run functional enrichment analysis on a set of genes using STRING. "
        "Returns significantly enriched GO terms (Biological Process, Molecular Function, "
        "Cellular Component), KEGG pathways, Reactome pathways, and Pfam domains. "
        "Returns p-values, false discovery rates, and matching gene counts. "
        "Use when users ask 'What pathways are enriched in these genes?', "
        "'What biological processes are shared by X, Y, Z?', or after identifying a set of "
        "interacting/co-regulated genes from perturbation data. "
        "Present results as a table using create_visualization."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "genes": types.Schema(
                type="ARRAY",
                items=types.Schema(type="STRING"),
                description="List of human gene symbols to analyze (e.g. ['TP53', 'BRCA2', 'MDM2'])",
            ),
        },
        required=["genes"],
    ),
)

INTERNAL_TOOL_DECLARATIONS = [
    SEARCH_DATASETS_DECLARATION,
    SEARCH_TARGET_SUMMARY_DECLARATION,
    FIND_DATASETS_FOR_TARGET_DECLARATION,
    QUERY_PERTURBATION_DATA_DECLARATION,
    GET_CATALOGUE_SUMMARY_DECLARATION,
    CREATE_VISUALIZATION_DECLARATION,
    CREATE_VOLCANO_PLOT_DECLARATION,
    CREATE_MAVE_HEATMAP_DECLARATION,
    CREATE_GENE_NETWORK_DECLARATION,
    GET_STRING_INTERACTIONS_DECLARATION,
    GET_FUNCTIONAL_ENRICHMENT_DECLARATION,
    LOOKUP_PROTEIN_DECLARATION,
    GET_PROTEIN_STRUCTURE_DECLARATION,
    MAP_IDENTIFIERS_DECLARATION,
    GET_PROTEIN_VARIANTS_DECLARATION,
    SEARCH_LITERATURE_DECLARATION,
    GET_DRUGGABILITY_DECLARATION,
    ANNOTATE_VARIANT_DECLARATION,
    GET_VARIANT_STRUCTURAL_CONTEXT_DECLARATION,
    PREDICT_VARIANT_CONSEQUENCE_DECLARATION,
    BATCH_VARIANT_CONSEQUENCES_DECLARATION,
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


async def _tool_find_datasets_for_target(args: dict) -> dict:
    """Find all datasets containing data for a gene across all modalities."""
    gene_name = args.get("gene_name", "").upper()
    modality_filter = args.get("modality")

    if not gene_name:
        return {"error": "gene_name is required"}

    pg_pool = db_pools["pg"]
    es = db_pools["es"]

    queries = {
        "perturb-seq": (
            "perturb_seq_dea",
            """SELECT dataset_id,
                      COUNT(*) AS total_genes,
                      COUNT(*) FILTER (WHERE padj < 0.05) AS significant_genes
               FROM perturb_seq_dea
               WHERE perturbed_target_symbol = $1 AND gene IS NOT NULL
               GROUP BY dataset_id
               ORDER BY significant_genes DESC""",
        ),
        "crispr-screen": (
            "crispr_data",
            """SELECT dataset_id,
                      COUNT(*) AS total_genes,
                      COUNT(*) FILTER (WHERE significant = 'True') AS significant_genes
               FROM crispr_data
               WHERE perturbed_target_symbol = $1
               GROUP BY dataset_id
               ORDER BY significant_genes DESC""",
        ),
        "mave": (
            "mave_data",
            """SELECT dataset_id,
                      COUNT(*) AS total_variants
               FROM mave_data
               WHERE perturbed_target_symbol = $1
               GROUP BY dataset_id
               ORDER BY total_variants DESC""",
        ),
    }

    if modality_filter:
        if modality_filter not in queries:
            return {"error": f"Unknown modality: {modality_filter}"}
        queries = {modality_filter: queries[modality_filter]}

    all_datasets = []
    all_dataset_ids = set()

    for modality, (table, query) in queries.items():
        rows = await pg_pool.fetch(query, gene_name)
        for row in rows:
            ds_id = row["dataset_id"]
            all_dataset_ids.add(ds_id)
            entry = {"dataset_id": ds_id, "modality": modality}
            if modality == "mave":
                entry["stats"] = {"total_variants": row["total_variants"]}
            else:
                entry["stats"] = {
                    "total_genes": row["total_genes"],
                    "significant_genes": row["significant_genes"],
                }
            all_datasets.append(entry)

    if not all_datasets:
        return {
            "total": 0,
            "datasets": [],
            "message": f"No perturbation data found for {gene_name} in any modality.",
        }

    # Fetch metadata from Elasticsearch for all dataset IDs
    if all_dataset_ids:
        es_result = await es.search(
            index=ES_DATASET_SUMMARY,
            query={"terms": {"dataset_id": list(all_dataset_ids)}},
            size=len(all_dataset_ids),
        )
        es_meta = {}
        for hit in es_result.get("hits", {}).get("hits", []):
            src = hit["_source"]
            es_meta[src.get("dataset_id")] = {
                "title": src.get("study_title", ""),
                "year": src.get("study_year", ""),
                "modality": src.get("data_modalities", ""),
                "license": src.get("license_labels", ""),
                "tissue": src.get("tissue_labels", ""),
                "cell_types": src.get("cell_type_labels", ""),
                "cell_lines": src.get("cell_line_labels", ""),
                "sex": src.get("sex_labels_labels", ""),
                "developmental_stage": src.get("developmental_stage_labels", ""),
                "disease": src.get("disease_labels", ""),
            }

        for ds in all_datasets:
            meta = es_meta.get(ds["dataset_id"], {})
            ds.update(meta)

    return {
        "gene": gene_name,
        "total": len(all_datasets),
        "datasets": all_datasets,
    }


async def _tool_query_perturbation_data(args: dict) -> dict:
    modality = args.get("modality", "crispr-screen")
    gene_name = args.get("gene_name")
    dataset_id = args.get("dataset_id")
    limit = min(args.get("limit", 20), 500)

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


async def _tool_create_volcano_plot(args: dict) -> dict:
    """Query perturb_seq_dea and return pre-formatted volcano plot data + summary.

    Two-step flow:
    1. If no dataset_id: return available datasets ranked by significance count
       so Gemini can ask the user which one to plot.
    2. If dataset_id provided: query DEA data sorted by padj ASC, build plot.
    """
    gene_name = args.get("gene_name", "")
    dataset_id = args.get("dataset_id")
    fc_threshold = float(args.get("fc_threshold", 0))
    padj_threshold = float(args.get("padj_threshold", 0.05))

    if not gene_name:
        return {"error": "gene_name is required"}

    pg_pool = db_pools["pg"]

    # Step 1: No dataset_id — list available datasets for this gene
    if not dataset_id:
        ds_query = """
            SELECT dataset_id,
                   COUNT(*) AS total_genes,
                   COUNT(*) FILTER (WHERE padj < 0.05) AS significant_genes
            FROM perturb_seq_dea
            WHERE perturbed_target_symbol = $1 AND gene IS NOT NULL
            GROUP BY dataset_id
            ORDER BY significant_genes DESC
        """
        ds_rows = await pg_pool.fetch(ds_query, gene_name.upper())

        if not ds_rows:
            return {"error": f"No Perturb-seq DEA data found for {gene_name.upper()}"}

        if len(ds_rows) == 1:
            # Only one dataset — proceed directly
            dataset_id = ds_rows[0]["dataset_id"]
        else:
            # Multiple datasets — return list for Gemini to present
            datasets = [
                {
                    "dataset_id": row["dataset_id"],
                    "total_genes": row["total_genes"],
                    "significant_genes": row["significant_genes"],
                }
                for row in ds_rows
            ]
            return {
                "action": "choose_dataset",
                "perturbed_gene": gene_name.upper(),
                "available_datasets": datasets,
                "message": (
                    f"{gene_name.upper()} has Perturb-seq DEA data in {len(datasets)} datasets. "
                    "Ask the user which dataset they want to see the volcano plot for."
                ),
            }

    # Step 2: Query DEA data for specific dataset, sorted by padj
    query = """
        SELECT gene, log2foldchange, padj
        FROM perturb_seq_dea
        WHERE perturbed_target_symbol = $1 AND dataset_id = $2 AND gene IS NOT NULL
        ORDER BY padj ASC
        LIMIT 500
    """
    rows = await pg_pool.fetch(query, gene_name.upper(), dataset_id)

    if not rows:
        return {"error": f"No Perturb-seq DEA data found for {gene_name.upper()} in dataset {dataset_id}"}

    genes = []
    log2fc = []
    padj_vals = []
    n_up = 0
    n_down = 0
    n_ns = 0

    for row in rows:
        g = row["gene"]
        fc = row["log2foldchange"]
        pv = row["padj"]
        if g is None or fc is None or pv is None:
            continue
        fc = float(fc)
        pv = float(pv)
        genes.append(str(g))
        log2fc.append(fc)
        padj_vals.append(pv)

        if pv < padj_threshold and fc > fc_threshold:
            n_up += 1
        elif pv < padj_threshold and fc < -fc_threshold:
            n_down += 1
        else:
            n_ns += 1

    if not genes:
        return {"error": f"No valid DEA data points for {gene_name.upper()} in dataset {dataset_id}"}

    # Build the visualization payload (sent directly to frontend via SSE)
    viz_data = {
        "type": "volcano_plot",
        "title": f"Volcano Plot: {gene_name.upper()} Perturbation",
        "data": {
            "genes": genes,
            "log2fc": log2fc,
            "padj": padj_vals,
            "fc_threshold": fc_threshold,
            "padj_threshold": padj_threshold,
            "perturbed_gene": gene_name.upper(),
        },
    }

    # Summary returned to Gemini so it can write an informative text response
    summary = {
        "status": "success",
        "perturbed_gene": gene_name.upper(),
        "dataset_id": dataset_id,
        "total_genes": len(genes),
        "significantly_upregulated": n_up,
        "significantly_downregulated": n_down,
        "not_significant": n_ns,
        "fc_threshold": fc_threshold,
        "padj_threshold": padj_threshold,
    }

    return {"viz_data": viz_data, "summary": summary}


async def _tool_create_mave_heatmap(args: dict) -> dict:
    """Query mave_data and return pre-formatted heatmap data + summary.

    Two-step flow (same as volcano plot):
    1. If no dataset_id: return available datasets ranked by variant count
       so Gemini can ask the user which one to plot.
    2. If dataset_id provided: query MAVE data, build position × amino acid
       matrix, and return viz payload.
    """
    gene_name = args.get("gene_name", "")
    dataset_id = args.get("dataset_id")
    score_name = args.get("score_name", "score")
    position_start = args.get("position_start")
    position_end = args.get("position_end")

    if not gene_name:
        return {"error": "gene_name is required"}

    pg_pool = db_pools["pg"]

    # Step 1: No dataset_id — list available datasets for this gene
    if not dataset_id:
        ds_query = """
            SELECT dataset_id,
                   COUNT(*) AS total_variants,
                   COUNT(DISTINCT perturbation_position) AS positions_covered
            FROM mave_data
            WHERE perturbed_target_symbol = $1 AND score_name = $2
                  AND perturbation_position IS NOT NULL
            GROUP BY dataset_id
            ORDER BY total_variants DESC
        """
        ds_rows = await pg_pool.fetch(ds_query, gene_name.upper(), score_name)

        if not ds_rows:
            # Discover what score names are available for this gene
            score_query = """
                SELECT DISTINCT score_name, COUNT(*) AS n
                FROM mave_data
                WHERE perturbed_target_symbol = $1 AND perturbation_position IS NOT NULL
                GROUP BY score_name ORDER BY n DESC
            """
            score_rows = await pg_pool.fetch(score_query, gene_name.upper())
            available = [r["score_name"] for r in score_rows]
            error_msg = f"No MAVE data found for {gene_name.upper()} with score_name='{score_name}'"
            if available:
                error_msg += f". Available score names: {', '.join(available)}"
            return {"error": error_msg}

        if len(ds_rows) == 1:
            dataset_id = ds_rows[0]["dataset_id"]
        else:
            datasets = [
                {
                    "dataset_id": row["dataset_id"],
                    "total_variants": row["total_variants"],
                    "positions_covered": row["positions_covered"],
                }
                for row in ds_rows
            ]
            return {
                "action": "choose_dataset",
                "gene_name": gene_name.upper(),
                "available_datasets": datasets,
                "message": (
                    f"{gene_name.upper()} has MAVE data in {len(datasets)} datasets. "
                    "Ask the user which dataset they want to see the heatmap for."
                ),
            }

    # Step 2: Query MAVE data for specific dataset
    # First, determine position range if not provided
    if position_start is None:
        range_query = """
            SELECT MIN(perturbation_position) AS min_pos,
                   MAX(perturbation_position) AS max_pos
            FROM mave_data
            WHERE perturbed_target_symbol = $1 AND dataset_id = $2 AND score_name = $3
                  AND perturbation_position IS NOT NULL
        """
        range_row = await pg_pool.fetchrow(
            range_query, gene_name.upper(), dataset_id, score_name
        )
        if not range_row or range_row["min_pos"] is None:
            return {
                "error": f"No MAVE data found for {gene_name.upper()} in dataset {dataset_id}"
            }
        position_start = range_row["min_pos"]
        total_positions = range_row["max_pos"] - range_row["min_pos"] + 1
        # Default window: 30 positions
        if position_end is None:
            position_end = position_start + min(29, total_positions - 1)
    elif position_end is None:
        position_end = position_start + 29

    query = """
        SELECT perturbation_position, perturbation_aa_wt, perturbation_aa_change, score_value
        FROM mave_data
        WHERE perturbed_target_symbol = $1 AND dataset_id = $2 AND score_name = $3
              AND perturbation_position >= $4 AND perturbation_position <= $5
              AND perturbation_position IS NOT NULL
        ORDER BY perturbation_position, perturbation_aa_change
    """
    rows = await pg_pool.fetch(
        query, gene_name.upper(), dataset_id, score_name, position_start, position_end
    )

    if not rows:
        return {
            "error": f"No MAVE data found for {gene_name.upper()} in dataset {dataset_id} "
            f"at positions {position_start}-{position_end}"
        }

    # Build position × amino acid matrix
    # Collect all data points and unique amino acids
    positions_data = {}  # {position: {aa: score}}
    wt_residues = {}  # {position: wt_aa}
    all_aas = set()

    for row in rows:
        pos = row["perturbation_position"]
        aa_wt = row["perturbation_aa_wt"]
        aa_change = row["perturbation_aa_change"]
        score = row["score_value"]

        if pos is None or score is None:
            continue

        is_ref = aa_change == "="
        aa = aa_wt if is_ref else aa_change

        if aa is None:
            continue

        if pos not in positions_data:
            positions_data[pos] = {}
        positions_data[pos][aa] = float(score)
        all_aas.add(aa)

        if aa_wt:
            wt_residues[pos] = aa_wt

    if not positions_data:
        return {
            "error": f"No valid MAVE data points for {gene_name.upper()} in dataset {dataset_id}"
        }

    # Sort amino acids and positions
    sorted_aas = sorted(all_aas)
    sorted_positions = sorted(positions_data.keys())

    # Build 2D matrix: rows = amino acids, columns = positions
    # Each cell is the score value (or null if not measured)
    z_matrix = []
    wt_annotations = []  # [{row, col}] for WT residue markers

    for aa_idx, aa in enumerate(sorted_aas):
        row_values = []
        for pos_idx, pos in enumerate(sorted_positions):
            score = positions_data[pos].get(aa)
            row_values.append(score)
            # Mark WT residue positions (only if cell has a score)
            if wt_residues.get(pos) == aa and score is not None:
                wt_annotations.append({"row": aa_idx, "col": pos_idx})
        z_matrix.append(row_values)

    position_labels = [str(p) for p in sorted_positions]

    viz_data = {
        "type": "mave_heatmap",
        "title": f"MAVE Heatmap: {gene_name.upper()} (positions {position_start}-{position_end})",
        "data": {
            "z": z_matrix,
            "amino_acids": sorted_aas,
            "positions": position_labels,
            "wt_annotations": wt_annotations,
            "gene_name": gene_name.upper(),
            "position_start": position_start,
            "position_end": position_end,
        },
    }

    # Summary for Gemini
    total_variants = sum(
        1 for row in z_matrix for v in row if v is not None
    )
    n_positions = len(sorted_positions)
    n_aas = len(sorted_aas)

    summary = {
        "status": "success",
        "gene_name": gene_name.upper(),
        "dataset_id": dataset_id,
        "score_name": score_name,
        "position_range": f"{position_start}-{position_end}",
        "positions_with_data": n_positions,
        "amino_acids_observed": n_aas,
        "total_variants_shown": total_variants,
        "message": (
            f"Heatmap displayed for {gene_name.upper()} showing {total_variants} variant scores "
            f"across {n_positions} positions ({position_start}-{position_end}) "
            f"and {n_aas} amino acid substitutions."
        ),
    }

    return {"viz_data": viz_data, "summary": summary}


async def _tool_create_gene_interaction_network(args: dict) -> dict:
    """Query perturb_seq_dea and return a gene interaction network payload.

    Same two-step dataset selection flow as volcano plot and MAVE heatmap.
    """
    gene_name = args.get("gene_name", "")
    dataset_id = args.get("dataset_id")
    padj_threshold = args.get("padj_threshold", 0.05)
    fc_threshold = args.get("fc_threshold", 0.5)
    max_genes = args.get("max_genes", 30)

    if not gene_name:
        return {"error": "gene_name is required"}

    # Clamp max_genes
    if max_genes is None or max_genes < 1:
        max_genes = 30
    elif max_genes > 50:
        max_genes = 50

    pg_pool = db_pools["pg"]

    # Step 1: Determine dataset_id (same pattern as volcano plot)
    if not dataset_id:
        ds_query = """
            SELECT dataset_id,
                   COUNT(*) AS total_genes,
                   COUNT(*) FILTER (WHERE padj < 0.05) AS significant_genes
            FROM perturb_seq_dea
            WHERE perturbed_target_symbol = $1
            GROUP BY dataset_id
            ORDER BY significant_genes DESC
        """
        ds_rows = await pg_pool.fetch(ds_query, gene_name.upper())

        if not ds_rows:
            return {
                "error": f"No Perturb-seq DEA data found for {gene_name.upper()}"
            }

        if len(ds_rows) == 1:
            dataset_id = ds_rows[0]["dataset_id"]
        else:
            datasets = [
                {
                    "dataset_id": row["dataset_id"],
                    "total_genes": row["total_genes"],
                    "significant_genes": row["significant_genes"],
                }
                for row in ds_rows
            ]
            return {
                "action": "choose_dataset",
                "perturbed_gene": gene_name.upper(),
                "available_datasets": datasets,
                "message": (
                    f"{gene_name.upper()} has Perturb-seq DEA data in {len(datasets)} datasets. "
                    "Ask the user which dataset they want to see the interaction network for."
                ),
            }

    # Step 2: Query significant DEGs sorted by significance
    query = """
        SELECT gene, log2foldchange, padj
        FROM perturb_seq_dea
        WHERE perturbed_target_symbol = $1 AND dataset_id = $2
              AND gene IS NOT NULL AND padj IS NOT NULL AND log2foldchange IS NOT NULL
              AND padj < $3 AND ABS(log2foldchange) >= $4
        ORDER BY padj ASC
        LIMIT $5
    """
    rows = await pg_pool.fetch(
        query,
        gene_name.upper(),
        dataset_id,
        padj_threshold,
        fc_threshold,
        max_genes,
    )

    if not rows:
        return {
            "error": (
                f"No significant genes found for {gene_name.upper()} in dataset {dataset_id} "
                f"with padj < {padj_threshold} and |log2FC| >= {fc_threshold}. "
                "Try relaxing the thresholds."
            )
        }

    # Build network data: center node + DEG nodes + edges
    nodes = []
    edges = []
    n_up = 0
    n_down = 0

    # Center node: the perturbed gene
    nodes.append(
        {
            "id": gene_name.upper(),
            "label": gene_name.upper(),
            "type": "perturbed",
            "log2fc": 0,
            "padj": 0,
        }
    )

    for row in rows:
        g = str(row["gene"])
        fc = float(row["log2foldchange"])
        pv = float(row["padj"])

        direction = "up" if fc > 0 else "down"
        if direction == "up":
            n_up += 1
        else:
            n_down += 1

        nodes.append(
            {
                "id": g,
                "label": g,
                "type": direction,
                "log2fc": round(fc, 4),
                "padj": pv,
            }
        )

        edges.append(
            {
                "source": gene_name.upper(),
                "target": g,
                "weight": round(abs(fc), 4),
            }
        )

    viz_data = {
        "type": "gene_interaction_network",
        "title": f"Gene Interaction Network: {gene_name.upper()} Perturbation",
        "data": {
            "nodes": nodes,
            "edges": edges,
            "perturbed_gene": gene_name.upper(),
            "fc_threshold": fc_threshold,
            "padj_threshold": padj_threshold,
        },
    }

    summary = {
        "status": "success",
        "perturbed_gene": gene_name.upper(),
        "dataset_id": dataset_id,
        "total_deg_shown": len(rows),
        "upregulated": n_up,
        "downregulated": n_down,
        "fc_threshold": fc_threshold,
        "padj_threshold": padj_threshold,
        "top_upregulated": [
            {"gene": str(r["gene"]), "log2fc": round(float(r["log2foldchange"]), 3)}
            for r in rows
            if float(r["log2foldchange"]) > 0
        ][:5],
        "top_downregulated": [
            {"gene": str(r["gene"]), "log2fc": round(float(r["log2foldchange"]), 3)}
            for r in rows
            if float(r["log2foldchange"]) < 0
        ][:5],
    }

    return {"viz_data": viz_data, "summary": summary}


STRING_BASE = "https://string-db.org/api/json"
STRING_CALLER = "perturbation-catalogue"


async def _tool_get_string_interactions(args: dict) -> dict:
    """Fetch STRING protein interaction partners and build a network visualization."""
    gene_name = args.get("gene_name", "").strip()
    if not gene_name:
        return {"error": "gene_name is required"}

    required_score = max(0, min(1000, int(args.get("required_score", 400))))
    limit = min(int(args.get("limit", 20)), 50)
    network_type = args.get("network_type", "functional")

    url = f"{STRING_BASE}/interaction_partners"
    params = {
        "identifiers": gene_name,
        "species": 9606,
        "limit": limit,
        "required_score": required_score,
        "network_type": network_type,
        "caller_identity": STRING_CALLER,
    }

    try:
        async with aiohttp.ClientSession() as http:
            async with http.get(
                url, params=params, timeout=aiohttp.ClientTimeout(total=20)
            ) as resp:
                if resp.status == 400:
                    body = await resp.text()
                    return {"error": f"STRING API error: {body[:300]}"}
                if resp.status != 200:
                    return {"error": f"STRING API returned status {resp.status}"}
                # STRING returns Content-Type: text/json, not application/json
                data = await resp.json(content_type=None)
    except Exception as exc:
        return {"error": f"STRING API error: {str(exc)}"}

    if not data:
        return {
            "error": (
                f"No STRING interactions found for '{gene_name}' "
                f"with required_score >= {required_score}. "
                "Try lowering required_score to 400 (medium confidence)."
            )
        }

    # Build network: query gene at center, partners radiating outward
    query_gene = data[0].get("preferredName_A", gene_name)

    nodes = [
        {
            "id": query_gene,
            "label": query_gene,
            "type": "query",
            "score": 1.0,
            "degree": len(data),
        }
    ]
    edges = []
    seen = set()

    for row in data:
        partner = row.get("preferredName_B", "")
        if not partner or partner in seen:
            continue
        seen.add(partner)

        score = float(row.get("score", 0))

        evidence = {}
        for key, label in [
            ("escore", "experimental"),
            ("dscore", "database"),
            ("ascore", "coexpression"),
            ("tscore", "textmining"),
            ("nscore", "neighborhood"),
            ("fscore", "fusion"),
            ("pscore", "phylogenetic"),
        ]:
            val = float(row.get(key, 0))
            if val > 0:
                evidence[label] = round(val, 3)

        nodes.append(
            {
                "id": partner,
                "label": partner,
                "type": "partner",
                "score": round(score, 4),
                "evidence": evidence,
            }
        )
        edges.append(
            {
                "source": query_gene,
                "target": partner,
                "weight": round(score, 4),
            }
        )

    if not edges:
        return {
            "error": (
                f"No interaction partners returned for '{gene_name}'. "
                "The gene symbol may not be recognized by STRING."
            )
        }

    viz_data = {
        "type": "string_interaction_network",
        "title": f"STRING Interaction Network: {query_gene}",
        "data": {
            "nodes": nodes,
            "edges": edges,
            "query_gene": query_gene,
            "required_score": required_score,
            "network_type": network_type,
        },
    }

    top_partners = sorted(edges, key=lambda e: e["weight"], reverse=True)[:5]
    summary = {
        "status": "success",
        "query_gene": query_gene,
        "total_interactions": len(edges),
        "required_score": required_score,
        "network_type": network_type,
        "top_partners": [
            {"gene": e["target"], "confidence": e["weight"]} for e in top_partners
        ],
        "string_url": f"https://string-db.org/network/{query_gene}",
        "message": (
            f"STRING network rendered for {query_gene} showing {len(edges)} interaction "
            f"partners with confidence >= {required_score / 1000:.1f}."
        ),
    }

    return {"viz_data": viz_data, "summary": summary}


async def _tool_get_functional_enrichment(args: dict) -> dict:
    """Run STRING functional enrichment on a gene set."""
    genes = args.get("genes", [])
    if not genes:
        return {"error": "genes list is required"}

    genes = [str(g).strip() for g in genes if g][:50]
    if not genes:
        return {"error": "genes list is empty after filtering"}

    url = f"{STRING_BASE}/enrichment"
    params = {
        "identifiers": "\r".join(genes),
        "species": 9606,
        "caller_identity": STRING_CALLER,
    }

    try:
        async with aiohttp.ClientSession() as http:
            async with http.get(
                url, params=params, timeout=aiohttp.ClientTimeout(total=30)
            ) as resp:
                if resp.status != 200:
                    body = await resp.text()
                    return {
                        "error": f"STRING enrichment API returned status {resp.status}: {body[:300]}"
                    }
                # STRING returns Content-Type: text/json, not application/json
                data = await resp.json(content_type=None)
    except Exception as exc:
        return {"error": f"STRING enrichment API error: {str(exc)}"}

    if not data:
        gene_preview = ", ".join(genes[:5]) + ("..." if len(genes) > 5 else "")
        return {"error": f"No significant enrichment found for genes ({gene_preview})."}

    # Group by category, sort by FDR, keep top 10 per category
    categories = {}
    for item in data:
        cat = item.get("category", "Other")
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(item)

    top_terms = []
    for cat, terms in categories.items():
        terms_sorted = sorted(
            terms, key=lambda t: t.get("fdr") if t.get("fdr") is not None else 1.0
        )[:10]
        for t in terms_sorted:
            fdr = t.get("fdr")
            pval = t.get("p_value")
            matching = t.get("preferredNames", [])
            if isinstance(matching, str):
                matching = [m.strip() for m in matching.split(",") if m.strip()]
            top_terms.append(
                {
                    "category": cat,
                    "term_id": t.get("term", ""),
                    "description": t.get("description", ""),
                    "gene_count": t.get("number_of_genes", 0),
                    "fdr": round(fdr, 6) if fdr is not None else None,
                    "p_value": round(pval, 6) if pval is not None else None,
                    "matching_genes": ", ".join(matching[:8]),
                }
            )

    top_terms.sort(key=lambda t: t["fdr"] if t["fdr"] is not None else 1.0)

    return {
        "gene_set": genes,
        "total_enriched_terms": len(data),
        "categories_found": list(categories.keys()),
        "top_terms": top_terms[:40],
        "string_url": "https://string-db.org/cgi/network?identifiers=" + "%0d".join(genes),
    }


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

# --- Ensembl VEP ---

VEP_BASE = "https://rest.ensembl.org/vep/human"
VEP_PARAMS = {
    "CADD": "1",
    "SpliceAI": "1",
    "AlphaMissense": "1",
    "LoF": "1",
    "Conservation": "1",
    "domains": "1",
    "hgvs": "1",
    "canonical": "1",
    "uniprot": "1",
    "protein": "1",
}


def _detect_variant_format(variant: str) -> str:
    """Detect variant input format: 'rsid', 'hgvs', or 'region'.

    - rsID: starts with 'rs' followed by digits (e.g. rs56116432)
    - Region/VCF-style: chr:pos:ref:alt or chr:start-end:strand/allele
    - HGVS: everything else (transcript:c.X, genomic g.X, protein p.X)
    """
    v = variant.strip()
    if v.lower().startswith("rs") and v[2:].isdigit():
        return "rsid"
    # VCF-style region: e.g. "9:22125504:G:C" or "9:22125503-22125504:1/C"
    parts = v.split(":")
    if len(parts) >= 3 and parts[0].replace("chr", "").replace("X", "1").replace("Y", "1").replace("MT", "1").replace("M", "1").isdigit():
        return "region"
    return "hgvs"


def _extract_canonical_consequence(data: list) -> dict:
    """Extract the most relevant consequence from a VEP response.

    Prioritises the canonical transcript, then falls back to the transcript
    with the most severe impact.
    """
    if not data:
        return {}

    variant = data[0]
    result = {
        "input": variant.get("input", variant.get("id", "")),
        "allele_string": variant.get("allele_string", ""),
        "location": f"{variant.get('seq_region_name', '')}:{variant.get('start', '')}-{variant.get('end', '')}",
        "most_severe_consequence": variant.get("most_severe_consequence", ""),
        "strand": variant.get("strand"),
    }

    # Pick canonical transcript, or most severe
    tx_cons = variant.get("transcript_consequences", [])
    chosen = None
    for tc in tx_cons:
        if tc.get("canonical") == 1:
            chosen = tc
            break
    if not chosen and tx_cons:
        impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
        tx_cons_sorted = sorted(tx_cons, key=lambda t: impact_order.get(t.get("impact", "MODIFIER"), 3))
        chosen = tx_cons_sorted[0]

    if chosen:
        result["gene_symbol"] = chosen.get("gene_symbol", "")
        result["gene_id"] = chosen.get("gene_id", "")
        result["transcript_id"] = chosen.get("transcript_id", "")
        result["biotype"] = chosen.get("biotype", "")
        result["impact"] = chosen.get("impact", "")
        result["consequence_terms"] = chosen.get("consequence_terms", [])
        result["amino_acids"] = chosen.get("amino_acids")
        result["codons"] = chosen.get("codons")
        result["protein_start"] = chosen.get("protein_start")
        result["hgvsc"] = chosen.get("hgvsc")
        result["hgvsp"] = chosen.get("hgvsp")

        # In-silico predictions
        predictions = {}
        if chosen.get("sift_prediction"):
            predictions["sift"] = f"{chosen['sift_prediction']}({chosen.get('sift_score', '')})"
        if chosen.get("polyphen_prediction"):
            predictions["polyphen"] = f"{chosen['polyphen_prediction']}({chosen.get('polyphen_score', '')})"
        if chosen.get("cadd_phred") is not None:
            predictions["cadd_phred"] = chosen["cadd_phred"]
        if chosen.get("cadd_raw") is not None:
            predictions["cadd_raw"] = chosen["cadd_raw"]
        # AlphaMissense
        if chosen.get("am_class"):
            predictions["alphamissense"] = f"{chosen['am_class']}({chosen.get('am_pathogenicity', '')})"
        # SpliceAI
        splice_keys = [k for k in chosen if k.startswith("spliceai_pred_ds_")]
        if splice_keys:
            splice_scores = {k.replace("spliceai_pred_", ""): chosen[k] for k in sorted(chosen.keys()) if k.startswith("spliceai_pred_")}
            predictions["spliceai"] = splice_scores
        # LOFTEE
        if chosen.get("lof"):
            predictions["loftee"] = chosen["lof"]
            if chosen.get("lof_filter"):
                predictions["loftee_filter"] = chosen["lof_filter"]
        # Conservation
        if chosen.get("conservation") is not None:
            predictions["conservation"] = chosen["conservation"]

        if predictions:
            result["predictions"] = predictions

        # Protein domains
        domains = chosen.get("domains")
        if domains:
            result["protein_domains"] = [
                f"{d.get('db', '')}:{d.get('name', '')}" for d in domains[:5]
            ]

        # UniProt
        if chosen.get("swissprot"):
            result["uniprot_id"] = chosen["swissprot"][0] if isinstance(chosen["swissprot"], list) else chosen["swissprot"]

    # Regulatory consequences
    reg_cons = variant.get("regulatory_feature_consequences", [])
    if reg_cons:
        result["regulatory_consequences"] = [
            {
                "regulatory_feature_id": rc.get("regulatory_feature_id", ""),
                "biotype": rc.get("biotype", ""),
                "consequence_terms": rc.get("consequence_terms", []),
                "impact": rc.get("impact", ""),
            }
            for rc in reg_cons[:5]
        ]

    # Motif consequences
    motif_cons = variant.get("motif_feature_consequences", [])
    if motif_cons:
        result["motif_consequences"] = [
            {
                "motif_feature_id": mc.get("motif_feature_id", ""),
                "consequence_terms": mc.get("consequence_terms", []),
                "impact": mc.get("impact", ""),
            }
            for mc in motif_cons[:5]
        ]

    # Colocated variants (population frequencies)
    colocated = variant.get("colocated_variants", [])
    if colocated:
        coloc_summary = []
        for cv in colocated[:3]:
            entry = {"id": cv.get("id", "")}
            freqs = cv.get("frequencies", {})
            if freqs:
                # Get global frequency from first allele
                for allele, allele_freqs in freqs.items():
                    if "gnomade" in allele_freqs:
                        entry["gnomad_global"] = allele_freqs["gnomade"]
                    break
            coloc_summary.append(entry)
        result["colocated_variants"] = coloc_summary

    return result


async def _tool_predict_variant_consequence(args: dict) -> dict:
    """Predict variant consequence using Ensembl VEP REST API."""
    variant = args.get("variant", "").strip()
    if not variant:
        return {"error": "variant is required"}

    fmt = _detect_variant_format(variant)

    # Build URL based on format
    if fmt == "rsid":
        url = f"{VEP_BASE}/id/{variant}"
    elif fmt == "region":
        # Convert "9:22125504:G:C" to "9:22125504-22125504:1/C" for the region endpoint
        parts = variant.split(":")
        if len(parts) == 4:
            chrom, pos, ref, alt = parts
            try:
                pos_int = int(pos)
            except ValueError:
                return {"error": f"Invalid position in variant '{variant}': '{pos}' is not numeric"}
            url = f"{VEP_BASE}/region/{chrom}:{pos_int}-{pos_int + len(ref) - 1}:1/{alt}"
        else:
            # Already in region format or close enough
            url = f"{VEP_BASE}/region/{variant}"
    else:
        url = f"{VEP_BASE}/hgvs/{quote(variant, safe='')}"

    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = dict(VEP_PARAMS)

    try:
        async with aiohttp.ClientSession() as http:
            async with http.get(
                url, params=params, headers=headers, timeout=aiohttp.ClientTimeout(total=30)
            ) as resp:
                if resp.status == 400:
                    body = await resp.text()
                    return {"error": f"Ensembl VEP error (invalid input): {body[:300]}"}
                if resp.status == 429:
                    return {"error": "Ensembl VEP rate limit reached. Please try again in a moment."}
                if resp.status != 200:
                    body = await resp.text()
                    return {"error": f"Ensembl VEP returned status {resp.status}: {body[:300]}"}
                data = await resp.json()
    except Exception as exc:
        return {"error": f"Ensembl VEP API error: {str(exc)}"}

    result = _extract_canonical_consequence(data)
    result["ensembl_vep_url"] = "https://www.ensembl.org/Homo_sapiens/Tools/VEP"
    return result


async def _tool_batch_variant_consequences(args: dict) -> dict:
    """Batch-predict variant consequences using Ensembl VEP POST endpoints."""
    variants = args.get("variants", [])
    if not variants:
        return {"error": "variants list is required"}
    if len(variants) > 200:
        return {"error": f"Maximum 200 variants per batch, got {len(variants)}"}

    # Group variants by format for batch POST
    groups: dict[str, list[str]] = {"rsid": [], "hgvs": [], "region": []}
    for v in variants:
        fmt = _detect_variant_format(str(v).strip())
        groups[fmt].append(str(v).strip())

    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = dict(VEP_PARAMS)

    all_results = []

    async with aiohttp.ClientSession() as http:
        # POST rsIDs
        if groups["rsid"]:
            url = f"{VEP_BASE}/id"
            body = {"ids": groups["rsid"]}
            try:
                async with http.post(
                    url,
                    json=body,
                    params=params,
                    headers=headers,
                    timeout=aiohttp.ClientTimeout(total=60),
                ) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        all_results.extend(data)
                    else:
                        body_text = await resp.text()
                        all_results.append(
                            {
                                "error": f"VEP rsID batch error ({resp.status}): {body_text[:200]}"
                            }
                        )
            except Exception as exc:
                all_results.append(
                    {"error": f"VEP rsID batch error: {str(exc)}"}
                )

        # POST HGVS
        if groups["hgvs"]:
            url = f"{VEP_BASE}/hgvs"
            body = {"hgvs_notations": groups["hgvs"]}
            try:
                async with http.post(
                    url,
                    json=body,
                    params=params,
                    headers=headers,
                    timeout=aiohttp.ClientTimeout(total=60),
                ) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        all_results.extend(data)
                    else:
                        body_text = await resp.text()
                        all_results.append(
                            {
                                "error": f"VEP HGVS batch error ({resp.status}): {body_text[:200]}"
                            }
                        )
            except Exception as exc:
                all_results.append(
                    {"error": f"VEP HGVS batch error: {str(exc)}"}
                )

        # POST regions — need to convert "9:22125504:G:C" to VCF-like format "9 22125504 . G C . . ."
        if groups["region"]:
            vcf_lines = []
            for v in groups["region"]:
                parts = v.split(":")
                if len(parts) == 4:
                    chrom, pos, ref, alt = parts
                    vcf_lines.append(f"{chrom} {pos} . {ref} {alt} . . .")
                else:
                    vcf_lines.append(v)
            url = f"{VEP_BASE}/region"
            body = {"variants": vcf_lines}
            try:
                async with http.post(
                    url,
                    json=body,
                    params=params,
                    headers=headers,
                    timeout=aiohttp.ClientTimeout(total=60),
                ) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        all_results.extend(data)
                    else:
                        body_text = await resp.text()
                        all_results.append(
                            {
                                "error": f"VEP region batch error ({resp.status}): {body_text[:200]}"
                            }
                        )
            except Exception as exc:
                all_results.append(
                    {"error": f"VEP region batch error: {str(exc)}"}
                )

    # Build summary table
    summaries = []
    errors = []
    for item in all_results:
        if "error" in item:
            errors.append(item["error"])
            continue
        row = _extract_canonical_consequence([item])
        summaries.append(row)

    output = {
        "variant_count": len(summaries),
        "annotations": summaries,
    }
    if errors:
        output["errors"] = errors

    return output


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
    "find_datasets_for_target": _tool_find_datasets_for_target,
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
    "get_functional_enrichment": _tool_get_functional_enrichment,
    "predict_variant_consequence": _tool_predict_variant_consequence,
    "batch_variant_consequences": _tool_batch_variant_consequences,
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
    "create_volcano_plot": "Building volcano plot...",
    "create_gene_interaction_network": "Building gene interaction network...",
    "get_string_interactions": "Fetching STRING interaction network...",
    "get_functional_enrichment": "Running STRING functional enrichment...",
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
- "What is the consequence of variant X?" / "VEP for X" / "Predict effect of rs123" / "Is this a splice variant?" / "Regulatory impact of variant" → call predict_variant_consequence (Ensembl VEP)
- "Annotate these variants" / "VEP for this list" / "Consequences of these SNPs" → call batch_variant_consequences (Ensembl VEP)
- "Show MAVE data for X" / "Variant effects for X" / "Heatmap for X" / "Deep mutational scanning" → call create_mave_heatmap
- "Network for X" / "Gene interactions for X" / "What genes are affected by X?" / "Show interaction network" → call create_gene_interaction_network
- "What proteins interact with X?" / "STRING network for X" / "Interaction partners of X" / "PPI network" → call get_string_interactions (STRING)
- "What pathways are enriched?" / "Enrichment analysis for these genes" / "GO terms for X, Y, Z" → call get_functional_enrichment (STRING), then display as table
- "What data do you have for X?" / "Show me data for X" / "Datasets for X" / "What experiments studied X?" → call find_datasets_for_target

IMPORTANT — NEVER FABRICATE DATA:
- NEVER invent or guess dataset IDs, titles, or statistics. ALL dataset information MUST come from tool calls.
- When users ask what data is available for a gene, ALWAYS call find_datasets_for_target first.
- Present the returned datasets as a table using create_visualization, then summarize.

Visualization guidelines:
- Use pie charts for distributions with 2-6 categories
- Use bar charts for comparisons or >6 categories
- Use tables for detailed data rows (limit to 20 rows for readability)
- For volcano plots of Perturb-seq DEA data, use the dedicated create_volcano_plot tool (NOT create_visualization)
- For MAVE variant effect heatmaps, use the dedicated create_mave_heatmap tool (NOT create_visualization)
- For gene interaction networks, use the dedicated create_gene_interaction_network tool (NOT create_visualization)
- When showing datasets, include dataset_id, title, modality, and key metadata
- When showing perturbation data, highlight significant results

VOLCANO PLOT:
When a user asks about Perturb-seq differential expression for a perturbation (e.g., "Show me the effect of knocking out BRCA2", "Volcano plot for TP53"), call create_volcano_plot with the gene name. The backend queries the database and renders the plot directly — do NOT query with query_perturbation_data separately.

Two-step flow:
1. If the gene has data in MULTIPLE datasets, the tool returns a list of available datasets with their significance counts. You MUST present these datasets to the user and ask which one they want to visualize. Then call create_volcano_plot again with the chosen dataset_id.
2. If the gene has data in only ONE dataset, the plot is rendered immediately and you receive a summary with counts of up/down/not-significant genes — use this to write an informative interpretation.

MAVE HEATMAP:
When a user asks about MAVE variant effects, functional scores, or deep mutational scanning results (e.g., "Show me MAVE data for BRCA1", "What are the variant effects for TP53?", "Heatmap for VKORC1"), call create_mave_heatmap with the gene name. The backend queries the database and renders the heatmap directly — do NOT query with query_perturbation_data separately.

Two-step flow (same as volcano plot):
1. If the gene has MAVE data in MULTIPLE datasets, the tool returns a list of available datasets with their variant counts. You MUST present these datasets to the user and ask which one they want to visualize. Then call create_mave_heatmap again with the chosen dataset_id.
2. If the gene has data in only ONE dataset, the heatmap is rendered immediately and you receive a summary — use this to write an informative interpretation about the variant effect landscape.

The heatmap shows 30 positions by default. If the user wants to see a specific region, pass position_start and position_end parameters. Tell the user the displayed position range and total available positions so they can request other regions.

GENE INTERACTION NETWORK:
When a user asks for a network view, gene interactions, or wants to visualize which genes are affected by a perturbation (e.g., "Show me a network for TP53", "What genes interact with BRCA2 perturbation?", "Network of TP53 effects"), call create_gene_interaction_network with the gene name. The backend queries Perturb-seq DEA data and renders an interactive Cytoscape.js network — do NOT query with query_perturbation_data separately.

Two-step flow (same as volcano plot and MAVE heatmap):
1. If the gene has data in MULTIPLE datasets, the tool returns a list of available datasets. Present them to the user and ask which one to visualize. Then call create_gene_interaction_network again with the chosen dataset_id.
2. If the gene has data in only ONE dataset, the network is rendered immediately. You receive a summary with counts of up/down-regulated genes and top affected genes — use this to write an informative interpretation.

The network shows the perturbed gene at the center with significant DEGs radiating outward. Edge width reflects effect magnitude (|log2FC|). Nodes are colored red (upregulated) or blue (downregulated). Default filters: padj < 0.05, |log2FC| >= 0.5, max 30 genes.

Available data modalities:
- Perturb-seq: Single-cell transcriptomic readout of gene perturbations. Key fields: perturbation gene, effect gene, log2FC, padj
- CRISPR screen: Fitness/viability screens. Key fields: perturbation gene, score name, score value, significant
- MAVE: Multiplexed Assay of Variant Effect. Key fields: perturbation gene, variant, position, score

When users ask what data or datasets are available for a gene, call find_datasets_for_target to get real dataset IDs. Use search_target_summary for aggregated overview (modalities, tissues, diseases studied). Never guess dataset IDs.

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
6. Use search_literature for relevant publications
7. Use get_string_interactions for protein-protein interaction context

STRING PROTEIN INTERACTIONS:
You have access to STRING database tools for protein-protein interaction networks and functional enrichment.

- get_string_interactions: Fetch protein-protein interactions for a gene from STRING. Renders an interactive network in the data portal with the query gene at center and partners radiating outward, edge width encoding confidence. Use when users ask "What proteins interact with X?", "Show STRING network for X", "interactors of X". Default required_score=400 (medium confidence); suggest 700 for high confidence or 900 for highest.

- get_functional_enrichment: Run GO/KEGG/Reactome/Pfam enrichment on a gene set using STRING. Returns enriched terms with FDR-corrected p-values. Use when users ask "What pathways are enriched in these genes?", "What biological processes do X, Y, Z share?", or to interpret a set of perturbation hits. After receiving results, display them as a table using create_visualization with columns: Category, Description, Gene Count, FDR, Matching Genes.

STRING workflow:
- For a single gene: call get_string_interactions to show the PPI network
- For a gene set (e.g., top DEGs from a Perturb-seq experiment): call get_functional_enrichment
- Distinguish STRING interaction data (known PPI from databases/literature) from Perturb-seq DEA data (experimental perturbation effects) in your response

ENSEMBL VEP (VARIANT CONSEQUENCE PREDICTION):
You have access to Ensembl VEP for comprehensive variant consequence prediction.

- predict_variant_consequence: Predict the functional consequence of ANY variant type (SNPs, indels, frameshifts, splice variants, regulatory). Returns consequence terms, impact severity (HIGH/MODERATE/LOW/MODIFIER), affected gene/transcript, protein change, and in-silico predictions: SIFT, PolyPhen, CADD, SpliceAI (splicing), AlphaMissense (missense pathogenicity), LOFTEE (loss-of-function), conservation. Also reports regulatory feature consequences (enhancer/promoter/CTCF disruption) and colocated known variants with gnomAD frequencies. Accepts rsIDs (rs56116432), HGVS (ENST00000366667:c.803C>T), or VCF-style (9:22125504:G:C).

- batch_variant_consequences: Annotate up to 200 variants in one call. Efficient for processing variant lists from screening results. Returns a summary per variant. Present results as a table using create_visualization.

VEP vs ProtVar guidance:
- Use VEP (predict_variant_consequence) for: consequence type prediction, splicing impact (SpliceAI), loss-of-function (LOFTEE), regulatory/motif consequences, any variant type (indels, frameshifts, not just missense), novel variants
- Use ProtVar (annotate_variant) for: protein stability (FoldX ddG), structural context (binding pockets, PPI interfaces), EVE/ESM-1b scores, PTM disruption
- For missense variants, both tools are complementary — VEP gives consequence + CADD + SpliceAI + regulatory context, ProtVar gives stability + structural context"""

# --- SSE streaming endpoint ---


def _sse_event(event: str, data: dict) -> str:
    return f"event: {event}\ndata: {json.dumps(data)}\n\n"


async def _persist_and_emit_viz(session_id: str, viz_data: dict) -> str:
    """Persist a visualization to DB and return SSE event string with viz_id included."""
    viz_id = await db_pools["pg"].fetchval(
        """INSERT INTO chat_visualizations (session_id, viz_type, title, data)
           VALUES ($1, $2, $3, $4) RETURNING id""",
        uuid.UUID(session_id),
        viz_data.get("type", ""),
        viz_data.get("title", ""),
        json.dumps(viz_data.get("data", {})),
    )
    viz_data["viz_id"] = viz_id
    return _sse_event("visualization", viz_data)


@router.post("/stream")
async def chat_stream(request: ChatRequest, user: dict = Depends(get_current_user)):
    message = request.message
    user_id = user["id"]

    async def generate():
        nonlocal message
        yield _sse_event("thinking", {"status": "Understanding your question..."})

        try:
            # --- Session creation / ownership verification ---
            session_id = request.session_id
            is_new_session = False

            if not session_id:
                session_id = str(uuid.uuid4())
                await db_pools["pg"].execute(
                    "INSERT INTO chat_sessions (id, user_id) VALUES ($1, $2)",
                    uuid.UUID(session_id),
                    user_id,
                )
                is_new_session = True
            else:
                row = await db_pools["pg"].fetchrow(
                    "SELECT user_id FROM chat_sessions WHERE id = $1",
                    uuid.UUID(session_id),
                )
                if not row or row["user_id"] != user_id:
                    yield _sse_event("error", {"message": "Session not found"})
                    yield _sse_event("done", {"session_id": session_id})
                    return

            # --- Load history from DB ---
            msg_rows = await db_pools["pg"].fetch(
                "SELECT role, content FROM chat_messages WHERE session_id = $1 ORDER BY id",
                uuid.UUID(session_id),
            )
            history = [
                {"role": r["role"], "parts": [{"text": r["content"]}]}
                for r in msg_rows
            ]

            # --- Save user message to DB ---
            await db_pools["pg"].execute(
                "INSERT INTO chat_messages (session_id, role, content) VALUES ($1, $2, $3)",
                uuid.UUID(session_id),
                "user",
                message,
            )

            # --- Auto-title on first message ---
            if is_new_session:
                title = message[:100].strip()
                if len(message) > 100:
                    title = title.rsplit(" ", 1)[0] + "..."
                await db_pools["pg"].execute(
                    "UPDATE chat_sessions SET title = $1 WHERE id = $2",
                    title,
                    uuid.UUID(session_id),
                )

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

            # Build conversation contents from history
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
                        yield await _persist_and_emit_viz(session_id, viz_data)
                        function_response_parts.append(
                            types.Part.from_function_response(
                                name=tool_name,
                                response={
                                    "status": "success",
                                    "message": "Visualization displayed to user",
                                },
                            )
                        )
                    elif tool_name == "create_volcano_plot":
                        # create_volcano_plot: backend queries DB, emits viz, returns summary to LLM
                        yield _sse_event(
                            "tool_call",
                            {"tool": tool_name, "description": "Building volcano plot..."},
                        )
                        yield _sse_event(
                            "thinking", {"status": "Building volcano plot..."}
                        )
                        try:
                            result = await _tool_create_volcano_plot(tool_args)
                        except Exception as exc:
                            logger.exception("Tool %s failed", tool_name)
                            result = {"error": str(exc)}

                        if "error" in result:
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result,
                                )
                            )
                        elif result.get("action") == "choose_dataset":
                            # Multiple datasets — pass list to Gemini to ask the user
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result,
                                )
                            )
                        else:
                            # Plot built — emit visualization and return summary
                            yield await _persist_and_emit_viz(session_id, result["viz_data"])
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result["summary"],
                                )
                            )
                    elif tool_name == "create_mave_heatmap":
                        # create_mave_heatmap: backend queries DB, emits viz, returns summary to LLM
                        yield _sse_event(
                            "tool_call",
                            {"tool": tool_name, "description": "Building MAVE heatmap..."},
                        )
                        yield _sse_event(
                            "thinking", {"status": "Building MAVE heatmap..."}
                        )
                        try:
                            result = await _tool_create_mave_heatmap(tool_args)
                        except Exception as exc:
                            logger.exception("Tool %s failed", tool_name)
                            result = {"error": str(exc)}

                        if "error" in result:
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result,
                                )
                            )
                        elif result.get("action") == "choose_dataset":
                            # Multiple datasets — pass list to Gemini to ask the user
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result,
                                )
                            )
                        else:
                            # Heatmap built — emit visualization and return summary
                            yield await _persist_and_emit_viz(session_id, result["viz_data"])
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result["summary"],
                                )
                            )
                    elif tool_name == "create_gene_interaction_network":
                        # create_gene_interaction_network: backend queries DB, emits viz, returns summary to LLM
                        yield _sse_event(
                            "tool_call",
                            {"tool": tool_name, "description": "Building gene interaction network..."},
                        )
                        yield _sse_event(
                            "thinking", {"status": "Building gene interaction network..."}
                        )
                        try:
                            result = await _tool_create_gene_interaction_network(tool_args)
                        except Exception as exc:
                            logger.exception("Tool %s failed", tool_name)
                            result = {"error": str(exc)}

                        if "error" in result:
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result,
                                )
                            )
                        elif result.get("action") == "choose_dataset":
                            # Multiple datasets — pass list to Gemini to ask the user
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result,
                                )
                            )
                        else:
                            # Network built — emit visualization and return summary
                            yield await _persist_and_emit_viz(session_id, result["viz_data"])
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result["summary"],
                                )
                            )
                    elif tool_name == "get_string_interactions":
                        yield _sse_event(
                            "tool_call",
                            {"tool": tool_name, "description": "Fetching STRING interaction network..."},
                        )
                        yield _sse_event(
                            "thinking", {"status": "Fetching STRING interaction network..."}
                        )
                        try:
                            result = await _tool_get_string_interactions(tool_args)
                        except Exception as exc:
                            logger.exception("Tool %s failed", tool_name)
                            result = {"error": str(exc)}

                        if "error" in result:
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result,
                                )
                            )
                        else:
                            yield await _persist_and_emit_viz(session_id, result["viz_data"])
                            function_response_parts.append(
                                types.Part.from_function_response(
                                    name=tool_name,
                                    response=result["summary"],
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
                            yield await _persist_and_emit_viz(session_id, viz_data)

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

            # Save assistant text response to DB (last model turn)
            assistant_text = ""
            for c in reversed(contents):
                if c.role == "model":
                    for p in c.parts:
                        if p.text:
                            assistant_text += p.text
                    break
            if assistant_text:
                await db_pools["pg"].execute(
                    "INSERT INTO chat_messages (session_id, role, content) VALUES ($1, $2, $3)",
                    uuid.UUID(session_id),
                    "model",
                    assistant_text,
                )

            # Update session timestamp
            await db_pools["pg"].execute(
                "UPDATE chat_sessions SET updated_at = NOW() WHERE id = $1",
                uuid.UUID(session_id),
            )

        except Exception as exc:
            logger.exception("Chat stream error")
            yield _sse_event("error", {"message": str(exc)})

        yield _sse_event("done", {"session_id": session_id})

    return StreamingResponse(generate(), media_type="text/event-stream")
