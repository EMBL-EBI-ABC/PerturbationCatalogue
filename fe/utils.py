"""Shared utilities for the Dash app"""

import requests
from typing import Dict, List, Any, Optional, Tuple
import json
from dash import html, dcc
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
import os
import time

# Backend API URL - can be set via environment variable for production
BACKEND_URL = os.getenv("PERTURBATION_CATALOGUE_BE", "https://perturbation-catalogue-be-328296435987.europe-west2.run.app")

if not BACKEND_URL:
    raise ValueError("PERTURBATION_CATALOGUE_BE environment variable is not set.")

# Facet fields to exclude from stats display
FACET_FIELDS = [
    "data_modalities",
    "tissues_tested",
    "cell_types_tested",
    "cell_lines_tested",
    "sex_tested",
    "developmental_stages_tested",
    "diseases_tested",
]

# Designer color scheme
COLORS = {
    "primary": "#007B53",
    "secondary": "#193F90",
    "accent1": "#563D82",
    "accent2": "#3B6FB6",
    "gray": "#54585A",
    "dark_green": "#0A5032",
    "red": "#A6093D",
}

DATA_MODALITIES_COLOURS = {
    "Perturb-seq": COLORS["primary"],
    "CRISPR screen": COLORS["secondary"],
    "MAVE": COLORS["red"],
}

# Color mapping for each filter category (for bar charts)
FACET_COLORS = {
    "data_modalities": "#007B53",  # primary green
    "tissues_tested": "#193F90",  # secondary blue
    "cell_types_tested": "#563D82",  # accent1 purple
    "cell_lines_tested": "#3B6FB6",  # accent2 light blue
    "sex_tested": "#A6093D",  # red
    "developmental_stages_tested": "#0A5032",  # dark green
    "diseases_tested": "#54585A",  # gray
}

# Lighter versions for unselected bars (for better contrast)
FACET_COLORS_LIGHT = {
    "data_modalities": "#4CAF7A",  # lighter green
    "tissues_tested": "#5B7FC7",  # lighter blue
    "cell_types_tested": "#8B6FA8",  # lighter purple
    "cell_lines_tested": "#6B9BD4",  # lighter light blue
    "sex_tested": "#C84A6D",  # lighter red
    "developmental_stages_tested": "#3D7A5C",  # lighter dark green
    "diseases_tested": "#7A7E80",  # lighter gray
}

SUMMARY_CACHE_TTL = 60  # seconds

_summary_cache = {"data": None, "timestamp": 0.0}


# Store for search results to access when clicking View More
results_store = {}


def get_landing_page_summary(
    ttl: int = SUMMARY_CACHE_TTL,
) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """Fetch and cache landing page summary from backend."""
    now = time.time()
    cached_summary = _summary_cache.get("data")
    if cached_summary is not None and now - _summary_cache.get("timestamp", 0.0) < ttl:
        return cached_summary, None

    try:
        response = requests.get(f"{BACKEND_URL}/summary", timeout=10)
        response.raise_for_status()
        summary = response.json()
        _summary_cache["data"] = summary
        _summary_cache["timestamp"] = now
        return summary, None
    except Exception as exc:
        error_message = f"Error fetching summary: {exc}"
        print(error_message)
        return None, error_message


# Helper function to fetch data from backend
def fetch_search_results(
    query: Optional[str] = None,
    filters: Optional[Dict[str, List[str]]] = None,
    page: int = 1,
    size: int = 6,
) -> Dict[str, Any]:
    """Fetch search results from backend API"""
    try:
        params = {"page": page, "size": size}

        if query:
            params["query"] = query

        if filters:
            for field, values in filters.items():
                if values:
                    params[field] = ",".join(values)

        response = requests.get(f"{BACKEND_URL}/search", params=params, timeout=10)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error fetching data: {e}")
        return {
            "total": 0,
            "page": 1,
            "size": size,
            "total_pages": 0,
            "results": [],
            "facets": {
                "data_modalities": [],
                "tissues_tested": [],
                "cell_types_tested": [],
                "cell_lines_tested": [],
                "sex_tested": [],
                "developmental_stages_tested": [],
                "diseases_tested": [],
            },
        }


# Helper functions to fetch table data from backend
def fetch_perturb_seq(
    perturbed_target_symbol: str,
    gene: Optional[str] = None,
    sort_by: Optional[str] = None,
    sort_order: str = "asc",
    page: int = 1,
    size: int = 10,
) -> Dict[str, Any]:
    """Fetch perturb_seq data from backend API"""
    try:
        params = {
            "perturbed_target_symbol": perturbed_target_symbol,
            "page": page,
            "size": size,
        }
        if gene:
            params["gene"] = gene
        if sort_by:
            params["sort_by"] = sort_by
            params["sort_order"] = sort_order

        response = requests.get(f"{BACKEND_URL}/perturb_seq", params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error fetching perturb_seq data: {e}")
        return {"total": 0, "page": 1, "size": size, "total_pages": 0, "results": []}


def fetch_mave_data(
    perturbed_target_symbol: str, page: int = 1, size: int = 10
) -> Dict[str, Any]:
    """Fetch mave_data from backend API"""
    try:
        params = {
            "perturbed_target_symbol": perturbed_target_symbol,
            "page": page,
            "size": size,
        }
        response = requests.get(f"{BACKEND_URL}/mave_data", params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error fetching mave_data: {e}")
        return {"total": 0, "page": 1, "size": size, "total_pages": 0, "results": []}


def fetch_crispr_data(
    perturbed_target_symbol: str, page: int = 1, size: int = 10
) -> Dict[str, Any]:
    """Fetch crispr_data from backend API"""
    try:
        params = {
            "perturbed_target_symbol": perturbed_target_symbol,
            "page": page,
            "size": size,
        }
        response = requests.get(f"{BACKEND_URL}/crispr_data", params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error fetching crispr_data: {e}")
        return {"total": 0, "page": 1, "size": size, "total_pages": 0, "results": []}


def fetch_target_details(perturbation_gene_name: str) -> Dict[str, Any]:
    """Fetch detailed perturb-seq datasets for a target"""
    if not perturbation_gene_name:
        return {"datasets": []}

    try:
        params = {"perturbation_gene_name": perturbation_gene_name}
        response = requests.get(
            f"{BACKEND_URL}/v1/perturb-seq/search", params=params, timeout=30
        )
        response.raise_for_status()
        return response.json()
    except Exception as exc:
        error_message = (
            f"Error fetching target details for {perturbation_gene_name}: {exc}"
        )
        print(error_message)
        return {"datasets": [], "error": error_message}


def fetch_modality_datasets(
    modality: str,
    filters: Optional[Dict[str, Any]] = None,
    dataset_limit: int = 4,
    dataset_offset: int = 0,
    rows_per_dataset_limit: int = 5,
    sort: Optional[str] = None,
) -> Dict[str, Any]:
    """Call /v1/{modality}/search with the provided filters."""
    try:
        params: Dict[str, Any] = {
            "dataset_limit": dataset_limit,
            "dataset_offset": dataset_offset,
            "rows_per_dataset_limit": rows_per_dataset_limit,
        }
        if sort:
            params["sort"] = sort

        if filters:
            for key, value in filters.items():
                if value is None:
                    continue
                params[key] = value

        response = requests.get(
            f"{BACKEND_URL}/v1/{modality}/search", params=params, timeout=30
        )
        response.raise_for_status()
        return response.json()
    except Exception as exc:
        error_message = f"Error fetching {modality} datasets: {exc}"
        print(error_message)
        return {"datasets": [], "error": error_message}


def fetch_dataset_rows(
    modality: str,
    dataset_id: str,
    filters: Optional[Dict[str, Any]] = None,
    limit: int = 5,
    offset: int = 0,
    sort: Optional[str] = None,
) -> Dict[str, Any]:
    """Call /v1/{modality}/{dataset_id}/search to fetch more rows."""
    try:
        params: Dict[str, Any] = {"limit": limit, "offset": offset}
        if sort:
            params["sort"] = sort

        if filters:
            for key, value in filters.items():
                if value is None:
                    continue
                params[key] = value

        response = requests.get(
            f"{BACKEND_URL}/v1/{modality}/{dataset_id}/search",
            params=params,
            timeout=30,
        )
        response.raise_for_status()
        return response.json()
    except Exception as exc:
        error_message = (
            f"Error fetching dataset {dataset_id} for modality {modality}: {exc}"
        )
        print(error_message)
        return {"results": [], "error": error_message, "total_rows_count": 0}


# Helper function to format value for display
def format_value(value: Any) -> str:
    """Format a value for display"""
    if isinstance(value, list):
        return ", ".join(str(v) for v in value) if value else "N/A"
    elif isinstance(value, dict):
        return json.dumps(value, indent=2)
    elif value is None:
        return "N/A"
    else:
        return str(value)


# Helper function to format number
def format_number(num):
    """Format number for display"""
    if isinstance(num, (int, float)):
        if abs(num) < 0.001:
            return f"{num:.2e}"
        elif abs(num) < 1:
            return f"{num:.4f}"
        else:
            return f"{num:.2f}"
    return str(num) if num is not None else "N/A"


# Helper function to create card content with organized stats
def create_card_content(result: Dict[str, Any]) -> List:
    """Create card content with organized stats and tables"""

    content = []

    # Get n_experiments at the beginning - full width on mobile
    n_experiments = result.get("n_experiments", "N/A")
    content.append(
        html.Div(
            [
                html.Strong("Total Experiments: "),
                html.Span(format_value(n_experiments)),
            ],
            className="mb-3 w-100",
            style={"fontSize": "16px"},
        )
    )

    # Group 1: RNA-seq related stats
    rnaseq_stats = []
    rnaseq_fields = [
        "n_sig_perturb_pairs_up",
        "n_sig_perturb_pairs_down",
        "lfc_q25",
        "lfc_median",
        "lfc_q75",
    ]
    for field in rnaseq_fields:
        if field in result:
            display_key = field.replace("_", " ").title()
            rnaseq_stats.append(
                html.Div(
                    [
                        html.Strong(f"{display_key}: "),
                        html.Span(format_value(result[field])),
                    ],
                    className="mb-1",
                )
            )

    # Helper function to parse and format a gene pair entry
    def parse_gene_pair(pair):
        """Parse a gene pair entry - can be dict or string"""
        if isinstance(pair, dict):
            return {
                "gene": pair.get("gene", "N/A"),
                "padj": pair.get("padj", "N/A"),
                "log2foldchange": pair.get("log2foldchange", "N/A"),
            }
        elif isinstance(pair, str):
            try:
                parsed = json.loads(pair)
                return {
                    "gene": parsed.get("gene", "N/A"),
                    "padj": parsed.get("padj", "N/A"),
                    "log2foldchange": parsed.get("log2foldchange", "N/A"),
                }
            except:
                return {"gene": pair, "padj": "N/A", "log2foldchange": "N/A"}
        return {"gene": str(pair), "padj": "N/A", "log2foldchange": "N/A"}

    # Significant Gene Pairs data
    sig_pairs_up = result.get("sig_gene_pairs_up", [])
    sig_pairs_down = result.get("sig_gene_pairs_down", [])

    # Build RNA-seq section with stats and gene pairs table
    if rnaseq_stats or sig_pairs_up or sig_pairs_down:
        rnaseq_section = []
        rnaseq_section.append(
            html.Div(
                [
                    html.H6(
                        "RNA-seq Statistics",
                        style={
                            "color": COLORS["secondary"],
                            "marginBottom": "10px",
                            "fontWeight": "600",
                        },
                    ),
                    html.Div(rnaseq_stats),
                ],
                className="mb-3 stats-box",
                style={"padding": "16px"},
            )
        )

        # Add Significant Gene Pairs table under RNA-seq Statistics
        if sig_pairs_up or sig_pairs_down:
            # Create table for up genes
            up_table_rows = []
            up_table_rows.append(
                html.Tr(
                    [
                        html.Th(
                            "Gene",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["primary"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            "p-adj",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["primary"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            "log2FoldChange",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["primary"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            html.I(className="bi bi-arrow-up"),
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["primary"],
                                "color": "white",
                                "textAlign": "center",
                            },
                        ),
                    ]
                )
            )

            for pair in sig_pairs_up[:10]:  # Limit to 10 for display
                parsed = parse_gene_pair(pair)
                up_table_rows.append(
                    html.Tr(
                        [
                            html.Td(
                                parsed["gene"],
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                format_number(parsed["padj"]),
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                format_number(parsed["log2foldchange"]),
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                [
                                    html.I(
                                        className="bi bi-arrow-up",
                                        style={
                                            "color": COLORS["primary"],
                                            "fontSize": "18px",
                                            "fontWeight": "bold",
                                        },
                                    )
                                ],
                                style={
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "textAlign": "center",
                                },
                            ),
                        ]
                    )
                )

            if len(sig_pairs_up) > 10:
                up_table_rows.append(
                    html.Tr(
                        [
                            html.Td(
                                f"... ({len(sig_pairs_up) - 10} more)",
                                colSpan=4,
                                style={
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "textAlign": "center",
                                    "fontStyle": "italic",
                                },
                            )
                        ]
                    )
                )

            # Create table for down genes
            down_table_rows = []
            down_table_rows.append(
                html.Tr(
                    [
                        html.Th(
                            "Gene",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["red"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            "p-adj",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["red"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            "log2FoldChange",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["red"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            html.I(className="bi bi-arrow-down"),
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["red"],
                                "color": "white",
                                "textAlign": "center",
                            },
                        ),
                    ]
                )
            )

            for pair in sig_pairs_down[:10]:  # Limit to 10 for display
                parsed = parse_gene_pair(pair)
                down_table_rows.append(
                    html.Tr(
                        [
                            html.Td(
                                parsed["gene"],
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                format_number(parsed["padj"]),
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                format_number(parsed["log2foldchange"]),
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                [
                                    html.I(
                                        className="bi bi-arrow-down",
                                        style={
                                            "color": COLORS["red"],
                                            "fontSize": "18px",
                                            "fontWeight": "bold",
                                        },
                                    )
                                ],
                                style={
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "textAlign": "center",
                                },
                            ),
                        ]
                    )
                )

            if len(sig_pairs_down) > 10:
                down_table_rows.append(
                    html.Tr(
                        [
                            html.Td(
                                f"... ({len(sig_pairs_down) - 10} more)",
                                colSpan=4,
                                style={
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "textAlign": "center",
                                    "fontStyle": "italic",
                                },
                            )
                        ]
                    )
                )

            # Add Significant Gene Pairs table - stacked on mobile
            rnaseq_section.append(
                html.Div(
                    [
                        html.H6(
                            "Significant Gene Pairs",
                            style={"color": COLORS["gray"], "marginBottom": "10px"},
                        ),
                        html.Div(
                            [
                                html.H6(
                                    "Up",
                                    style={
                                        "color": COLORS["primary"],
                                        "fontSize": "14px",
                                        "marginBottom": "5px",
                                    },
                                ),
                                html.Table(
                                    [
                                        html.Thead(up_table_rows[0]),
                                        html.Tbody(up_table_rows[1:]),
                                    ],
                                    style={
                                        "width": "100%",
                                        "borderCollapse": "collapse",
                                        "fontSize": "11px",
                                    },
                                    className="table-responsive mb-3",
                                ),
                            ],
                            className="d-md-none",
                        ),  # Stacked version for mobile
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Div(
                                            [
                                                html.H6(
                                                    "Up",
                                                    style={
                                                        "color": COLORS["primary"],
                                                        "fontSize": "14px",
                                                        "marginBottom": "5px",
                                                        "display": "none",
                                                    },
                                                    className="d-md-block",
                                                ),
                                                html.Table(
                                                    [
                                                        html.Thead(up_table_rows[0]),
                                                        html.Tbody(up_table_rows[1:]),
                                                    ],
                                                    style={
                                                        "width": "100%",
                                                        "borderCollapse": "collapse",
                                                        "fontSize": "11px",
                                                    },
                                                    className="table-responsive",
                                                ),
                                            ]
                                        )
                                    ],
                                    width={"size": 6, "md": 6, "sm": 12},
                                    className="mb-3 mb-md-0",
                                ),
                                dbc.Col(
                                    [
                                        html.Div(
                                            [
                                                html.H6(
                                                    "Down",
                                                    style={
                                                        "color": COLORS["red"],
                                                        "fontSize": "14px",
                                                        "marginBottom": "5px",
                                                        "display": "none",
                                                    },
                                                    className="d-md-block",
                                                ),
                                                html.Table(
                                                    [
                                                        html.Thead(down_table_rows[0]),
                                                        html.Tbody(down_table_rows[1:]),
                                                    ],
                                                    style={
                                                        "width": "100%",
                                                        "borderCollapse": "collapse",
                                                        "fontSize": "11px",
                                                    },
                                                    className="table-responsive",
                                                ),
                                            ]
                                        )
                                    ],
                                    width={"size": 6, "md": 6, "sm": 12},
                                ),
                            ],
                            className="d-none d-md-flex",
                        ),  # Side by side for desktop
                        html.Div(
                            [
                                html.H6(
                                    "Down",
                                    style={
                                        "color": COLORS["red"],
                                        "fontSize": "14px",
                                        "marginBottom": "5px",
                                    },
                                ),
                                html.Table(
                                    [
                                        html.Thead(down_table_rows[0]),
                                        html.Tbody(down_table_rows[1:]),
                                    ],
                                    style={
                                        "width": "100%",
                                        "borderCollapse": "collapse",
                                        "fontSize": "11px",
                                    },
                                    className="table-responsive mb-3",
                                ),
                            ],
                            className="d-md-none",
                        ),  # Stacked version for mobile
                    ],
                    className="mb-3",
                )
            )

        content.extend(rnaseq_section)

    # Group 2: CRISPR stats
    crispr_stats = []
    if "n_sig_crisps" in result:
        crispr_stats.append(
            html.Div(
                [
                    html.Strong("N Sig Crisps: "),
                    html.Span(format_value(result["n_sig_crisps"])),
                ],
                className="mb-1",
            )
        )
    n_sig_crispr = result.get("n_sig_crispr", 0)
    crispr_stats.append(
        html.Div(
            [html.Strong("N Sig Crispr: "), html.Span(format_value(n_sig_crispr))],
            className="mb-1",
        )
    )

    if crispr_stats:
        content.append(
            html.Div(
                [
                    html.H6(
                        "CRISPR Screen Statistics",
                        style={
                            "color": COLORS["accent1"],
                            "marginBottom": "10px",
                            "fontWeight": "600",
                        },
                    ),
                    html.Div(crispr_stats),
                ],
                className="mb-3 stats-box",
                style={"padding": "16px"},
            )
        )

    # Group 3: MAVE stats with MAVE scores table
    mave_scores_up = result.get("mave_scores_up", [])
    mave_scores_down = result.get("mave_scores_down", [])

    # Helper function to parse and format a MAVE score entry
    def parse_mave_score(score):
        """Parse a MAVE score entry - can be dict or string"""
        if isinstance(score, dict):
            return {
                "variant": score.get("perturbation_name", "N/A"),
                "score": score.get("score_value", "N/A"),
            }
        elif isinstance(score, str):
            try:
                parsed = json.loads(score)
                return {
                    "variant": parsed.get("perturbation_name", "N/A"),
                    "score": parsed.get("score_value", "N/A"),
                }
            except:
                return {"variant": score, "score": "N/A"}
        return {"variant": str(score), "score": "N/A"}

    # Format numbers for display (for scores)
    def format_score(num):
        """Format score number for display"""
        if isinstance(num, (int, float)):
            if abs(num) < 0.001:
                return f"{num:.2e}"
            elif abs(num) < 1:
                return f"{num:.4f}"
            else:
                return f"{num:.2f}"
        return str(num)

    if "n_mave" in result or mave_scores_up or mave_scores_down:
        mave_section = []
        mave_section.append(
            html.Div(
                [
                    html.H6(
                        "MAVE Statistics",
                        style={
                            "color": COLORS["accent2"],
                            "marginBottom": "10px",
                            "fontWeight": "600",
                        },
                    ),
                    html.Div(
                        [
                            html.Strong("N Mave: "),
                            html.Span(format_value(result.get("n_mave", "N/A"))),
                        ]
                    ),
                ],
                className="mb-3 stats-box",
                style={"padding": "16px"},
            )
        )

        # Add MAVE Scores table under MAVE Statistics
        if mave_scores_up or mave_scores_down:
            # Create table for up scores
            up_table_rows = []
            up_table_rows.append(
                html.Tr(
                    [
                        html.Th(
                            "Variant",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["primary"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            "Score",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["primary"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            html.I(className="bi bi-arrow-up"),
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["primary"],
                                "color": "white",
                                "textAlign": "center",
                            },
                        ),
                    ]
                )
            )

            for score in mave_scores_up[:10]:  # Limit to 10 for display
                parsed = parse_mave_score(score)
                up_table_rows.append(
                    html.Tr(
                        [
                            html.Td(
                                parsed["variant"],
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                format_score(parsed["score"]),
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                [
                                    html.I(
                                        className="bi bi-arrow-up",
                                        style={
                                            "color": COLORS["primary"],
                                            "fontSize": "18px",
                                            "fontWeight": "bold",
                                        },
                                    )
                                ],
                                style={
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "textAlign": "center",
                                },
                            ),
                        ]
                    )
                )

            if len(mave_scores_up) > 10:
                up_table_rows.append(
                    html.Tr(
                        [
                            html.Td(
                                f"... ({len(mave_scores_up) - 10} more)",
                                colSpan=3,
                                style={
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "textAlign": "center",
                                    "fontStyle": "italic",
                                },
                            )
                        ]
                    )
                )

            # Create table for down scores
            down_table_rows = []
            down_table_rows.append(
                html.Tr(
                    [
                        html.Th(
                            "Variant",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["red"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            "Score",
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["red"],
                                "color": "white",
                            },
                        ),
                        html.Th(
                            html.I(className="bi bi-arrow-down"),
                            style={
                                "padding": "8px",
                                "border": "1px solid #ddd",
                                "backgroundColor": COLORS["red"],
                                "color": "white",
                                "textAlign": "center",
                            },
                        ),
                    ]
                )
            )

            for score in mave_scores_down[:10]:  # Limit to 10 for display
                parsed = parse_mave_score(score)
                down_table_rows.append(
                    html.Tr(
                        [
                            html.Td(
                                parsed["variant"],
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                format_score(parsed["score"]),
                                style={"padding": "8px", "border": "1px solid #ddd"},
                            ),
                            html.Td(
                                [
                                    html.I(
                                        className="bi bi-arrow-down",
                                        style={
                                            "color": COLORS["red"],
                                            "fontSize": "18px",
                                            "fontWeight": "bold",
                                        },
                                    )
                                ],
                                style={
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "textAlign": "center",
                                },
                            ),
                        ]
                    )
                )

            if len(mave_scores_down) > 10:
                down_table_rows.append(
                    html.Tr(
                        [
                            html.Td(
                                f"... ({len(mave_scores_down) - 10} more)",
                                colSpan=3,
                                style={
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "textAlign": "center",
                                    "fontStyle": "italic",
                                },
                            )
                        ]
                    )
                )

            # Add MAVE Scores table - stacked on mobile
            mave_section.append(
                html.Div(
                    [
                        html.H6(
                            "MAVE Scores",
                            style={"color": COLORS["gray"], "marginBottom": "10px"},
                        ),
                        html.Div(
                            [
                                html.H6(
                                    "Up",
                                    style={
                                        "color": COLORS["primary"],
                                        "fontSize": "14px",
                                        "marginBottom": "5px",
                                    },
                                ),
                                html.Table(
                                    [
                                        html.Thead(up_table_rows[0]),
                                        html.Tbody(up_table_rows[1:]),
                                    ],
                                    style={
                                        "width": "100%",
                                        "borderCollapse": "collapse",
                                        "fontSize": "11px",
                                    },
                                    className="table-responsive mb-3",
                                ),
                            ],
                            className="d-md-none",
                        ),  # Stacked version for mobile
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Div(
                                            [
                                                html.H6(
                                                    "Up",
                                                    style={
                                                        "color": COLORS["primary"],
                                                        "fontSize": "14px",
                                                        "marginBottom": "5px",
                                                        "display": "none",
                                                    },
                                                    className="d-md-block",
                                                ),
                                                html.Table(
                                                    [
                                                        html.Thead(up_table_rows[0]),
                                                        html.Tbody(up_table_rows[1:]),
                                                    ],
                                                    style={
                                                        "width": "100%",
                                                        "borderCollapse": "collapse",
                                                        "fontSize": "11px",
                                                    },
                                                    className="table-responsive",
                                                ),
                                            ]
                                        )
                                    ],
                                    width={"size": 6, "md": 6, "sm": 12},
                                    className="mb-3 mb-md-0",
                                ),
                                dbc.Col(
                                    [
                                        html.Div(
                                            [
                                                html.H6(
                                                    "Down",
                                                    style={
                                                        "color": COLORS["red"],
                                                        "fontSize": "14px",
                                                        "marginBottom": "5px",
                                                        "display": "none",
                                                    },
                                                    className="d-md-block",
                                                ),
                                                html.Table(
                                                    [
                                                        html.Thead(down_table_rows[0]),
                                                        html.Tbody(down_table_rows[1:]),
                                                    ],
                                                    style={
                                                        "width": "100%",
                                                        "borderCollapse": "collapse",
                                                        "fontSize": "11px",
                                                    },
                                                    className="table-responsive",
                                                ),
                                            ]
                                        )
                                    ],
                                    width={"size": 6, "md": 6, "sm": 12},
                                ),
                            ],
                            className="d-none d-md-flex",
                        ),  # Side by side for desktop
                        html.Div(
                            [
                                html.H6(
                                    "Down",
                                    style={
                                        "color": COLORS["red"],
                                        "fontSize": "14px",
                                        "marginBottom": "5px",
                                    },
                                ),
                                html.Table(
                                    [
                                        html.Thead(down_table_rows[0]),
                                        html.Tbody(down_table_rows[1:]),
                                    ],
                                    style={
                                        "width": "100%",
                                        "borderCollapse": "collapse",
                                        "fontSize": "11px",
                                    },
                                    className="table-responsive mb-3",
                                ),
                            ],
                            className="d-md-none",
                        ),  # Stacked version for mobile
                    ],
                    className="mb-3",
                )
            )

        content.extend(mave_section)

    return content


# Helper function to create stats display (excluding tables) for details page
def create_stats_display(result: Dict[str, Any]) -> List:
    """Create stats display without tables for details page"""

    content = []

    # Get n_experiments at the beginning
    n_experiments = result.get("n_experiments", "N/A")
    content.append(
        html.Div(
            [
                html.Strong("Total Experiments: "),
                html.Span(format_value(n_experiments)),
            ],
            className="mb-3",
            style={"fontSize": "16px"},
        )
    )

    # Group 1: RNA-seq related stats (excluding tables)
    rnaseq_stats = []
    rnaseq_fields = [
        "n_sig_perturb_pairs_up",
        "n_sig_perturb_pairs_down",
        "lfc_q25",
        "lfc_median",
        "lfc_q75",
    ]
    for field in rnaseq_fields:
        if field in result:
            display_key = field.replace("_", " ").title()
            rnaseq_stats.append(
                html.Div(
                    [
                        html.Strong(f"{display_key}: "),
                        html.Span(format_value(result[field])),
                    ],
                    className="mb-1",
                )
            )

    if rnaseq_stats:
        content.append(
            html.Div(
                [
                    html.H6(
                        "RNA-seq Statistics",
                        style={
                            "color": COLORS["secondary"],
                            "marginBottom": "10px",
                            "fontWeight": "600",
                        },
                    ),
                    html.Div(rnaseq_stats),
                ],
                className="mb-3 stats-box",
                style={"padding": "16px"},
            )
        )

    # Group 2: CRISPR stats
    crispr_stats = []
    if "n_sig_crisps" in result:
        crispr_stats.append(
            html.Div(
                [
                    html.Strong("N Sig Crisps: "),
                    html.Span(format_value(result["n_sig_crisps"])),
                ],
                className="mb-1",
            )
        )
    n_sig_crispr = result.get("n_sig_crispr", 0)
    crispr_stats.append(
        html.Div(
            [html.Strong("N Sig Crispr: "), html.Span(format_value(n_sig_crispr))],
            className="mb-1",
        )
    )

    if crispr_stats:
        content.append(
            html.Div(
                [
                    html.H6(
                        "CRISPR Screen Statistics",
                        style={"color": COLORS["accent1"], "marginBottom": "10px"},
                    ),
                    html.Div(crispr_stats),
                ],
                className="mb-3",
                style={
                    "padding": "10px",
                    "backgroundColor": "#f8f9fa",
                    "borderRadius": "5px",
                },
            )
        )

    # Group 3: MAVE stats (excluding tables)
    if "n_mave" in result:
        content.append(
            html.Div(
                [
                    html.H6(
                        "MAVE Statistics",
                        style={
                            "color": COLORS["accent2"],
                            "marginBottom": "10px",
                            "fontWeight": "600",
                        },
                    ),
                    html.Div(
                        [
                            html.Strong("N Mave: "),
                            html.Span(format_value(result.get("n_mave", "N/A"))),
                        ]
                    ),
                ],
                className="mb-3 stats-box",
                style={"padding": "16px"},
            )
        )

    return content


# Helper function to create filter component (kept for backward compatibility if needed)
def create_filter_component(
    facet_name: str, facet_values: List[Dict[str, Any]], selected: List[str]
) -> html.Div:
    """Create a filter component for a facet - showing only top 5 by count"""
    # Format facet name for display
    display_name = facet_name.replace("_", " ").title()

    # Sort by count descending and take top 5
    sorted_values = sorted(facet_values, key=lambda x: x.get("count", 0), reverse=True)[
        :5
    ]

    options = [
        {"label": f"{item['value']} ({item['count']})", "value": item["value"]}
        for item in sorted_values
    ]

    return html.Div(
        [
            html.Label(
                display_name,
                style={
                    "fontWeight": "bold",
                    "marginTop": "10px",
                    "marginBottom": "5px",
                },
            ),
            dcc.Checklist(
                id={"type": "filter", "field": facet_name},
                options=options,
                value=selected,
                inputStyle={"marginRight": "5px"},
                labelStyle={"display": "block", "marginBottom": "5px"},
            ),
        ],
        style={"marginBottom": "20px"},
        className="col-12",
    )  # Full width on mobile, one per row


# Helper function to create chart component (pie or bar chart)
def create_chart_component(
    facet_name: str, facet_values: List[Dict[str, Any]], selected: List[str]
) -> html.Div:
    """Create a chart component for a facet - pie chart if ≤5 values, bar chart otherwise"""
    # Format facet name for display
    display_name = facet_name.replace("_", " ").title()

    if not facet_values:
        return html.Div(
            [
                html.H6(
                    display_name,
                    style={
                        "fontWeight": "600",
                        "marginBottom": "5px",
                        "color": COLORS["primary"],
                        "fontSize": "0.9rem",
                    },
                ),
                html.P(
                    "No data available", style={"color": "#999", "fontSize": "0.8rem"}
                ),
            ],
            style={"width": "100%", "marginBottom": "0.5rem"},
        )

    # Sort by count descending
    sorted_values = sorted(facet_values, key=lambda x: x.get("count", 0), reverse=True)

    # Determine if we should use pie chart (≤5 values) or bar chart (>5 values)
    use_pie = len(sorted_values) <= 5

    if use_pie:
        # Create pie chart
        labels = [item["value"] for item in sorted_values]
        values = [item["count"] for item in sorted_values]
        colors_list = px.colors.qualitative.Set3[: len(labels)]

        # Highlight selected values
        pie_colors = []
        for i, label in enumerate(labels):
            if label in selected:
                pie_colors.append(COLORS["primary"])
            else:
                pie_colors.append(colors_list[i] if i < len(colors_list) else "#cccccc")

        fig = go.Figure(
            data=[
                go.Pie(
                    labels=labels,
                    values=values,
                    hole=0.3,
                    marker=dict(colors=pie_colors),
                    textinfo="percent",
                    textposition="inside",
                    hovertemplate="<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>",
                    textfont=dict(size=10),
                )
            ]
        )

        fig.update_layout(
            title=display_name,
            title_x=0.5,
            showlegend=True,
            height=250,
            margin=dict(t=35, b=10, l=10, r=10),
            font=dict(family="Arial, sans-serif", size=10),
        )

        chart = dcc.Graph(
            figure=fig,
            id={"type": "chart", "field": facet_name},
            config={"displayModeBar": False},
            style={"width": "100%", "height": "100%"},
        )

        # For pie charts showing all values (≤5), no dropdown needed
        # All values are visible in the chart itself
        # Note: To enable filtering for pie charts, chart click events would need to be implemented
        return html.Div([chart], style={"width": "100%", "marginBottom": "0.5rem"})
    else:
        # Create bar chart with top 5 values (reduced from 10 for compactness)
        top_5 = sorted_values[:5]
        rest_values = sorted_values[5:]

        labels = [item["value"] for item in top_5]
        values = [item["count"] for item in top_5]

        # Use category-specific colors for bar charts
        # Each filter category has its own color scheme
        category_color = FACET_COLORS.get(facet_name, COLORS["accent1"])
        category_color_light = FACET_COLORS_LIGHT.get(facet_name, COLORS["accent2"])

        # Selected values use the main category color, unselected use lighter version
        bar_colors = [
            category_color if label in selected else category_color_light
            for label in labels
        ]

        fig = go.Figure(
            data=[
                go.Bar(
                    x=labels,
                    y=values,
                    marker=dict(color=bar_colors),
                    text=values,
                    textposition="outside",
                    hovertemplate="<b>%{x}</b><br>Count: %{y}<extra></extra>",
                )
            ]
        )

        fig.update_layout(
            title=display_name,
            title_x=0.5,
            xaxis_title="",
            yaxis_title="Count",
            height=250,
            margin=dict(t=35, b=60, l=40, r=10),
            font=dict(family="Arial, sans-serif", size=10),
            xaxis=dict(tickangle=-45),
        )

        chart = dcc.Graph(
            figure=fig,
            id={"type": "chart", "field": facet_name},
            config={"displayModeBar": False},
            style={"width": "100%", "height": "100%"},
        )

        # Add dropdown for all values (top 5 shown in chart, rest in dropdown)
        # For bar charts, include all values in dropdown so users can select from any
        all_options = [
            {"label": f"{item['value']} ({item['count']})", "value": item["value"]}
            for item in sorted_values
        ]
        dropdown = dcc.Dropdown(
            id={"type": "filter-dropdown", "field": facet_name},
            options=all_options,
            value=selected,
            multi=True,
            placeholder=f"Select {display_name}...",
            style={"marginTop": "5px", "fontSize": "0.85rem"},
        )

        if rest_values:
            return html.Div(
                [
                    chart,
                    html.P(
                        f"Showing top 5. Select from all {len(sorted_values)} values:",
                        style={
                            "fontSize": "0.75rem",
                            "color": "#666",
                            "marginTop": "5px",
                            "marginBottom": "3px",
                        },
                    ),
                    dropdown,
                ],
                style={"width": "100%", "marginBottom": "0.5rem"},
            )
        else:
            return html.Div(
                [chart, dropdown], style={"width": "100%", "marginBottom": "0.5rem"}
            )


# Helper function to format statistics for table display
def format_rnaseq_stats(result: Dict[str, Any]) -> str:
    """Format RNA-seq statistics for table display"""
    stats = []
    if "n_sig_perturb_pairs_up" in result:
        stats.append(f"Up: {format_value(result['n_sig_perturb_pairs_up'])}")
    if "n_sig_perturb_pairs_down" in result:
        stats.append(f"Down: {format_value(result['n_sig_perturb_pairs_down'])}")
    if "lfc_median" in result:
        stats.append(f"LFC median: {format_number(result['lfc_median'])}")
    return ", ".join(stats) if stats else "N/A"


def format_crispr_stats(result: Dict[str, Any]) -> str:
    """Format CRISPR Screen statistics for table display"""
    stats = []
    if "n_sig_crisps" in result:
        stats.append(f"Sig Crisps: {format_value(result['n_sig_crisps'])}")
    if "n_sig_crispr" in result:
        stats.append(f"Sig Crispr: {format_value(result['n_sig_crispr'])}")
    return ", ".join(stats) if stats else "N/A"


def format_mave_stats(result: Dict[str, Any]) -> str:
    """Format MAVE statistics for table display"""
    if "n_mave" in result:
        return f"N Mave: {format_value(result['n_mave'])}"
    return "N/A"
