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
BACKEND_URL = os.getenv("PERTURBATION_CATALOGUE_BE", "http://localhost:8000")

if not BACKEND_URL:
    raise ValueError("PERTURBATION_CATALOGUE_BE environment variable is not set.")

# Facet fields to exclude from stats display
FACET_FIELDS = [
    "license",
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


def fetch_all_search_results(
    query: Optional[str] = None,
    filters: Optional[Dict[str, List[str]]] = None,
    page_size: int = 100,
) -> List[Dict[str, Any]]:
    """Fetch all search results by paginating through the API.

    Args:
        query: Search query string
        filters: Filter dictionary
        page_size: Number of results per page (max 100)

    Returns:
        List of all result records
    """
    all_results = []
    page = 1

    while True:
        data = fetch_search_results(query=query, filters=filters, page=page, size=page_size)
        results = data.get("results", [])
        all_results.extend(results)

        total_pages = data.get("total_pages", 1)
        if page >= total_pages or not results:
            break
        page += 1

    return all_results


def fetch_modality_datasets(
    modality: str,
    filters: Optional[Dict[str, Any]] = None,
    dataset_limit: int = 4,
    dataset_offset: int = 0,
    rows_per_dataset_limit: Optional[int] = 5,
    sort: Optional[str] = None,
) -> Dict[str, Any]:
    """Call /v1/{modality}/search with the provided filters."""
    try:
        params: Dict[str, Any] = {
            "dataset_limit": dataset_limit,
            "dataset_offset": dataset_offset,
        }
        # Only add rows_per_dataset_limit if it's not None (MAVE doesn't use it)
        if rows_per_dataset_limit is not None:
            params["rows_per_dataset_limit"] = rows_per_dataset_limit
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
    limit: Optional[int] = 5,
    offset: Optional[int] = 0,
    sort: Optional[str] = None,
) -> Dict[str, Any]:
    """Call /v1/{modality}/{dataset_id}/search to fetch more rows."""
    try:
        params: Dict[str, Any] = {}
        # Only add limit and offset if they're not None (MAVE doesn't use them)
        if limit is not None:
            params["limit"] = limit
        if offset is not None:
            params["offset"] = offset
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


def fetch_dataset(dataset_id: str) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """Fetch dataset metadata from backend API."""
    try:
        response = requests.get(
            f"{BACKEND_URL}/dataset/{dataset_id}",
            timeout=10,
        )
        response.raise_for_status()
        return response.json(), None
    except Exception as exc:
        error_message = f"Error fetching dataset {dataset_id}: {exc}"
        print(error_message)
        return None, error_message


def fetch_perturb_seq_gsea(
    dataset_id: str,
    perturbed_gene_name: str,
) -> Dict[str, Any]:
    """Fetch GSEA results for a perturbed gene in a dataset."""
    try:
        params = {
            "dataset_id": dataset_id,
            "perturbed_gene_name": perturbed_gene_name,
        }
        response = requests.get(
            f"{BACKEND_URL}/v1/perturb-seq-gsea",
            params=params,
            timeout=30,
        )
        response.raise_for_status()
        return {"results": response.json(), "error": None}
    except Exception as exc:
        error_message = f"Error fetching GSEA data for {perturbed_gene_name}: {exc}"
        print(error_message)
        return {"results": [], "error": error_message}


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
