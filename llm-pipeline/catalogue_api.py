import requests
import pandas as pd
import numpy as np
from scipy import stats
import logging
import time

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)

BASE_URL = "https://perturbation-catalogue-be-328296435987.europe-west2.run.app"

# These are the score names we recognise as primary effect scores
# Order matters — we prefer the first match found
EFFECT_SCORE_NAMES = [
    "CRISPR Score (CS)",
    "Log2FC",
    "Gamma (normalized log2e/t)",
    "Rho (Log2e Treated vs. Untreated)",
    "MAGeCK neg score",
    "LFC",
]

FDR_SCORE_NAMES = [
    "FDR",
    "fdr",
    "q-value",
    "adjusted p-value",
]


def query_crispr_screen(dataset_id=None, limit=100, max_records=5000):
    """
    Query CRISPR screen data from the Perturbation Catalogue API.

    Parameters
    ----------
    dataset_id : str or None
        If provided, query a specific dataset (e.g. "biogrid_5").
        If None, query across all CRISPR screen datasets.
    limit : int
        Records per page. 100 is a reasonable default.
    max_records : int
        Maximum total records to retrieve. Prevents runaway queries
        on datasets with hundreds of thousands of rows.

    Returns
    -------
    list of raw result dicts from the API
    """
    if dataset_id:
        endpoint = f"{BASE_URL}/v1/crispr-screen/{dataset_id}/search"
    else:
        endpoint = f"{BASE_URL}/v1/crispr-screen/search"

    all_results = []
    offset = 0

    while True:
        params = {"limit": limit, "offset": offset}

        try:
            response = requests.get(endpoint, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
        except requests.exceptions.RequestException as e:
            log.error(f"API request failed: {e}")
            break

        # Handle both single dataset response and multi-dataset response
        # Single dataset: {"total_rows_count": N, "results": [...]}
        # Multi dataset: [{"dataset": {...}, "results": [...]}, ...]
        if isinstance(data, list):
            # Multi-dataset response
            for dataset_block in data:
                results = dataset_block.get("results", [])
                dataset_meta = dataset_block.get("dataset", {})
                for r in results:
                    r["_dataset_meta"] = dataset_meta
                all_results.extend(results)
            break
        elif isinstance(data, dict):
            results = data.get("results", [])
            all_results.extend(results)
            total = data.get("total_rows_count", 0)
            if offset + limit >= total:
                break
        else:
            break

        offset += limit
        if len(all_results) >= max_records:
            log.info(f"Reached max_records limit ({max_records})")
            break

        time.sleep(0.1)

    all_results = all_results[:max_records]
    log.info(f"Retrieved {len(all_results)} raw records")
    return all_results


def query_perturb_seq(dataset_id=None, limit=100, max_records=5000):
    """
    Query scPerturb-seq data from the Perturbation Catalogue API.

    Same structure as query_crispr_screen but hits the perturb-seq
    endpoint. The response format is identical — perturbation + effect
    per record.
    """
    if dataset_id:
        endpoint = f"{BASE_URL}/v1/perturb-seq/{dataset_id}/search"
    else:
        endpoint = f"{BASE_URL}/v1/perturb-seq/search"

    all_results = []
    offset = 0

    while True:
        params = {"limit": limit, "offset": offset}

        try:
            response = requests.get(endpoint, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
        except requests.exceptions.RequestException as e:
            log.error(f"API request failed: {e}")
            break

        if isinstance(data, list):
            for dataset_block in data:
                results = dataset_block.get("results", [])
                dataset_meta = dataset_block.get("dataset", {})
                for r in results:
                    r["_dataset_meta"] = dataset_meta
                all_results.extend(results)
            break
        elif isinstance(data, dict):
            results = data.get("results", [])
            all_results.extend(results)
            total = data.get("total_rows_count", 0)
            if offset + limit >= total:
                break
        else:
            break

        offset += limit
        if len(all_results) >= max_records:
            break

        time.sleep(0.1)

    all_results = all_results[:max_records]
    log.info(f"Retrieved {len(all_results)} raw records")
    return all_results


def query_mave(dataset_id=None, limit=100, max_records=5000):
    """
    Query MAVE data from the Perturbation Catalogue API.

    MAVE data has a different biological meaning but the same API
    structure. Instead of fitness scores, values represent variant
    effect scores — how much a specific mutation changes protein
    function.
    """
    if dataset_id:
        endpoint = f"{BASE_URL}/v1/mave/{dataset_id}/search"
    else:
        endpoint = f"{BASE_URL}/v1/mave/search"

    all_results = []
    offset = 0

    while True:
        params = {"limit": limit, "offset": offset}

        try:
            response = requests.get(endpoint, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
        except requests.exceptions.RequestException as e:
            log.error(f"API request failed: {e}")
            break

        if isinstance(data, list):
            for dataset_block in data:
                results = dataset_block.get("results", [])
                dataset_meta = dataset_block.get("dataset", {})
                for r in results:
                    r["_dataset_meta"] = dataset_meta
                all_results.extend(results)
            break
        elif isinstance(data, dict):
            results = data.get("results", [])
            all_results.extend(results)
            total = data.get("total_rows_count", 0)
            if offset + limit >= total:
                break
        else:
            break

        offset += limit
        if len(all_results) >= max_records:
            break

        time.sleep(0.1)

    all_results = all_results[:max_records]
    log.info(f"Retrieved {len(all_results)} raw records")
    return all_results


def identify_primary_score(score_names_in_dataset):
    """
    Identify which score name is the primary effect score for a dataset.

    Different datasets use different score names. We check against
    our known list in priority order and return the first match.

    If nothing matches we fall back to whatever score name appears
    most frequently — a reasonable heuristic.

    Parameters
    ----------
    score_names_in_dataset : list of str
        All unique score names found in this dataset.

    Returns
    -------
    str — the score name to use as primary effect score
    """
    for known_score in EFFECT_SCORE_NAMES:
        if known_score in score_names_in_dataset:
            return known_score

    # Fallback — use most common score name that isn't FDR-like
    non_fdr = [
        s
        for s in score_names_in_dataset
        if not any(fdr in s.lower() for fdr in ["fdr", "q-value", "p-value"])
    ]
    if non_fdr:
        return non_fdr[0]

    return score_names_in_dataset[0]


def pivot_gene_records(raw_results):
    """
    Pivot API records from (gene, score_type) rows to one row per gene.

    This handles the record structure problem — the API returns
    multiple rows per gene, one for each score type reported.
    We collapse these into one clean record per gene.

    Parameters
    ----------
    raw_results : list of dicts
        Raw API response records with perturbation and effect fields.

    Returns
    -------
    pd.DataFrame with one row per gene, columns for each score type.
    """
    if not raw_results:
        return pd.DataFrame()

    rows = []
    for r in raw_results:
        perturbation = r.get("perturbation", {})
        effect = r.get("effect", {})
        dataset_meta = r.get("_dataset_meta", {})

        cell_lines = dataset_meta.get("dataset_cell_lines", [])
        if not cell_lines and dataset_meta:
            log.warning(
                "No cell line label found in dataset metadata — using 'unknown'"
            )
        cell_line_value = cell_lines[0] if cell_lines else "unknown"

        rows.append(
            {
                "gene": perturbation.get("gene_name", "unknown"),
                "score_name": effect.get("score_name", "unknown"),
                "score_value": effect.get("score_value", np.nan),
                "significant": effect.get("significant", "False") == "True",
                "significance_criteria": effect.get("significance_criteria", ""),
                "dataset_id": dataset_meta.get("dataset_id", "unknown"),
                "cell_line": cell_line_value,
                "disease": (
                    dataset_meta.get("dataset_diseases", ["unknown"])[0]
                    if dataset_meta.get("dataset_diseases")
                    else "unknown"
                ),
                "perturbation_type": (
                    dataset_meta.get("dataset_perturbation_types", ["unknown"])[0]
                    if dataset_meta.get("dataset_perturbation_types")
                    else "unknown"
                ),
            }
        )

    df = pd.DataFrame(rows)

    if df.empty:
        return df

    # Find the primary effect score for this batch of records
    score_names = df["score_name"].unique().tolist()
    primary_score = identify_primary_score(score_names)
    fdr_score = next(
        (s for s in score_names if any(f in s.lower() for f in ["fdr", "q-value"])),
        None,
    )

    log.info(f"Primary effect score identified: '{primary_score}'")
    if fdr_score:
        log.info(f"FDR score identified: '{fdr_score}'")

    # Pivot — one row per gene
    # Get effect scores
    effect_df = df[df["score_name"] == primary_score][
        [
            "gene",
            "score_value",
            "significant",
            "dataset_id",
            "cell_line",
            "disease",
            "perturbation_type",
        ]
    ].rename(columns={"score_value": "effect_score"})

    # Get FDR scores if available
    if fdr_score:
        fdr_df = df[df["score_name"] == fdr_score][["gene", "score_value"]].rename(
            columns={"score_value": "fdr"}
        )
        result = effect_df.merge(fdr_df, on="gene", how="left")
    else:
        result = effect_df.copy()
        result["fdr"] = np.nan

    # Drop duplicate genes — keep the one with highest absolute effect
    result = result.copy()
    result["abs_effect"] = result["effect_score"].abs()
    result = result.sort_values("abs_effect", ascending=False)
    result = result.drop_duplicates(subset=["gene"], keep="first")
    result = result.drop(columns=["abs_effect"])

    log.info(f"Pivoted to {len(result)} unique genes")
    return result


def get_dataset_metadata(dataset_id):
    """
    Fetch dataset-level metadata from the Catalogue API.

    When querying a specific dataset by ID, the API returns only
    gene-level results without dataset metadata. This function
    makes a separate call to get cell line, disease, and other
    context that enriches the training records.

    Parameters
    ----------
    dataset_id : str
        Catalogue dataset ID. e.g. "biogrid_2373"

    Returns
    -------
    dict with clean metadata fields
    """
    endpoint = f"{BASE_URL}/dataset/{dataset_id}"

    try:
        response = requests.get(endpoint, timeout=30)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        log.warning(f"Could not fetch metadata for {dataset_id}: {e}")
        return {}

    # Extract the fields we care about
    # Each is a list — take first element if available
    def first(lst):
        return lst[0] if lst else "unknown"

    return {
        "cell_line": first(data.get("cell_line_labels", [])),
        "disease": first(data.get("disease_labels", [])),
        "tissue": first(data.get("tissue_labels", [])),
        "cell_type": first(data.get("cell_type_labels", [])),
        "perturbation_type": first(data.get("perturbation_type_labels", [])),
        "treatment": first(data.get("treatment_labels", [])),
        "score_interpretation": data.get("score_interpretation", ""),
        "experiment_summary": data.get("experiment_summary", ""),
        "library_perturbation_type": first(data.get("library_perturbation_type_labels", [])),
    }
