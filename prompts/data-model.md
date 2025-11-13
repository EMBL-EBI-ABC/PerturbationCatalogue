# Data Model

## 1. Overview

This document describes the API contract for the Perturbation Catalogue. The data is represented using a `Dataset/Perturbation/Effect` model, and the API is organized by modality.

The core semantic concept is: "According to a `DATASET`, a `PERTURBATION` has an `EFFECT`."

The available modalities are: `perturb-seq`, `crispr-screen`, `mave`.

## 2. API Specification

### 2.1. Common Principles

#### 2.1.1. Endpoints
The API is structured around modality-specific endpoints:
*   `/v1/{MODALITY}/search`: Search across all datasets within a modality.
*   `/v1/{MODALITY}/{DATASET_ID}/search`: Search within a specific dataset in a modality.

#### 2.1.2. Parameter Prefixes
All API parameters and response fields are prefixed to indicate which part of the data model they correspond to:
*   `dataset_`: For dataset metadata.
*   `perturbation_`: For perturbation attributes.
*   `effect_`: For effect attributes (a combination of the old "change" and "phenotype" concepts).

#### 2.1.3. Filtering
*   **String fields**: Filtered by exact match (e.g., `perturbation_gene_name=SPI1`).
*   **Numeric (int/float) fields**: Filtered by a single parameter with a special syntax:
    *   `param=value`: Exact match (e.g., `effect_log2fc=1.5`).
    *   `param=min_max`: Inclusive range (e.g., `effect_log2fc=0.5_2.0`).
    *   `param=min_`: Greater than or equal to (e.g., `effect_log2fc=0.5_`).
    *   `param=_max`: Less than or equal to (e.g., `effect_log2fc=_2.0`).

#### 2.1.4. Sorting
Sorting is controlled by the `sort` parameter, which accepts a comma-separated list of fields. Each field can have a direction specifier (`:asc` or `:desc`).
*   **Syntax**: `sort=field1:direction,field2:direction`
*   **Example**: `sort=effect_padj:asc,perturbation_gene_name:desc`

### 2.2. Modality Search: `/v1/{MODALITY}/search`

This endpoint searches across multiple datasets for a given modality.

**Parameters:**
*   **Dataset Filters** (see 2.1.3):
    *   `dataset_search`: Free text search across dataset metadata.
    *   `dataset_tissue`, `dataset_cell_type`, `dataset_cell_line`, `dataset_sex`, `dataset_developmental_stage`, `dataset_disease`, `dataset_library_perturbation_type`.
*   **Perturbation/Effect Filters** (see 2.1.3):
    *   Any `perturbation_*` or `effect_*` field can be used as a filter.
*   **Sorting**:
    *   `sort`: Sorts the `results` within each dataset. Accepts any filterable `perturbation_*` or `effect_*` field.
*   **Pagination**:
    *   `dataset_limit`: Max number of datasets to return (default: 10).
    *   `dataset_offset`: Offset for dataset pagination (default: 0).
    *   `rows_per_dataset_limit`: Max number of (perturbation, effect) rows to return per dataset (default: 10).

### 2.3. Single Dataset Search: `/v1/{MODALITY}/{DATASET_ID}/search`

This endpoint searches for data within a single specified dataset.

**Parameters:**
*   **Perturbation/Effect Filters** (see 2.1.3):
    *   Same as the modality search endpoint.
*   **Sorting**:
    *   Same as the modality search endpoint.
*   **Pagination**:
    *   `limit`: Max number of rows to return (default: 50).
    *   `offset`: Offset for row pagination (default: 0).

## 3. Data Types

All fields are prefixed as described above.

### 3.1. Dataset
*   `dataset_id`: string
*   `dataset_tissue`: string
*   `dataset_cell_type`: string
*   `dataset_cell_line`: string
*   `dataset_sex`: string
*   `dataset_developmental_stage`: string
*   `dataset_disease`: string
*   `dataset_library_perturbation_type`: string

### 3.2. Perturbation
*   `perturbation_gene_name`: string (All modalities)
*   `perturbation_name`: string (MAVE only, e.g., "p.Pro73Gln")
*   **For `perturb-seq` only:**
    *   `perturbation_n_total`: integer (How many phenotypes are affected by this perturbation)
    *   `perturbation_n_up`: integer
    *   `perturbation_n_down`: integer

### 3.3. Effect
This entity combines attributes of the outcome of a perturbation.

*   **For `perturb-seq`:**
    *   `effect_gene_name`: string
    *   `effect_direction`: string ("increased" or "decreased")
    *   `effect_log2fc`: float
    *   `effect_padj`: float
    *   `effect_base_mean`: float
    *   `effect_n_total`: integer (How many perturbations affect this phenotype)
    *   `effect_n_up`: integer
    *   `effect_n_down`: integer
*   **For `crispr-screen` and `mave`:**
    *   `effect_score_name`: string
    *   `effect_score_value`: float

## 4. Response Structure

### 4.1. Modality Search: `/v1/{MODALITY}/search`

```json
{
    "total_datasets_count": 123,
    "facet_counts": {
        "dataset_tissue": [
            { "value": "blood", "count": 50 }
        ],
        "dataset_cell_type": [
            { "value": "T cell", "count": 30 }
        ]
    },
    "datasets": [
        {
            "dataset": {
                // dataset object fields
            },
            "results": [
                {
                    "perturbation": {
                        // perturbation object fields for this modality
                    },
                    "effect": {
                        // effect object fields for this modality
                    }
                }
                // ... more results, limited by rows_per_dataset_limit
            ]
        }
        // ... more datasets, limited by dataset_limit
    ]
}
```

### 4.2. Single Dataset Search: `/v1/{MODALITY}/{DATASET_ID}/search`

```json
{
    "total_rows_count": 456,
    "offset": 0,
    "limit": 50,
    "results": [
        {
            "perturbation": {
                // perturbation object fields
            },
            "effect": {
                // effect object fields
            }
        }
        // ... more results, limited by limit
    ]
}
```