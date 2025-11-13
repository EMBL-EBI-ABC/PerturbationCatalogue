# Back-end Specification

First, study the updated data model specification in @prompts/data-model.md. This document provides the back-end implementation details for the new API structure.

## 1. Task Overview

Your job is to implement `be/data_query.py` to create the API endpoints defined in the data model. It will then be imported into the existing `be/main.py`. The implementation must be asynchronous.

*   `/v1/{MODALITY}/search`
*   `/v1/{MODALITY}/{DATASET_ID}/search`

Refer to @be/Dockerfile and @be/README.md for Docker setup and environment variables.

## 2. Data Sources

### 2.1. Databases
Data is split between Elastic (dataset metadata) and Postgres (perturbation/effect data). Use asynchronous libraries for all database interactions. Credentials are provided via environment variables.

### 2.2. Dataset Metadata (Elastic)

*   **Index**: `dataset-summary-v2`
*   **Modality Field**: `data_modalities` (string, e.g., "Perturb-seq")
    *   **Note on Casing**: The values in the `data_modalities` field are case-sensitive and may not match the lowercase modality names from the API URL. The implementation must handle this mapping. For example: `perturb-seq` -> `Perturb-seq`, `crispr-screen` -> `CRISPR screen`, `mave` -> `MAVE`.
*   **Query Fields**: Use the following fields for filtering and output, mapping them to the `dataset_*` prefixed fields in the API. For all `*_labels` fields, take the first element from the source list to produce a single string value in the response.
    *   `dataset_id` -> `dataset_id`
    *   `tissue_labels` -> `dataset_tissue`
    *   `cell_type_labels` -> `dataset_cell_type`
    *   `cell_line_labels` -> `dataset_cell_line`
    *   `sex_labels` -> `dataset_sex`
    *   `developmental_stage_labels` -> `dataset_developmental_stage`
    *   `disease_labels` -> `dataset_disease`
    *   `library_perturbation_type_labels` -> `dataset_library_perturbation_type`

### 2.3. Perturbation and Effect Data (Postgres)

Data for each modality is in a separate table within the `$PG_DB` database.

#### 2.3.1. `perturb-seq`
*   **Table**: `perturb_seq_2`
*   **Column Mapping**:
    *   `perturbed_target_symbol` -> `perturbation_gene_name`
    *   `gene` -> `effect_gene_name`
    *   `log2foldchange` -> `effect_log2fc` (used to derive `effect_direction`)
    *   `padj` -> `effect_padj`
    *   `basemean` -> `effect_base_mean`

#### 2.3.2. `crispr-screen`
*   **Table**: `crispr_data`
*   **Column Mapping**:
    *   `perturbed_target_symbol` -> `perturbation_gene_name`
    *   `score_name` -> `effect_score_name`
    *   `score_value` -> `effect_score_value`

#### 2.3.3. `mave`
*   **Table**: `mave_data`
*   **Column Mapping**:
    *   `perturbed_target_symbol` -> `perturbation_gene_name`
    *   `perturbation_name` -> `perturbation_name`
    *   `score_name` -> `effect_score_name`
    *   `score_value` -> `effect_score_value`

#### 2.3.4. `perturb-seq` Summary Views (Postgres)
To enrich the `perturb-seq` data with summary counts, two materialized views must be used:

*   **`perturb_seq_summary_perturbation`**: Provides counts for each perturbation.
    *   **Key**: `(dataset_id, perturbed_target_symbol)`
    *   **Columns**: `n_total`, `n_down`, `n_up` (map to `perturbation_n_total`, etc.)
*   **`perturb_seq_summary_effect`**: Provides counts for each effect/phenotype.
    *   **Key**: `(dataset_id, gene)`
    *   **Columns**: `n_total`, `n_down`, `n_up` (map to `effect_n_total`, etc.)

These lookups are essential for providing context to each `perturb-seq` result.

## 3. Query Logic

### 3.1. `/v1/{MODALITY}/search`

The query execution order is critical to ensure accurate facet counts.

1.  **Pre-filter Datasets (Postgres)**: First, execute a `SELECT DISTINCT dataset_id` query on the appropriate modality table.
    *   Apply all `perturbation_*` and `effect_*` filters from the API parameters. This will require parsing the numeric filter syntax (e.g., `0.5_2.0`) into the correct SQL conditions.
    *   This query returns the definitive list of `dataset_id`s that contain relevant data. If this list is empty, the search can terminate and return an empty response.

2.  **Filter Datasets & Get Facets (Elastic)**: Next, execute a single query against the `dataset-summary-v2` index.
    *   Filter by `data_modalities` based on the `{MODALITY}` from the URL.
    *   Apply all `dataset_*` filters from the API parameters, including the `dataset_search` free-text query.
    *   **Crucially, add a `terms` filter to restrict the query to the list of `dataset_id`s obtained from Postgres in the previous step.**
    *   Use an aggregation to compute `facet_counts` for all filterable `dataset_*` fields. The facets will now be correctly calculated based only on the datasets that have matching data in both databases.
    *   This query returns the final, correctly filtered list of datasets and the accurate `total_datasets_count`.

3.  **Paginate Datasets**: Apply `dataset_offset` and `dataset_limit` to the list of datasets returned from Elastic.

4.  **Fetch Data (Postgres)**: For each dataset in the paginated list:
    *   Construct a `SELECT` query on the appropriate modality table, filtering by `dataset_id`.
    *   Re-apply all `perturbation_*` and `effect_*` filters.
    *   Add `ORDER BY` clauses based on the `sort` parameter.
    *   Apply the `rows_per_dataset_limit`.
    *   **For `perturb-seq` only**: For each row returned, perform two additional fast lookups on the summary views (`perturb_seq_summary_perturbation` and `perturb_seq_summary_effect`) to fetch the `n_total`, `n_up`, and `n_down` counts.

5.  **Assemble Response**: Combine the dataset metadata from Elastic with the (potentially enriched) rows from Postgres to create the final response structure.

### 3.2. `/v1/{MODALITY}/{DATASET_ID}/search`

1.  **Count Rows (Postgres)**: Execute a `SELECT COUNT(*)` query on the appropriate modality table.
    *   Filter by `dataset_id` from the URL.
    *   Apply all `perturbation_*` and `effect_*` filters, parsing numeric syntax as needed.
    *   This provides the `total_rows_count`.
2.  **Fetch Rows (Postgres)**: Execute a `SELECT` query on the same table.
    *   Apply the same filters as the count query.
    *   Add `ORDER BY` clauses based on the `sort` parameter.
    *   Apply `LIMIT` and `OFFSET` for pagination.
    *   **For `perturb-seq` only**: Enrich each row with data from the summary views as described in section 3.1.
3.  **Assemble Response**: Format the results into the specified JSON structure.

## 4. Implementation Details

*   **Parameter Validation**: The API must reject any requests that include query parameters not explicitly defined in the data model. Return a `400 Bad Request`.
*   **Code Style**: Minimize code duplication. Write compact, succinct, and easy-to-understand code.
*   **Performance**: Ensure database queries are efficient. The initial `SELECT DISTINCT` query in Postgres must be fast. All sortable/filterable columns in Postgres should be indexed. The lookups on the summary views must be on their key columns to ensure high performance.