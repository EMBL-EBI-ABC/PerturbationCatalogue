# Data model

## 1. Overview

This document describes model for the data which is being provided by BE and used by FE. It provides a high level contract on data exchange between BE and FE and does not include implementation details, which are provided in separate documents to BE and to FE only.

The data being represented can be semantically described in the following manner: "According to this DATASET, introducing this PERTURBATION leads to this CHANGE in this PHENOTYPE".

The data is presented grouped by modality, and within each modality, by dataset. There are three modalities: "Perturb-seq", "CRISPR screen", and "MAVE".

## 2. API specification

API endpoint is /v1/search. It accepts the following optional parameters for filtering:

* `dataset_metadata`: string, used to free text search in the dataset metadata.
* `perturbation_gene_name`: string, used to filter by a perturbed gene name.
* `change_direction`: string, used to filter by the direction of the change. Must be "increased" or "decreased". Only applies to Perturb-seq modality. If any value is specified, only Perturb-seq data matching the filter should be returned.
* `phenotype_gene_name`: string, used to filter by a phenotype gene name.
* `modalities`: string, comma-separated list of modalities to return. Default is "perturb-seq,crispr-screen,mave".

In order to be returned, the data row must satisfy all specified filter conditions.

The API also accepts the following mandatory parameter for grouping:

* `group_by`: string, possible values are "perturbation_gene_name" or "phenotype_gene_name". This is used to aggregate the data rows for the "Perturb-seq" modality. For other modalities, no grouping is performed.

The API also accepts the following optional parameters for limiting the number of results:

* `max_datasets_per_modality`: integer, default 10. The maximum number of datasets to return for each modality.
* `max_top_level`: integer, default 10. For "Perturb-seq" modality, the maximum number of top-level groups (e.g., perturbations or phenotypes) to return per dataset.
* `max_rows`: integer, default 10. For "Perturb-seq" modality, the maximum number of underlying data rows (e.g., `change_phenotype` items) to return for each top-level group. For other modalities, the maximum number of data records to return per dataset.

## 3. Data types

### 3.1. Dataset
* dataset_id
* tissue
* cell_type
* cell_line
* sex
* developmental_stage
* disease
* library_perturbation_type
All fields are strings. dataset_id is used to link together Elastic and Postgres data. The rest of the fields are searchable metadata.

### 3.2. Perturbation

#### 3.2.1. For "Perturb-seq" modality
* perturbation_type: hardcoded to "gene_knockout"
* perturbation_gene_name, e.g. ABC123
* n_total, e.g. 120
* n_up, e.g. 20
* n_down, e.g. 100
n_total, n_up, n_down are integers. They show how many phenotypes in the dataset *are affected by* perturbing this gene.

#### 3.2.2. For "CRISPR screen" modality
* perturbation_type: hardcoded to "gene_knockout"
* perturbation_gene_name, e.g. ABC123

#### 3.2.3. For "MAVE" modality
* perturbation_type: hardcoded to "mave"
* perturbation_gene_name, e.g. ABC123
* perturbation_name, e.g. p.Pro73Gln

### 3.3. Change

#### 3.3.1. For "Perturb-seq" modality
* direction: string, either "increased" or "decreased" based on log2fc
* log2fc: double
* padj: double

#### 3.3.2. For "CRISPR screen" and "MAVE" modalities
* score_value: double

### 3.4. Phenotype

#### 3.4.1. For "Perturb-seq" modality
* phenotype_gene_name, e.g. DEF456
* base_mean, float
* n_total, e.g. 5
* n_down, e.g. 1
* n_up, e.g. 4
n_total, n_up, n_down are integers. They show how many perturbations in the dataset *affect* this phenotype.

#### 3.4.2. For "CRISPR screen" and "MAVE" modalities
* score_name: string

## 4. Response structure

The top-level response is an object containing a list of modalities and global facet counts for dataset metadata.

```json
{
    "modalities": [
        {
            "modality": "perturb-seq",
            "total_datasets_count": 123,
            "datasets": [
                // dataset objects for this modality, limited by max_datasets_per_modality
            ]
        },
        {
            "modality": "crispr-screen",
            "total_datasets_count": 45,
            "datasets": [
                // dataset objects for this modality
            ]
        },
        {
            "modality": "mave",
            "total_datasets_count": 67,
            "datasets": [
                // dataset objects for this modality
            ]
        }
    ],
    "facet_counts": {
        "tissue": [
            { "value": "blood", "count": 50 },
            { "value": "liver", "count": 20 }
        ],
        "cell_type": [
            { "value": "T cell", "count": 30 },
            { "value": "B cell", "count": 15 }
        ]
    }
}
```

### 4.1. For "Perturb-seq" modality (grouping applies)

The number of elements in `by_perturbation` or `by_phenotype` is limited by `max_top_level`. The number of elements in `change_phenotype` or `perturbation_change` is limited by `max_rows`.

#### 4.1.1. Grouped by perturbation
```json
{
    "dataset": {
        // dataset object fields
    },
    "by_perturbation": [
        {
            "perturbation": {
                // perturbation object fields for "Perturb-seq"
            },
            "change_phenotype": [
                {
                    "change": {
                        // change object fields for "Perturb-seq"
                    },
                    "phenotype": {
                        // phenotype object fields for "Perturb-seq"
                    }
                }
                // More change_phenotype elements
            ]
        }
        // More by_perturbation elements
    ]
}
```

#### 4.1.2. Grouped by phenotype
```json
{
    "dataset": {
        // dataset object fields
    },
    "by_phenotype": [
        {
            "phenotype": {
                // phenotype object fields for "Perturb-seq"
            },
            "perturbation_change": [
                {
                    "perturbation": {
                        // perturbation object fields for "Perturb-seq"
                    },
                    "change": {
                        // change object fields for "Perturb-seq"
                    }
                }
                // More perturbation_change elements
            ]
        }
        // More by_phenotype elements
    ]
}
```

### 4.2. For "CRISPR screen" and "MAVE" modalities (no grouping)
For these modalities, the structure is a simple list of records. The number of elements in `data` is limited by `max_rows`.

```json
{
    "dataset": {
        // dataset object fields
    },
    "data": [
        {
            "perturbation": {
                // perturbation object fields for the modality
            },
            "change": {
                // change object fields for the modality
            },
            "phenotype": {
                // phenotype object fields for the modality
            }
        }
        // More data records
    ]
}
```