# Data model

## 1. Overview

This document describes model for the data which is being provided by BE and used by FE. It provides a high level contract on data exchange between BE and FE and does not include implementation details, which are provided in separate documents to BE and to FE only.

The data being represented can be semantically described in the following manner: "According to this DATASET, introducing this PERTURBATION leads to this CHANGE in this PHENOTYPE".

The data is presented grouped by modality, and within each modality, by dataset. There are three modalities: "Perturb-seq", "CRISPR screen", and "MAVE".

## 2. API specification

API endpoint is /v1/search. It accepts the following optional parameters for filtering:

* dataset_metadata: string, used to free text search in the dataset metadata
* perturbation_gene_name: string, used to filter by a perturbed gene name
* change_direction: string, used to filter by the direction of the change. Must be "increased" or "decreased".
* phenotype_gene_name: string, used to filter by a phenotype gene name

In order to be returned, the data row must satisfy all specified filter conditions.

The API also accepts the following mandatory parameter for grouping:

* group_by: string, possible values are "perturbation_gene_name" or "phenotype_gene_name". This is used to aggregate the data rows for the "Perturb-seq" modality. For other modalities, no grouping is performed.

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
* phenotype_gene_name, e.g. DEF456
* perturbation_name, e.g. V123D

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

The top-level response is a list of modalities.

```json
[
    {
        "modality": "perturb-seq",
        "datasets": [
            // dataset objects for this modality
        ]
    },
    {
        "modality": "crispr-screen",
        "datasets": [
            // dataset objects for this modality
        ]
    },
    {
        "modality": "mave",
        "datasets": [
            // dataset objects for this modality
        ]
    }
]
```

### 4.1. For "Perturb-seq" modality (grouping applies)

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
For these modalities, the structure is a simple list of records.

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