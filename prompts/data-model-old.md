# Data model

## 1. Overview

This document describes model for the data which is being provided by BE and used by FE. It provides a high level contract on data exchange between BE and FE and does not include implementation details, which are provided in separate documents to BE and to FE only.

The data being represented can be semantically described in the following manner: "According to this DATASET, introducing this PERTURBATION leads to this CHANGE in this PHENOTYPE".

Hence, there are four individual object types being built by BE; and there are four corresponding columns in the FE interface.

## 2. API specification

API endpoint is /v1/search. It accepts the following optional parameters for filtering:

* dataset_metadata: string, used to free text search in the dataset metadata
* perturbation_gene_name: string, used to filter by a perturbed gene name
* change_direction: string, used to filter by the direction of the change. Must be "increased" or "decreased".
* phenotype_gene_name: string, used to filter by a phenotype gene name

In order to be returned, the data row must satisfy all specified filter conditions.

The API also accepts the following mandatory parameter for grouping:

* group_by: string, possible values are "perturbation_gene_name" or "phenotype_gene_name". This is used to aggregate the data rows.

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
* perturbation_gene_name, e.g. ABC123
* n_total, e.g. 120
* n_up, e.g. 20
* n_down, e.g. 100
n_total, n_up, n_down are integers. They show how many phenotypes in the dataset *are affected by* perturbing this gene.
Note that ALL these fields must be included into the perturbation object whether it's being grouped on, or whether it's a part of the data row in by_phenotype/perturbation_change/perturbation.

### 3.3. Change
* direction: string, either "increased" or "decreased" based on log2fc
* log2fc: double
* padj: double

### 3.4. Phenotype
* phenotype_gene_name, e.g. DEF456
* base_mean, float
* n_total, e.g. 5
* n_down, e.g. 1
* n_up, e.g. 4
n_total, n_up, n_down are integers. They show how many perturbations in the dataset *affect* this phenotype.
Note that ALL these fields must be included into the phenotype object whether it's being grouped on, or whether it's a part of the data row in by_perturbation/change_phenotype/phenotype.

## 4. Grouping

In addition, data is *always* aggregated by dataset on top level.

### 4.1. BE response if grouped by perturbation
```json
[
    {
        "dataset": {
            // dataset object fields
        },
        "by_perturbation": [
            {
                "perturbation": {
                    // perturbation object fields
                },
                "change_phenotype": [
                    {
                        "change": {
                            // change object fields
                        },
                        "phenotype": {
                            // phenotype object fields
                        }
                    }
                    // More change_phenotype elements
                ]
            }
            // More by_perturbation elements
        ]
    }
    // More top level dataset elements
]
```

### 4.2. BE response if aggregated by phenotype
```json
[
    {
        "dataset": {
            // dataset object fields
        },
        "by_phenotype": [
            {
                "phenotype": {
                    // phenotype object fields
                },
                "perturbation_change": [
                    {
                        "perturbation": {
                            // perturbation object fields
                        },
                        "change": {
                            // change object fields
                        }
                    }
                    // More perturbation_change elements
                ]
            }
            // More by_phenotype elements
        ]
    }
    // More top level dataset elements
]
```
