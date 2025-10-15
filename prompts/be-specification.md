# Back-end specification

## 1. Task overview

Study @be/Dockerfile and @be/README.md to see BE skeleton for the project. It describes how it's Dockerised and ran locally, as well as which environment variables it relies upon.

Your job is to implement be/main.py and add any corresponding requirements to be/requirements.txt.

## 2. High level API definition: /v1/search

The BE should have a single API on the URL "/v1/search". It accepts exactly four parameters, all optional:
* dataset_metadata: string, used to free text search in the dataset metadata
* perturbation_gene_name: string, used to filter by a perturbed gene name
* effect_gene_name: string, used to filter by an effect gene name
* group_by: string, possible values are "perturbation_gene_name" or "effect_gene_name"

The API retrieves information about datasets, perturbations and their effects. It then filters according to the first three parameters (applying the AND logic throughout), optionally groups by either perturbation or effect gene name, then compiles the output, and returns in JSON format.

## 3. Input data

### 3.1. Databases
There are two data sources: Elastic, which contains dataset metadata; and Postgres, which contains perturbation and effect data. All queries performed must be asynchronous for maximum performance.

Database access credentials are provided as environment variables; use them.

### 3.2. Dataset metadata in Elastic
Dataset metadata is stored in Elastic. The data core which you need to use is named dataset-summary-v1. This is an example of one row from it. It contains a lot of columns, but here only the relevant fields are kept:
```json
{
    "dataset_id": "datlinger_2017",
    "tissue_labels": [
      "blood"
    ],
    "cell_type_labels": [
      "T cell"
    ],
    "cell_line_labels": [
      "JURKAT cell"
    ],
    "sex_labels": [
      "male"
    ],
    "developmental_stage_labels": [
      "adolescent"
    ],
    "disease_labels": [
      "T-cell acute lymphoblastic leukemia"
    ],
    "library_perturbation_type_labels": [
      "knockout"
    ]
}
```

You must only use a fixed subset of fields both for search and for returning the results. Here is the subset: dataset_id, tissue_labels, cell_type_labels, cell_line_labels, sex_labels, developmental_stage_labels, disease_labels, library_perturbation_type_labels.

You will see that some of these fields are lists and are named in plural (tissue_labels). You must flatten the list, always assume it only contains a single entry (make an assert), and name output in singular, for example: tissue_label.

### 3.3. Perturbation and effect data in Postgres
Inside the $PS_DB database, use a table named "perturb_seq". An example of data row in Postgres:

```tsv
dataset_id	perturbed_target_symbol	gene	log2foldchange	padj	basemean	max_ingested_at
adamson_2016_pilot	SPI1	CHST11	0.23117926187171009	0.50391124972854473	27.0810834104398	2025-09-30T13:37:59.073273
```

Use "perturbation_gene_name" API parameter to filter by perturbed_target_symbol, and "effect_gene_name" parameter to filter by gene. As a reminder, all of the provided filters must match for the data piece to be returned.

### 3.4. Perturbation and effect summary views in Postgres
When group_by is set then, as explained below, we want to extract certain summary statistics of perturbations or effects, and sort data by those statistics.

In order to speed up performance, materialised views for such summary statistics are pre-computed:
```sql
## View perturb_seq_summary_perturbation - when grouping by perturbation
"dataset_id","perturbed_target_symbol","n_total","n_down","n_up"
"orion_2025_hek293t","MSRB1","1","1","0"
"orion_2025_hct116","SPRR4","10","4","6"
"orion_2025_hct116","NUB1","10","8","2"
"orion_2025_hek293t","GRK3","62","2","60"
"orion_2025_hct116","ID2","5","1","4"

## View perturb_seq_summary_effect - when grouping by effect
"dataset_id","gene","n_total","n_down","n_up","base_mean"
"adamson_2016_pilot","ABCB10","1","0","1","44.2341576206527"
"adamson_2016_pilot","AC005256.1","1","0","1","4.607358250167694"
"adamson_2016_pilot","AC006262.5","1","1","0","5.548977295822654"
"adamson_2016_pilot","AC096559.1","1","0","1","4.410971172950684"
"adamson_2016_pilot","AC108868.6","1","0","1","22.91355626193407"
```

## 4. Output format

### 4.1. Top level
Top level of the output is always a single list. Each element of the list is information related to one particular dataset. At this level, *all* datasets which remain after filtering are returned:
```json
{"dataset_id": "some_dataset_123", "tissue_label": "foo", "cell_type_label": "bar", "cell_line_labels": "baz", "sex_label": "female", "developmental_stage_label": "Adult", "disease_label": "Healthy", "library_perturbation_type_label": "knockout", "data": []}
```

The `data` field is always a list. It must never be empty: only datasets with some data found for them must be returned.

Structure of the elements of the `data` list depends on the value of the group_by field and is described in the next section.

### 4.2. `data` field elements if no grouping is specified
If no grouping is specified, we simply request the *top 20* data rows for this given dataset, subject to perturbation and effect gene filters (if specififed), with the smallest padj. Each element is a flat data point:
```json
{"perturbation_gene_name": "ABC123", "effect_gene_name": "DEF456", "base_mean": 0.13, "log2foldchange": -2.17, "padj": 3.2e-17}
```

### 4.3. `data` field elements if group_by = "perturbation_gene_name"
In this case, we want to group results - within each dataset - by the perturbed gene name (named "perturbed_target_symbol" in Postgres). Subject to all filters, we want to display the *top 3* perturbed genes, sorting by n_total descending using the perturb_seq_summary_perturbation view. Within each of those 3 perturbed genes, we want to display the *top 10* perturbation effects, sorted by padj increasing. Each element in the `data` field looks like this:
```json
{"perturbation_gene_name": "ABC123", "n_total": 1100, "n_up": 400, "n_down": 700, "effects": [...]}
```
Where each `effect` element looks like this:
```json
{"effect_gene_name": "DEF456", "base_mean": 0.13, "log2foldchange": -2.17, "padj": 3.2e-17}
```

### 4.4. `data` field elements if group_by = "effect_gene_name"
This is similar to the previous case, but we are grouping by effect gene rather than a perturbed gene.

Subject to all filters, we want to display the *top 3* effect genes, sorting by n_total descending using the perturb_seq_summary_effect view. Within each of those 3 effect genes, we want to display the *top 10* perturbations which affected them the most, sorting by padj increasing. Each element in the `data` field looks like this:
```json
{"effect_gene_name": "DEF456", "n_total": 1100, "n_up": 400, "n_down": 700, "base_mean": 213.4, "perturbations": [...]}
```
Where each `perturbation` element looks like this:
```json
{"perturbation_gene_name": "DEF456", "log2foldchange": -2.17, "padj": 3.2e-17}
```

## 5. Avoiding code duplication
When implementing the BE, it's important to minimise code duplication; to make it as compact, succinct and easy to understand as possible.

## 6. Query performance
When constructing the queries, it's extremely important to make sure they are efficient. You must only ever use the quick look up operations, which are:
* All lookups from Elastic
* Lookups on the main perturb_seq table which use the following indexes:
  ```
  "indexname","indexdef"
  "idx_dataset_padj","CREATE INDEX idx_dataset_padj ON public.perturb_seq USING btree (dataset_id, padj)"
  "idx_perturbed_padj","CREATE INDEX idx_perturbed_padj ON public.perturb_seq USING btree (perturbed_target_symbol, padj)"
  "idx_gene_padj","CREATE INDEX idx_gene_padj ON public.perturb_seq USING btree (gene, padj)"
  "idx_perturbed_gene_padj","CREATE INDEX idx_perturbed_gene_padj ON public.perturb_seq USING btree (perturbed_target_symbol, gene, padj)"
  ```
* All lookups on the summary views are fast.
