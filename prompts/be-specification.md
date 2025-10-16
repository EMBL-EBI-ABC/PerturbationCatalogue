# Back-end specification

First, study the data model specification in @prompts/data-model.md. It describes the API contract between the front-end and the back-end. This document provides the back-end implementation details.

## 1. Task overview

Study @be/Dockerfile and @be/README.md to see BE skeleton for the project. It describes how it's Dockerised and ran locally, as well as which environment variables it relies upon.

Your job is to implement be/main.py and add any corresponding requirements to be/requirements.txt.

## 2. Input data

### 2.1. Databases
There are two data sources: Elastic, which contains dataset metadata; and Postgres, which contains perturbation and effect data. All queries performed must be asynchronous for maximum performance.

Database access credentials are provided as environment variables; use them.

### 2.2. Dataset metadata in Elastic
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

### 2.3. Perturbation, change and phenotype data in Postgres
Inside the $PS_DB database, use a table named "perturb_seq". An example of data row in Postgres:

```tsv
dataset_id	perturbed_target_symbol	gene	log2foldchange	padj	basemean	max_ingested_at
adamson_2016_pilot	SPI1	CHST11	0.23117926187171009	0.50391124972854473	27.0810834104398	2025-09-30T13:37:59.073273
```

Use the API parameters to filter data from this table:
* "perturbation_gene_name" filters by "perturbed_target_symbol".
* "phenotype_gene_name" filters by "gene".
* "effect_direction" filters "log2foldchange" (> 0 for "increased", < 0 for "decreased").

As a reminder, all of the provided filters must match for the data piece to be returned.

### 2.4. Summary views in Postgres
The API requires grouping. In order to speed up performance, materialised views for such summary statistics are pre-computed.

When `group_by` is `perturbation_gene_name`, use the `perturb_seq_summary_perturbation` view:
```sql
## View perturb_seq_summary_perturbation
"dataset_id","perturbed_target_symbol","n_total","n_down","n_up"
"orion_2025_hek293t","MSRB1","1","1","0"
"orion_2025_hct116","SPRR4","10","4","6"
"orion_2025_hct116","NUB1","10","8","2"
"orion_2025_hek293t","GRK3","62","2","60"
"orion_2025_hct116","ID2","5","1","4"
```

When `group_by` is `effect_gene_name` (as specified in the data model), you should group by phenotype. Use the `perturb_seq_summary_effect` view for this, as it is grouped by `gene` (the phenotype gene):
```sql
## View perturb_seq_summary_effect
"dataset_id","gene","n_total","n_down","n_up","base_mean"
"adamson_2016_pilot","ABCB10","1","0","1","44.2341576206527"
"adamson_2016_pilot","AC005256.1","1","0","1","4.607358250167694"
"adamson_2016_pilot","AC006262.5","1","1","0","5.548977295822654"
"adamson_2016_pilot","AC096559.1","1","0","1","4.410971172950684"
"adamson_2016_pilot","AC108868.6","1","0","1","22.91355626193407"
```

## 3. Code style.
When implementing the BE, it's important to minimise code duplication; to make it as compact, succinct and easy to understand as possible.

Avoid adding too many comments or internal thoughts to the code - it should be compact, concise, simple and the point.

## 4. Query performance
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
