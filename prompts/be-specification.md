First, study general instructions for LLM agents in @prompts/README.md.

Then study all files in the @be directory which contains back-end implementation.

In addition to existing APIs (keep them and do not modify them whatsoever), I want you to implement a new one.

It should use a URL "/v1/search" and should accept exactly three parameters, all optional:
* dataset_metadata: string, used to free text search in the dataset metadata
* perturbation_gene_name: string, used to filter by a perturbed gene name
* effect_gene_name: string, used to filter by an effect gene name

What it does is retrieves information from databases according to the three filters specified, compiles it, and returns data in JSON format.

# Environment variables
The following environment variables are set in the back-end. These are for your reference, so that you know which ones you can use in the code.
## Back-end
export ES_URL=https://....es.europe-west2.gcp.elastic-cloud.com
export ES_USERNAME=elastic
export ES_PASSWORD=...
## Data warehousing
export ELASTIC_ENDPOINT="https://elastic:...@....es.europe-west2.gcp.elastic-cloud.com"
## PostgreSQL
export PS_HOST=...
export PS_PORT=...
export PS_USER=...
export PS_PASSWORD=...
export PS_DB=...

# Overall data organisation
Dataset metadata is stored in Elastic. Perturbation and effect data are stored in Postgres. The API will need to query both databases asynchronously, stitch the results together using the dataset_id field, and reply in JSON format.

# Data in Elastic
Data core which you need to use from Elastic is named dataset-summary-v1. This is an example of one row from it:
```json
{
  "_index": "dataset-summary-v1",
  "_id": "datlinger_2017",
  "_version": 1,
  "_source": {
    "dataset_id": "datlinger_2017",
    "timepoints": [
      "P0DT0H0M0S"
    ],
    "treatment_labels": [
      "anti-CD3 antibody|anti-CD28 antibody",
      "UNTREATED CONTROL"
    ],
    "treatment_ids": [
      "EFO:0003317|EFO:0003304",
      "NCIT:C184729"
    ],
    "model_system_labels": [
      "cell line"
    ],
    "model_system_ids": [],
    "tissue_labels": [
      "blood"
    ],
    "tissue_ids": [
      "UBERON:0000178"
    ],
    "cell_type_labels": [
      "T cell"
    ],
    "cell_type_ids": [
      "CL:0000084"
    ],
    "cell_line_labels": [
      "JURKAT cell"
    ],
    "cell_line_ids": [
      "CLO:0007043"
    ],
    "sex_labels": [
      "male"
    ],
    "sex_ids": [],
    "developmental_stage_labels": [
      "adolescent"
    ],
    "developmental_stage_ids": [],
    "disease_labels": [
      "T-cell acute lymphoblastic leukemia"
    ],
    "disease_ids": [
      "MONDO:0004963"
    ],
    "study_title": "Pooled CRISPR screening with single-cell transcriptome readout",
    "study_uri": "https://doi.org/10.1038/nmeth.4177",
    "study_year": 2017,
    "first_author": "Paul Datlinger",
    "last_author": "Christoph Bock",
    "experiment_title": "Transcriptomics measurements of 5905 Jurkat cells induced with anti-CD3 and anti-CD28 antibodies",
    "experiment_summary": "\n            Jurkat cells were transduced with a gRNA library targeting high-level\n            regulators of T cell receptor signaling and a set of transcription factors. After 10\n            days of antibiotic selection and expansion, cells were stimulated with anti-CD3 and\n            anti-CD28 antibodies or left untreated. Both conditions were analyzed using CROP-seq,\n            measuring TCR activation for each gene knockout. The dataset comprises 5,905 high-quality\n            single-cell transcriptomes with uniquely assigned gRNAs.\n            ",
    "number_of_perturbed_targets": [
      33
    ],
    "number_of_perturbed_samples": [
      5905
    ],
    "library_generation_type_labels": [
      "endogenous"
    ],
    "library_generation_type_ids": [
      "EFO:0022868"
    ],
    "library_generation_method_labels": [
      "SpCas9"
    ],
    "library_generation_method_ids": [
      "EFO:0022876"
    ],
    "enzyme_generation_delivery_labels": [
      "lentiviral transduction"
    ],
    "enzyme_generation_delivery_ids": [],
    "library_generation_delivery_labels": [
      "lentiviral transduction"
    ],
    "library_generation_delivery_ids": [],
    "enzyme_integration_state_labels": [
      "random locus integration"
    ],
    "enzyme_integration_state_ids": [],
    "library_integration_state_labels": [
      "random locus integration"
    ],
    "library_integration_state_ids": [],
    "enzyme_expression_control_labels": [
      "constitutive expression"
    ],
    "enzyme_expression_control_ids": [],
    "library_expression_control_labels": [
      "constitutive expression"
    ],
    "library_expression_control_ids": [],
    "library_names": [
      "custom"
    ],
    "library_uris": [],
    "library_format_labels": [
      "pooled"
    ],
    "library_format_ids": [],
    "library_scope_labels": [
      "focused"
    ],
    "library_scope_ids": [],
    "library_perturbation_type_labels": [
      "knockout"
    ],
    "library_perturbation_type_ids": [],
    "library_manufacturers": [
      "Bock"
    ],
    "library_lentivaral_generations": [
      "3"
    ],
    "library_grnas_per_targets": [
      "3"
    ],
    "library_total_grnas": [
      116
    ],
    "library_total_variants": [],
    "readout_dimensionality_labels": [
      "high-dimensional assay"
    ],
    "readout_dimensionality_ids": [],
    "readout_type_labels": [
      "transcriptomic"
    ],
    "readout_type_ids": [],
    "readout_technology_labels": [
      "single-cell rna-seq"
    ],
    "readout_technology_ids": [],
    "method_name_labels": [
      "Perturb-seq"
    ],
    "method_name_ids": [],
    "method_uri": [],
    "sequencing_library_kit_labels": [
      "Nextera XT"
    ],
    "sequencing_library_ki_ids": [],
    "sequencing_platform_labels": [
      "Illumina HiSeq 4000"
    ],
    "sequencing_platform_ids": [],
    "sequencing_strategy_labels": [
      "barcode sequencing"
    ],
    "sequencing_strategy_ids": [],
    "software_counts_labels": [
      "Drop-seq Tools"
    ],
    "software_counts_ids": [],
    "software_analysis_labels": [
      "Custom"
    ],
    "software_analysis_ids": [],
    "reference_genome_labels": [
      "GRCh38"
    ],
    "reference_genome_ids": [],
    "associated_datasets": [
      "[{\"dataset_accession\": \"GSE92872\", \"dataset_uri\": \"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92872\", \"dataset_description\": \"Digital expression matrix\", \"dataset_file_name\": \"GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz\"}, {\"dataset_accession\": \"DatlingerBock2017.h5ad\", \"dataset_uri\": \"https://zenodo.org/records/13350497/files/DatlingerBock2017.h5ad\", \"dataset_description\": \"Processed .h5ad file\", \"dataset_file_name\": \"DatlingerBock2017.h5ad\"}]"
    ]
  },
  "fields": {
    "method_name_labels.text": [
      "Perturb-seq"
    ],
    "experiment_title": [
      "Transcriptomics measurements of 5905 Jurkat cells induced with anti-CD3 and anti-CD28 antibodies"
    ],
    "software_counts_labels": [
      "drop-seq tools"
    ],
    "library_names": [
      "custom"
    ],
    "library_perturbation_type_labels.text": [
      "knockout"
    ],
    "associated_datasets": [
      "[{\"dataset_accession\": \"gse92872\", \"dataset_uri\": \"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse92872\", \"dataset_description\": \"digital expression matrix\", \"dataset_file_name\": \"gse92872_crop-seq_jurkat_tcr.digital_expression.csv.gz\"}, {\"dataset_accession\": \"datlingerbock2017.h5ad\", \"dataset_uri\": \"https://zenodo.org/records/13350497/files/datlingerbock2017.h5ad\", \"dataset_description\": \"processed .h5ad file\", \"dataset_file_name\": \"datlingerbock2017.h5ad\"}]"
    ],
    "enzyme_expression_control_labels.text": [
      "constitutive expression"
    ],
    "sequencing_library_kit_labels.text": [
      "Nextera XT"
    ],
    "timepoints": [
      "p0dt0h0m0s"
    ],
    "study_uri": [
      "https://doi.org/10.1038/nmeth.4177"
    ],
    "tissue_ids": [
      "uberon:0000178"
    ],
    "library_generation_delivery_labels.text": [
      "lentiviral transduction"
    ],
    "readout_technology_labels": [
      "single-cell rna-seq"
    ],
    "cell_line_labels": [
      "jurkat cell"
    ],
    "enzyme_generation_delivery_labels.text": [
      "lentiviral transduction"
    ],
    "library_generation_method_labels.text": [
      "SpCas9"
    ],
    "software_counts_labels.text": [
      "Drop-seq Tools"
    ],
    "method_name_labels": [
      "perturb-seq"
    ],
    "study_title": [
      "Pooled CRISPR screening with single-cell transcriptome readout"
    ],
    "sequencing_platform_labels": [
      "illumina hiseq 4000"
    ],
    "cell_type_ids": [
      "cl:0000084"
    ],
    "library_expression_control_labels": [
      "constitutive expression"
    ],
    "library_generation_delivery_labels": [
      "lentiviral transduction"
    ],
    "dataset_id": [
      "datlinger_2017"
    ],
    "library_integration_state_labels": [
      "random locus integration"
    ],
    "number_of_perturbed_samples": [
      5905
    ],
    "reference_genome_labels": [
      "grch38"
    ],
    "cell_type_labels.text": [
      "T cell"
    ],
    "experiment_title.raw": [
      "transcriptomics measurements of 5905 jurkat cells induced with anti-cd3 and anti-cd28 antibodies"
    ],
    "number_of_perturbed_targets": [
      33
    ],
    "tissue_labels.text": [
      "blood"
    ],
    "library_manufacturers": [
      "bock"
    ],
    "readout_type_labels": [
      "transcriptomic"
    ],
    "sequencing_strategy_labels.text": [
      "barcode sequencing"
    ],
    "readout_technology_labels.text": [
      "single-cell rna-seq"
    ],
    "model_system_labels.text": [
      "cell line"
    ],
    "experiment_summary": [
      "\n            Jurkat cells were transduced with a gRNA library targeting high-level\n            regulators of T cell receptor signaling and a set of transcription factors. After 10\n            days of antibiotic selection and expansion, cells were stimulated with anti-CD3 and\n            anti-CD28 antibodies or left untreated. Both conditions were analyzed using CROP-seq,\n            measuring TCR activation for each gene knockout. The dataset comprises 5,905 high-quality\n            single-cell transcriptomes with uniquely assigned gRNAs.\n            "
    ],
    "disease_labels": [
      "t-cell acute lymphoblastic leukemia"
    ],
    "enzyme_integration_state_labels.text": [
      "random locus integration"
    ],
    "sex_labels": [
      "male"
    ],
    "study_year": [
      2017
    ],
    "developmental_stage_labels.text": [
      "adolescent"
    ],
    "library_names.text": [
      "custom"
    ],
    "sex_labels.text": [
      "male"
    ],
    "library_format_labels": [
      "pooled"
    ],
    "readout_dimensionality_labels": [
      "high-dimensional assay"
    ],
    "library_total_grnas": [
      116
    ],
    "readout_dimensionality_labels.text": [
      "high-dimensional assay"
    ],
    "disease_ids": [
      "mondo:0004963"
    ],
    "tissue_labels": [
      "blood"
    ],
    "study_title.raw": [
      "pooled crispr screening with single-cell transcriptome readout"
    ],
    "library_generation_method_labels": [
      "spcas9"
    ],
    "library_generation_type_labels": [
      "endogenous"
    ],
    "software_analysis_labels.text": [
      "Custom"
    ],
    "enzyme_generation_delivery_labels": [
      "lentiviral transduction"
    ],
    "first_author.text": [
      "Paul Datlinger"
    ],
    "cell_line_ids": [
      "clo:0007043"
    ],
    "first_author": [
      "paul datlinger"
    ],
    "readout_type_labels.text": [
      "transcriptomic"
    ],
    "library_format_labels.text": [
      "pooled"
    ],
    "library_lentivaral_generations": [
      "3"
    ],
    "treatment_labels": [
      "anti-cd3 antibody|anti-cd28 antibody",
      "untreated control"
    ],
    "cell_type_labels": [
      "t cell"
    ],
    "reference_genome_labels.text": [
      "GRCh38"
    ],
    "library_generation_type_ids": [
      "efo:0022868"
    ],
    "library_generation_type_labels.text": [
      "endogenous"
    ],
    "enzyme_integration_state_labels": [
      "random locus integration"
    ],
    "library_expression_control_labels.text": [
      "constitutive expression"
    ],
    "software_analysis_labels": [
      "custom"
    ],
    "last_author": [
      "christoph bock"
    ],
    "library_perturbation_type_labels": [
      "knockout"
    ],
    "treatment_labels.text": [
      "anti-CD3 antibody|anti-CD28 antibody",
      "UNTREATED CONTROL"
    ],
    "library_generation_method_ids": [
      "efo:0022876"
    ],
    "disease_labels.text": [
      "T-cell acute lymphoblastic leukemia"
    ],
    "last_author.text": [
      "Christoph Bock"
    ],
    "cell_line_labels.text": [
      "JURKAT cell"
    ],
    "model_system_labels": [
      "cell line"
    ],
    "library_grnas_per_targets": [
      "3"
    ],
    "treatment_ids": [
      "efo:0003317|efo:0003304",
      "ncit:c184729"
    ],
    "developmental_stage_labels": [
      "adolescent"
    ],
    "enzyme_expression_control_labels": [
      "constitutive expression"
    ],
    "library_scope_labels.text": [
      "focused"
    ],
    "sequencing_strategy_labels": [
      "barcode sequencing"
    ],
    "library_scope_labels": [
      "focused"
    ],
    "sequencing_library_kit_labels": [
      "nextera xt"
    ],
    "sequencing_platform_labels.text": [
      "Illumina HiSeq 4000"
    ],
    "library_integration_state_labels.text": [
      "random locus integration"
    ]
  }
}
```

As you can see it's quite long, but you must only use a fixed subset of fields both for search and for returning the results. Here is the subset: dataset_id, tissue_labels, cell_type_labels.text, cell_line_labels, sex_labels, developmental_stage_labels.text, disease_labels, library_perturbation_type_labels.

Some of these fields are lists and are named in plural (tissue_labels). You must flatteh the list, always assume it only contains a single entry (make an assert), and name output in singular: tissue_label.

# Data in Postgres
Inside the $PS_DB database, you must use a table named "perturb_seq". Hardcode this into the API.

An example of data row in Postgres:

dataset_id	perturbed_target_symbol	gene	log2foldchange	padj	basemean	max_ingested_at
adamson_2016_pilot	SPI1	CHST11	0.23117926187171009	0.50391124972854473	27.0810834104398	2025-09-30T13:37:59.073273

Use "perturbation_gene_name" API parameter to filter by perturbed_target_symbol, and "effect_gene_name" parameter to filter by gene.

# Output format
Output should conform to this schema: [(list of datasets)]
Each dataset is: {"dataset_id": ..., "tissue_label": ..., (other dataset fields), "perturbations": [(list of perturbations of a given dataset after filtering)}
Each perturbation is essentially a row from Postgres with some fields renamed: {"perturbation_gene_name": "SPI1", "base_mean": 27.0810834104398, "effect_gene_name": "CHST11", "log2foldchange": ..., "padj": ...}

# Data truncation
There can be a lot of data, so truncate it as follows:
* Always return all datasets, sorted alphabetically by their ID
* For a given dataset, the total number of *perturbations* must not exceed 20. Sort them in the increasing order of padj.
