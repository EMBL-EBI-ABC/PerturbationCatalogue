# Data Exploration and Curation Notebooks

This directory contains resources for exploring, standardizing, and curating genetic perturbation datasets (Perturb-seq, CRISPR and MAVE screens) into a unified, clean, and schema-compliant format.

## Environment Setup

Data exploration and curation notebooks use a common environment. To set it up on your machine, run the following commands from this directory:

```bash
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

Install the Jupyter kernel:
```bash
python -m ipykernel install --user --name=data_exploration_env
```

Verify that the kernel has been successfully added:
```bash
jupyter kernelspec list | grep data_exploration_env
```

Start the Jupyter notebook server:
```bash
jupyter notebook
```

Open any notebook and change the kernel to use `data_exploration_env`: **Kernel** → **Change kernel** → **data_exploration_env**.

---

## Data Curation Pipeline

To build a unified and searchable Perturbation Catalogue, datasets from various publications are standardised using the repository's custom `curation_tools` library. For ultimate data provenance tractability, the curation process is documented and executed inside the notebooks located in `Perturbseq|CRISPR|MaveDB/curation_notebooks/`.

### Unified Metadata Schema
The unified metadata schema is stored in [`data_exploration/curation_tools/perturbseq_anndata_schema.py`](curation_tools/perturbseq_anndata_schema.py). This schema (defined using `pandera`) enforces strict constraints, data types, and ontology validations for cell/observation-level metadata (`adata.obs` via `ObsSchema`) and feature-level metadata (`adata.var` via `VarSchema`).

A typical curation workflow consists of the following phases:

```
[ Ingest Raw Data ] ──> [ Curate OBS (Cells) ] ──> [ Curate VAR (Genes) ] ──> [ Validate & Export ]
```

### 1. Dataset Ingestion & Initialization
* **Retrieval**: Raw or pre-processed datasets (commonly `.h5ad` format) are retrieved programmatically using `download_file()` and placed in the `../non_curated/` directory.
* **Wrapper Initialization**: A `CuratedDataset` object is initialized, binding the raw dataset with target schemas:
  ```python
  cur_data = CuratedDataset(
      obs_schema=ObsSchema,
      var_schema=VarSchema,
      noncurated_path=noncurated_path
  )
  cur_data.load_data()
  ```

### 2. Observation (OBS) Slot Curation
The cells/observations metadata (`adata.obs`) undergo strict clean-up and enrichment:
* **Cell Barcode Standardization**: Assigning or formatting unique cell barcodes and linking them to technical replication units (e.g., lanes or GEM groups).
* **Quality Filtering**: Drop invalid cells, unassigned guides, or complex multi-guide RNA combinations where guide identities cannot be unambiguously resolved.
* **Guide & Spacer Mapping**: Linking guide IDs included within `.h5ad` files (e.g., `sgID_AB`) to their actual nucleotide sequences often provided separately in the supplemental metadata spreadsheets or Excel sheets from the source publications.
* **Perturbed Target Standardization**:
  * Original target gene symbols/Ensembl IDs are resolved and standardized using `cur_data.standardize_genes()`.
  * Non-targeting controls are standardized (e.g., mapped to `control_nontargeting`).
  * Gene positions and chromosome assignments are resolved and stored.
  * Chromosomes are mapped to standard integer values via `cur_data.chromosome_encoding()`.
* **Ontology Standardization**: Free-text metadata values are matched against controlled vocabularies using `cur_data.standardize_ontology()`:
  * **Tissues** mapped to [UBERON](https://www.ebi.ac.uk/ols4/ontologies/uberon) terms.
  * **Cell Lines** mapped to [Cell Line Ontology (CLO)](https://www.ebi.ac.uk/ols4/ontologies/clo)/[EFO](https://www.ebi.ac.uk/ols4/ontologies/efo) terms.
  * **Cell Types** mapped to [Cell Ontology (CL)](https://www.ebi.ac.uk/ols4/ontologies/cl) terms.
  * **Diseases** mapped to [MONDO](https://www.ebi.ac.uk/ols4/ontologies/mondo) / [NCIT](https://www.ebi.ac.uk/ols4/ontologies/ncit) / [EFO](https://www.ebi.ac.uk/ols4/ontologies/efo) disease terms.
* **Experiment Metadata Enrichment**: Over 50 experimental attributes (such as study title, DOI, sequencing platform, library kit, vector delivery, and enzyme delivery methods) are appended to every single-cell row for downstream querying.
* **Schema Matching**: Ordering and typing columns to match target database formats via `cur_data.match_schema_columns(slot='obs')`.

### 3. Feature (VAR) Slot Curation
The genes/features metadata (`adata.var`) are curated to ensure downstream expression matrices are easily joinable:
* **Gene ID Standardization**: Original indices or IDs are passed to `cur_data.standardize_genes(slot='var', ...)` to map features to Ensembl IDs and official HGNC symbols.
* **Unmapped Features & Controls**: Custom reporter genes (e.g., `PuroR`, green fluorescent protein/GFP reporter barcodes) are handled gracefully, ensuring their original indexes are retained without breaking validation.

### 4. Validation, Storage & Deployment
* **Enforced Schema Validation**: Validating the final table formats against strict schemas:
  ```python
  cur_data.validate_data(slot='obs')
  cur_data.validate_data(slot='var')
  ```
* **H5AD Serialization**: The final high-dimensional expression dataset with clean `obs` and `var` annotations is saved to the `curated/` folder using `cur_data.save_curated_data_h5ad()`.
* **Metadata Export**: Metadata tables are split off and exported into a lightweight `.parquet` format using `cur_data.save_curated_data_parquet(split_metadata=True, save_metadata_only=True)`.
* **Cloud Upload**:
  * Parquet metadata is uploaded to Google BigQuery for global relational queries:
    ```python
    upload_parquet_to_bq(
        parquet_path='../curated/parquet/[dataset]_metadata.parquet',
        bq_dataset_id='[project].perturb_seq',
        bq_table_name='metadata',
        key_columns=['dataset_id', 'sample_id']
    )
    ```
  * Curated `.h5ad` files are uploaded to Google Cloud Storage (GCS) lakes for scalable file-system access.

---

### Notebook examples
For real examples of curation notebooks, feel free to explore the notebooks located in e.g. [`data_exploration/Perturbseq/curation_notebooks`](data_exploration/Perturbseq/curation_notebooks).

## Metadata Data Dictionary

This section provides a clear, verbal reference of the metadata fields defined in the schema for non-technical users. These fields are defined in the schema code in [`data_exploration/curation_tools/perturbseq_anndata_schema.py`](curation_tools/perturbseq_anndata_schema.py).

### Observation Metadata (`adata.obs`)

The observation metadata table describes cell-level properties, experimental setup, and the specific genetic perturbations applied.

| Field Name | Description | Required / Optional | Example / Allowed Values |
| :--- | :--- | :--- | :--- |
| `dataset_id` | Unique identifier for the dataset, following the format `<firstauthor_year>` | **Required** | `smith_2024` |
| `sample_id` | Unique identifier for the sample | **Required** | `sample_1` |
| `cell_barcode` | Unique cell barcode | **Required** | `AAACCTGAGCTAGGCA-1` |
| `data_modality` | The assay's overall experimental approach | **Required** | `Perturb-seq`, `CRISPR screen`, `MAVE` |
| `significant` | Indicates whether the perturbation had a statistically significant effect | Optional | `True`, `False` |
| `significance_criteria` | Criteria used to determine significance | Optional | e.g., `FDR < 0.05` |
| `perturbation_name` | Name of the perturbation (e.g. gene symbol or genomic coordinates) | **Required** | `BRCA1` |
| `perturbed_target_coord` | Genomic coordinates of the perturbed target | Optional | `chr:start-end;strand` |
| `perturbed_target_chromosome` | Chromosome of the perturbed target | Optional | `chr17` |
| `perturbed_target_chromosome_encoding` | Numeric chromosome encoding for database partitioning | Optional | e.g. `17` |
| `perturbed_target_number` | Number of targets perturbed simultaneously in a cell | **Required** | `1` |
| `perturbed_target_ensg` | Ensembl gene ID(s) of the target (multiple joined by `+`) | Optional | `ENSG00000012048` or  `ENSG00000012048+ENSG00000139618` |
| `perturbed_target_symbol` | Official gene symbol(s) of the target (multiple joined by `+`) | Optional | `BRCA1` or `BRCA1+BRCA2` |
| `perturbed_target_id` | Canonical ID in `symbol\|Ensembl` format | Optional | `BRCA1\|ENSG00000012048` or `BRCA1\|ENSG00000012048+BRCA2\|ENSG00000139618` |
| `perturbed_target_biotype` | Gene biotype(s) of the perturbed target | Optional | `protein_coding` |
| `guide_sequence` | Guide RNA sequence (A, C, G, T, N characters only) in 5' to 3' direction | Optional | `ATGC...` |
| `perturbation_type_label` | Ontology term label describing the perturbation mechanism | **Required** | `CRISPRn`, `CRISPRi`, `CRISPRa`, `DMS` |
| `perturbation_type_id` | Ontology term ID for the perturbation mechanism | Optional | *yet to be defined* |
| `timepoint` | Timepoint of the sample in ISO 8601 duration format | Optional | e.g. `P1DT12H30M15S` |
| `treatment_label` | Name of the treatment/compound (ChEMBL label, or 'untreated control') | Optional | e.g. `untreated control` or `Interferon gamma` |
| `treatment_id` | ChEMBL compound ID for the treatment | Optional | e.g. `CHEMBL3286073` |
| `technical_replicate` | Identifier for technical replicate | Optional | `replicate_1` |
| `biological_replicate` | Identifier for biological replicate | Optional | `biorep_A` |
| `model_system_label` | High-level type of model system investigated | **Required** | `cell_line`, `primary_cell`, `organoid`, `yeast` |
| `model_system_id` | Ontology term ID for the model system | Optional | e.g. `CLO:0000031` |
| `species` | Organism species name | **Required** | `Homo sapiens` |
| `tissue_label` | Tissue ontology term label (UBERON) | Optional | e.g. `blood` |
| `tissue_id` | Tissue ontology term ID (UBERON) | Optional | e.g. `UBERON:0000178` |
| `cell_type_label` | Cell type ontology term label (Cell Ontology - CL) | Optional | e.g. `T cell` |
| `cell_type_id` | Cell type ontology term ID (Cell Ontology - CL) | Optional | e.g. `CL:0000084` |
| `cell_line_label` | Cell line ontology term label (Cell Line Ontology - CLO) | Optional | e.g. `HEK293H cell` |
| `cell_line_id` | Cell line ontology term ID (Cell Line Ontology - CLO) | Optional | e.g. `CLO:0037346` |
| `sex_label` | Biological sex of the sample | Optional | `female`, `male`, `mixed`, `unknown` |
| `sex_id` | Biological sex ontology ID | Optional | e.g. `PATO:0000383` |
| `developmental_stage_label`| Developmental/age stage label | Optional | `embryonic`, `fetal`, `neonatal`, `child`, `adolescent`, `adult`, `senior adult` |
| `developmental_stage_id` | Developmental stage ontology ID | Optional | e.g. `PATO:0001189` |
| `disease_label` | Disease ontology term label (MONDO) | Optional | e.g. `breast cancer` |
| `disease_id` | Disease ontology term ID (MONDO) | Optional | e.g. `MONDO:0007254` |
| `study_title` | Title of the publication/study | **Required** | e.g. `Genome-wide CRISPR screens...` |
| `study_uri` | Web link or DOI of the study | **Required** | e.g. `https://doi.org/...` |
| `study_year` | Publication year | **Required** | e.g. `2024` |
| `first_author` | Full name of the first author | Optional | `Jane Smith` |
| `last_author` | Full name of the last author | Optional | `John Doe` |
| `experiment_title` | Descriptive title of the experiment | **Required** | `K562 CRISPRi Screen` |
| `experiment_summary` | Verbose summary explaining the experiment | Optional | `A high-throughput screen to...` |
| `number_of_perturbed_targets` | Total number of unique perturbed targets in the experiment | **Required** | `50` |
| `number_of_perturbed_samples` | Total count of cells or samples perturbed in this experiment | Optional | `10000` |
| `library_generation_type_id` | EFO ID for library generation type (under `EFO:0022867`) | Optional | e.g. `EFO:0022868` |
| `library_generation_type_label` | EFO label for library generation type | Optional | `Endogenous genetic perturbation method`, `Exogenous genetic perturbation method` |
| `library_generation_method_id` | EFO ID for library generation method (subtype of `EFO:0022868` or `EFO:0022869`) | Optional | e.g. `EFO:0022876` |
| `library_generation_method_label` | EFO label for library generation method | Optional | e.g. `SpCas9` |
| `enzyme_delivery_method_id` | Ontology ID for enzyme delivery (under `EFO:0920064`) | Optional | e.g. `EFO:0920067` |
| `enzyme_delivery_method_label` | Method used to deliver the editing enzyme | Optional | `electroporation`, `lipofection`, `microinjection`, `nanoparticle-based transfection`,`nucleofection`,`optical transfection`,`sonoporation`,`lentivirus transduction`,`retrovirus transduction`,`adeno-associated virus transduction`,`adenovirus transduction`,`herpes virus transduction` |
| `library_delivery_method_id` | Ontology ID for library delivery (under `EFO:0920064`) | Optional | e.g. `EFO:0920067` |
| `library_delivery_method_label` | Method used to deliver the perturbation library | Optional | `electroporation`, `lipofection`, `microinjection`, `nanoparticle-based transfection`,`nucleofection`,`optical transfection`,`sonoporation`,`lentivirus transduction`,`retrovirus transduction`,`adeno-associated virus transduction`,`adenovirus transduction`,`herpes virus transduction` |
| `enzyme_integration_state_id` | Ontology ID for integration state of enzyme (under `EFO:0920081`) | Optional | e.g. `EFO:0920082` |
| `enzyme_integration_state_label` | Integration state of editing enzyme transgene | Optional | `random locus integration`, `targeted locus integration`, `native locus replacement`, `non-integrative transgene expression` |
| `library_integration_state_id` | Ontology ID for integration state of library (under `EFO:0920081`) | Optional | e.g. `EFO:0920082` |
| `library_integration_state_label` | Integration state of library construct | Optional | `random locus integration`, `targeted locus integration`, `native locus replacement`, `non-integrative transgene expression` |
| `enzyme_expression_control_id` | Ontology ID for expression control of enzyme | Optional | *yet to be defined* |
| `enzyme_expression_control_label` | Transcriptional control of the enzyme | Optional | `constitutive transgene expression`, `inducible transgene expression`, `native promoter-driven transgene expression`, `degradation domain-based transgene control` |
| `library_expression_control_id` | Ontology ID for expression control of library | Optional | *yet to be defined* |
| `library_expression_control_label` | Transcriptional control of the guide/library | Optional | `constitutive transgene expression`, `inducible transgene expression`, `native promoter-driven transgene expression`, `degradation domain-based transgene control` |
| `library_name` | Full name of the perturbation library | Optional | `Bassik Human CRISPR Knockout Library` |
| `library_uri` | Web link or accession for the library | Optional | e.g. `https://www.addgene.org/...` |
| `library_format_id` | Ontology ID for library format | Optional | *yet to be defined* |
| `library_format_label` | Presentation format of the library | Optional | `pooled`, `arrayed`, `arrayed\|pooled`, `in vivo` |
| `library_scope_id` | Ontology ID for library scope | Optional | *yet to be defined* |
| `library_scope_label` | Biological scope of the library | Optional | `focused`, `genome-wide` |
| `library_perturbation_type_id` | Ontology ID for library perturbation type | Optional | *yet to be defined* |
| `library_perturbation_type_label` | Functional type of library perturbation | Optional | `knockout`, `inhibition`, `activation`, `base editing`, `prime editing`, `mutagenesis` |
| `library_manufacturer` | Vendor, manufacturer, or origin lab of the library | Optional | `Bassik` |
| `library_lentiviral_generation` | Lentiviral packaging generation number | Optional | `3` |
| `library_grnas_per_target` | Average number of gRNAs designed per target gene | Optional | `4`, `5-7` |
| `library_total_grnas` | Total number of gRNAs contained in the library | Optional | e.g. `20,000` |
| `library_total_variants` | Total number of unique variants (only applicable for MAVE) | Optional | e.g. `5,000` |
| `readout_dimensionality_id` | Ontology ID for assay dimensionality | Optional | *yet to be defined* |
| `readout_dimensionality_label` | Dimensionality classification of the readout assay | Optional | `single-dimensional assay`, `high-dimensional assay` |
| `readout_type_id` | Ontology ID for readout type | Optional | *yet to be defined* |
| `readout_type_label` | Analyzed molecule category or system phenotypic readout | Optional | `transcriptomic`, `proteomic`, `phenotypic` |
| `readout_technology_id` | Ontology ID for readout technology | Optional | *yet to be defined* |
| `readout_technology_label`| General class of readout technology | Optional | `single-cell rna-seq`, `population growth assay`, `flow cytometry` |
| `method_name_id` | Ontology ID for the methodology name | Optional | *yet to be defined* |
| `method_name_label` | Name of the experimental profiling assay | Optional | `Perturb-seq`, `Perturb-CITE-seq`, `scRNA-seq`, `proliferation CRISPR screen`, `DMS-TileSeq`, `DMS-BarSeq`, `Joined and refined DMS-BarSeq and DMS-TileSeq`, `Combined DMS-BarSeq and DMS-TileSeq` |
| `method_uri` | Web link detailing the assay methodology | Optional | e.g. `https://doi.org/10.1038/...` |
| `sequencing_library_kit_id` | Ontology ID for the library prep kit | Optional | *yet to be defined* |
| `sequencing_library_kit_label` | Commercial kit name used for sequencing library preparation | Optional | `10x Genomics Chromium GEM-X Single Cell 5-prime kit v3`, `10x Genomics Chromium Next GEM Single Cell 5-prime HT Kit v2`, `10x Genomics Single Cell 3-prime`, `10x Genomics Single Cell 3-prime v2`, `10x Genomics Single Cell 3-prime v3`, `Nextera XT DNA Library Preparation Kit`, `GEM-X Flex Gene Expression Human n-plex kit`, `Parse Biosciences Evercode Whole Transcriptome Mega v1 kit` |
| `sequencing_platform_id` | Ontology ID for the sequencing machine platform | Optional | *yet to be defined* |
| `sequencing_platform_label` | Specific instrument model used to generate sequence data | Optional | `Illumina NovaSeq X`, `Illumina NovaSeq X Plus`, `Illumina HiSeq 4000`, `Illumina HiSeq 2500`, `Illumina HiSeq 2000`, `Illumina NovaSeq 6000`, `Illumina NextSeq 500`, `Ultima Genomics UG100` |
| `sequencing_strategy_id` | Ontology ID for sequencing strategy | Optional | *yet to be defined* |
| `sequencing_strategy_label` | High-level sequencing protocol approach | Optional | `barcode sequencing`, `direct sequencing`, `barcode sequencing\|direct sequencing` |
| `software_counts_id` | Ontology ID for counting software | Optional | *yet to be defined* |
| `software_counts_label` | Software pipeline used to align and generate raw counts | Optional | `custom`, `MaGeCK`, `CellRanger`, `Drop-seq Tools` |
| `software_analysis_id` | Ontology ID for analysis software | Optional | *yet to be defined* |
| `software_analysis_label` | Software tools used for downstream analysis/filtering | Optional | `custom`, `MAGeCK`, `Achilles`, `TRADE`, `Seurat`, `MAST`, `scanpy` |
| `score_interpretation` | Verbal interpretation of how the perturbation score is read | Optional | e.g., `positive score indicates increased abundance` |
| `reference_genome_id` | Ontology ID for reference genome | Optional | *yet to be defined* |
| `reference_genome_label` | Target reference assembly version used | Optional | `GRCh38`, `GRCh37` |
| `associated_datasets` | Associated accession numbers or metadata links | Optional | e.g. `GSE123456` |
| `license_label` | Software/Data distribution license name | **Required** | e.g. `Creative Commons Attribution 4.0` |
| `license_id` | SWO ontology term ID for the license | Optional | e.g. `SWO:0000002` |

### Feature Metadata (`adata.var`)

The feature metadata table describes the properties of the analyzed genetic readouts (typically genes or transcripts).

| Field Name | Description | Required / Optional | Example / Format |
| :--- | :--- | :--- | :--- |
| `index` | Unique identifier for each gene feature (typically Ensembl ID) | **Required** | `ENSG0000012048` |
| `ensembl_gene_id` | Standardized Ensembl gene ID | Optional | Starts with `ENSG` or `control` |
| `gene_symbol` | Official gene symbol | Optional | `BRCA1` |
