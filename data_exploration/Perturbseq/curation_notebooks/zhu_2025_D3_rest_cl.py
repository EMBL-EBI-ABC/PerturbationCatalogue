# %% Import
import pandas as pd
import json
from datetime import datetime

from curation_tools.curation_tools import (
    CuratedDataset,
    ObsSchema,
    VarSchema,
    Experiment,
    download_file,
    upload_parquet_to_bq
)

import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=[
        logging.FileHandler("curation.log"),
        logging.StreamHandler(),
    ],
    force=True,
)
pd.set_option('display.max_columns', None)

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

log("=== zhu_2025_D3_rest_cl curation started ===")

# %% Download data
# Download the data from AWS as such:
# !aws s3 cp --no-sign-request s3://genome-scale-tcell-perturb-seq/marson2025_data/{name_of_the_file}.h5ad ..aa/non_curated/h5ad/{name_of_the_file}.h5ad

# %% Initialise the dataset object
noncurated_path = '/hps/nobackup/mfreeberg/marson_downloads/zhu_2025_D3_rest_cl.h5ad'
cur_data = CuratedDataset(
    obs_schema=ObsSchema,
    var_schema=VarSchema,
    exp_metadata_schema=Experiment,
    noncurated_path=noncurated_path
)

log(f"Loading data from {noncurated_path}...")
cur_data.load_data()
log(f"Data loaded: {cur_data.adata.n_obs} cells, {cur_data.adata.n_vars} vars")

# %% OBS slot curation
# Since for multi-guide perturbations the identities of said guides are unknown,
# we are filtering out these cells.
log(f"Filtering multi_sgRNA cells (before: {cur_data.adata.n_obs} cells)...")
cur_data.adata = cur_data.adata[cur_data.adata.obs['guide_id'] != 'multi_sgRNA']
log(f"Filtering done (after: {cur_data.adata.n_obs} cells)")

# %% Add index as perturbation_name
log("Adding perturbation_name...")
cur_data.adata.obs['cell_barcode'] = cur_data.adata.obs.index.str.split('_').str[0] + '_' + 'D3-REST'
cur_data.adata.obs['perturbation_name'] = cur_data.adata.obs['cell_barcode'] + '_' + cur_data.adata.obs['lane_id'].astype(str)

# %% Add guide RNA information
log("Downloading guide RNA spreadsheet...")
download_file(
    url="https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/refs/heads/master/metadata/suppl_tables/sgrna_library_metadata.suppl_table.csv",
    dest_path="../supplementary/zhu_2025_guide_info.csv"
)
# Read in the guide RNA info csv
guide_info_df = pd.read_csv("../supplementary/zhu_2025_guide_info.csv")

guide_info_df = guide_info_df[['sgRNA', 'seq']]

guide_info_df['sgRNA'] = guide_info_df['sgRNA'].str.replace('1-Jun', 'JUN-1').str.replace('2-Jun', 'JUN-2')

guide_info_df = guide_info_df.rename(columns={'seq': 'guide_sequence', 'sgRNA': 'guide_id'})

# %% Check guide overlap
log(f"Checking guide overlap ({len(guide_info_df)} guides in reference)...")
print(cur_data.adata.obs['guide_id'].isin(guide_info_df['guide_id'].to_list()).value_counts())

# %% Merge guide info into obs
log("Merging guide info into obs...")
cur_data.adata.obs = cur_data.adata.obs.merge(guide_info_df, on='guide_id', how='left')
log(f"Merge done. Missing guide sequences: {cur_data.adata.obs['guide_sequence'].isna().sum()}")

# %% Fix perturbed gene id/name types and set control labels
log("Setting control labels...")
cur_data.adata.obs['perturbed_gene_id'] = cur_data.adata.obs['perturbed_gene_id'].astype("string")
cur_data.adata.obs['perturbed_gene_name'] = cur_data.adata.obs['perturbed_gene_name'].astype("string")

cur_data.adata.obs.loc[cur_data.adata.obs['perturbed_gene_id'].isin(['NTC']), 'perturbed_gene_id'] = 'control_nontargeting'
cur_data.adata.obs.loc[cur_data.adata.obs['perturbed_gene_name'].isin(['NTC']), 'perturbed_gene_name'] = 'control_nontargeting'

cur_data.adata.obs.loc[cur_data.adata.obs['perturbed_gene_id'].isna(), 'perturbed_gene_id'] = 'control_casonly'
cur_data.adata.obs.loc[cur_data.adata.obs['perturbed_gene_name'].isna(), 'perturbed_gene_name'] = 'control_casonly'

# %% Standardise perturbation targets
log("Standardising perturbation target genes...")
cur_data.standardize_genes(
    slot='obs',
    input_column='perturbed_gene_id',
    input_column_type='ensembl_gene_id',
    multiple_entries=False,
)

# %% Manually replace some genes
genes_to_replace = cur_data.adata.obs[cur_data.adata.obs['perturbed_target_symbol'].isna()]['perturbed_gene_name'].drop_duplicates().to_list()
log(f"Manually replacing {len(genes_to_replace)} unmapped genes...")
cols_to_replace = ['perturbed_target_ensg', 'perturbed_target_symbol', 'perturbed_target_biotype', 'perturbed_target_coord', 'perturbed_target_chromosome']
gene_ont_col_mapping = {
    'ensembl_gene_id': 'perturbed_target_ensg',
    'gene_symbol': 'perturbed_target_symbol',
    'biotype': 'perturbed_target_biotype',
    'gene_coord': 'perturbed_target_coord',
    'chromosome_name': 'perturbed_target_chromosome'
}
for replacement_gene in genes_to_replace:
    if replacement_gene in cur_data.gene_ont['synonym'].values:
        replacement = (
            cur_data.gene_ont.rename(columns=gene_ont_col_mapping)
            .loc[cur_data.gene_ont['synonym'].isin([replacement_gene]), cols_to_replace]
        )[cols_to_replace].values[0]
        cur_data.adata.obs.loc[cur_data.adata.obs['perturbed_gene_name'].isin([replacement_gene]), cols_to_replace] = replacement
        log(f"  Replaced {replacement_gene} with {replacement}")
    else:
        log(f"  No replacement found for {replacement_gene}")

# %% Replace non-mapped perturbed_target_symbol with original gene symbols
cur_data.adata.obs.loc[cur_data.adata.obs['perturbed_target_symbol'].isna(), 'perturbed_target_symbol'] = cur_data.adata.obs.loc[cur_data.adata.obs['perturbed_target_symbol'].isna(), 'perturbed_gene_name'].values

# %% Replace non-ENSEMBL perturbed_target_ensg entries with NA
cur_data.adata.obs.loc[(~cur_data.adata.obs['perturbed_target_ensg'].str.startswith('ENSG')) & (cur_data.adata.obs['perturbed_target_ensg'] != 'control_nontargeting'), 'perturbed_target_ensg'] = pd.NA

# %% Add perturbed_target_number column
cur_data.adata.obs['perturbed_target_number'] = 1

# %% Encode chromosomes as integers
log("Encoding chromosomes...")
cur_data.chromosome_encoding()

# %% Curate replicates
cur_data.adata.obs = cur_data.adata.obs.rename(columns={'lane_id': 'technical_replicate'})

# %% Add metadata
log("Adding obs metadata columns...")
cur_data.create_columns(
    overwrite=True,
    slot="obs",
    col_dict={
        "dataset_id": cur_data.dataset_id,
        "sample_id": range(1, cur_data.adata.obs.shape[0] + 1),
        # perturbation type
        "perturbation_type_label": "CRISPRi",
        "perturbation_type_id": None,
        "data_modality": "Perturb-seq",
        "significant": None,
        "significance_criteria": None,
        "score_interpretation": None,

        # treatment
        "treatment_label": None, #"anti-CD3 antibody|anti-CD28 antibody|anti-CD2 antibody",
        "treatment_id": None, #"EFO:0003317|EFO:0003304|NCIT:C184729",
        # replicates
        "biological_replicate": "D3_CE0008678",
        # model system
        "model_system_label": "primary_cell",
        "model_system_id": None,
        "tissue": "lymphoid tissue",
        "cell_line_label": None,
        "cell_line_id": None,
        "cell_type_label": "CD4-positive, alpha-beta T cell",
        "disease_label": "healthy",
        "disease_id": None,

        "timepoint": "P12DT8H0M0S",
        "species": "Homo sapiens",
        "sex_label": "male",
        "sex_id": None,
        "developmental_stage_label": "adult",
        "developmental_stage_id": None,

        "study_title": "Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits",
        "study_uri": "https://doi.org/10.64898/2025.12.23.696273",
        "study_year": 2025,
        "first_author": "Ronghui Zhu",
        "last_author": "Alexander Marson",

        "experiment_title": "Perturb-seq of primary human CD4-positive T cells for patient D3_CE0008678 under resting conditions",
        "experiment_summary": """
            Isolated human CD4-positive cells from four healthy donors were stimulated with ImmunoCult CD3/CD28/CD2 activator and sequentially transduced with dCas9-KRAB-Zim3 lentivirus (next morning after stimulation) and a Perturb-seq guide library (next afternoon after stimulation, MOI 0.2).
            The library consisted all genes expressed in human CD4+ T cells, all transcription factors annotated in the Lambert et al. (2018) and non-targeting controls totalling 12,748 genes.
            The gRNA sequences were selected from hCRISPRiv2 and Dolcetto libraries.
            On day 12, cells were split into three conditions:
            (1) Rest - cells were left for 8 hr without stimulation;
            (2) Stim8hr - cells were stimulated with ImmunoCult CD3/CD28/CD2 activator for 8 hr;
            (3) Stim48hr - cells were stimulated with ImmunoCult CD3/CD28/CD2 activator for 48 hr.
            Cells were harvested, fixed, stored using GEM-X Flex Sample Preparation v2 Kit and sequenced using GEM-X Flex Gene Expression Human n-plex kit converted into Ultima compatible libraries for sequencing on the Ultima Genomics UG100.
            """,

        "number_of_perturbed_targets": len(set(cur_data.adata.obs['perturbed_target_symbol'])),
        "number_of_perturbed_samples": cur_data.adata.obs.shape[0],

        "library_generation_type_id": "EFO:0022868",
        "library_generation_type_label": "endogenous",

        "library_generation_method_id": None,
        "library_generation_method_label": "dCas9-KRAB-Zim3",

        "enzyme_delivery_method_id": None,
        "enzyme_delivery_method_label": "lentivirus transduction",

        "library_delivery_method_id": None,
        "library_delivery_method_label": "lentivirus transduction",

        "enzyme_integration_state_id": None,
        "enzyme_integration_state_label": "random locus integration",

        "library_integration_state_id": None,
        "library_integration_state_label": "random locus integration",

        "enzyme_expression_control_id": None,
        "enzyme_expression_control_label": "constitutive transgene expression",

        "library_expression_control_id": None,
        "library_expression_control_label": "constitutive transgene expression",

        "library_name": "custom",
        "library_uri": None,

        "library_format_id": None,
        "library_format_label": "pooled",

        "library_scope_id": None,
        "library_scope_label": "focused",

        "library_perturbation_type_id": None,
        "library_perturbation_type_label": "inhibition",

        "library_manufacturer": "Marson lab",
        "library_lentiviral_generation": "2",
        "library_grnas_per_target": "2",
        "library_total_grnas": str(cur_data.adata.obs['guide_sequence'].str.split('|').explode().nunique()),
        "library_total_variants": None,

        "readout_dimensionality_id": None,
        "readout_dimensionality_label": "high-dimensional assay",

        "readout_type_id": None,
        "readout_type_label": "transcriptomic",

        "readout_technology_id": None,
        "readout_technology_label": "single-cell rna-seq",

        "method_name_id": None,
        "method_name_label": "Perturb-seq",

        "method_uri": None,

        "sequencing_library_kit_id": None,
        "sequencing_library_kit_label": "GEM-X Flex Gene Expression Human n-plex kit",

        "sequencing_platform_id": None,
        "sequencing_platform_label": "Ultima Genomics UG100",

        "sequencing_strategy_id": None,
        "sequencing_strategy_label": "barcode sequencing",

        "software_counts_id": None,
        "software_counts_label": "CellRanger",

        "software_analysis_id": None,
        "software_analysis_label": "scanpy",

        "reference_genome_id": None,
        "reference_genome_label": "GRCh38",

        "license_label": "MIT License",
        "license_id": "SWO:9000074",

        "associated_datasets": json.dumps([
            {
                "dataset_accession": "Primary Human CD4+ T Cell Perturb-seq",
                "dataset_uri": "s3://genome-scale-tcell-perturb-seq/marson2025_data/D3_Rest.assigned_guide.h5ad",
                "dataset_description": "Cell expression profiles for cells from donor D3_CE0008678, stimulated for 8hr with ImmunoCult CD3/CD28/CD2 activator.",
                "dataset_file_name": "D3_Rest.assigned_guide.h5ad",
            }
        ])
    }
)
log("Obs metadata added.")

# %% Curate tissue information
log("Standardising ontology: tissue...")
cur_data.standardize_ontology(
    input_column='tissue',
    column_type='term_name',
    ontology_type='tissue',
    overwrite=True
)

# %% Curate cell type information
log("Standardising ontology: cell_type...")
cur_data.standardize_ontology(
    input_column='cell_type_label',
    column_type='term_name',
    ontology_type='cell_type',
    overwrite=True
)

# %% Curate disease information
log("Standardising ontology: disease...")
cur_data.standardize_ontology(
    input_column='disease_label',
    column_type='term_name',
    ontology_type='disease',
    overwrite=True
)

# %% Match schema column order
log("Matching schema column order...")
cur_data.match_schema_columns(slot='obs')

# %% Validate obs metadata
log("Validating obs...")
cur_data.validate_data(slot='obs', verbose=True)
log("Obs validation done.")

# %% Standardise genes in var
log("Standardising var genes...")
cur_data.standardize_genes(
    slot="var",
    input_column="gene_ids",
    input_column_type="ensembl_gene_id",
    remove_version=False,
    multiple_entries=False
)

# %% Replace non-ENSG gene symbols with original gene name
cur_data.adata.var.loc[
    (cur_data.adata.var['gene_symbol'].isna()) &
    (~cur_data.adata.var['gene_name'].str.startswith('ENSG')),
    'gene_symbol'
] = cur_data.adata.var.loc[
    (cur_data.adata.var['gene_symbol'].isna()) &
    (~cur_data.adata.var['gene_name'].str.startswith('ENSG')),
    'gene_name'
]

# %% Fix PuroR ensembl_gene_id
cur_data.adata.var.loc[cur_data.adata.var['ensembl_gene_id'] == 'CUSTOM001_PuroR', 'ensembl_gene_id'] = None

# %% Validate var metadata
log("Validating var...")
cur_data.validate_data(slot='var')
log("Var validation done.")

# %% Save the dataset
log("Saving h5ad...")
cur_data.save_curated_data_h5ad()
log("h5ad saved.")

log("Saving parquet...")
cur_data.save_curated_data_parquet(split_metadata=True, save_metadata_only=True)
log("Parquet saved.")

# %% Upload to BigQuery (commented out)
upload_parquet_to_bq(
    parquet_path='/hps/nobackup/mfreeberg/marson_downloads/zhu_2025_D3_rest_cl_curated_metadata.parquet',
    bq_dataset_id='prj-ext-dev-pertcat-437314.perturb_seq',
    bq_table_name='metadata',
    key_columns=['dataset_id', 'sample_id'],
    verbose=True
)

# %% Upload to GC Storage (commented out)
# !gcloud storage cp /hps/nobackup/mfreeberg/marson_downloads/zhu_2025_D3_rest_cl_curated.h5ad gs://perturbation-catalogue-lake/perturbseq/curated/

log("=== Curation complete ===")
