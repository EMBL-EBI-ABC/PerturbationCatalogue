import pandas as pd
import pandera as pa
from pandera import Field, DataFrameModel
from pandera.typing import Series, Index, String, Int64, Float32
from typing import Optional
from pathlib import Path

# Get the absolute path to the current module file
module_path = Path(__file__).resolve()

# Navigate up to the project root and then to 'ontologies'
ont_dir = module_path.parent / "ontologies"

gene_ont = pd.read_parquet(ont_dir / "genes.parquet").drop_duplicates()
ctype_ont = pd.read_parquet(ont_dir / "cell_types.parquet").drop_duplicates()
cline_ont = pd.read_parquet(ont_dir / "cell_lines.parquet").drop_duplicates()
tis_ont = pd.read_parquet(ont_dir / "tissues.parquet").drop_duplicates()
dis_ont = pd.read_parquet(ont_dir / "diseases.parquet").drop_duplicates()


# adata.obs schema
class ObsSchema(DataFrameModel):
    dataset_id: Series[String] = Field(nullable=False)
    sample_id: Series[Int64] = Field(nullable=False)
    perturbation_name: Series[String] = Field(nullable=False)
    perturbed_target_coord: Series[String] = Field(nullable=True)
    perturbed_target_chromosome: Series[String] = Field(nullable=True)
    perturbed_target_chromosome_encoding: Int64 = Field(nullable=True, ge=0)
    perturbed_target_number: Series[Int64] = Field(nullable=False, ge=0)
    perturbed_target_ensg: Series[String] = Field(nullable=True)
    perturbed_target_symbol: Series[String] = Field(nullable=True)
    perturbed_target_biotype: Series[String] = Field(nullable=True)
    guide_sequence: Series[String] = Field(
        nullable=True,
        regex=r"^[ACGTN]+$",
        coerce=True,
    )
    perturbation_type_label: Series[String] = Field(nullable=False)
    perturbation_type_id: Series[String] = Field(nullable=True, str_contains=":")
    timepoint: Series[String] = Field(
        nullable=True, regex=r"^P\d+DT\d{1,2}H\d{1,2}M\d{1,2}S$"
    )
    treatment_label: Series[String] = Field(nullable=True)
    treatment_id: Series[String] = Field(nullable=True, str_contains=":")
    # model system details
    model_system_label: Series[String] = Field(nullable=False)
    model_system_id: Series[String] = Field(nullable=True, str_contains=":")
    species: Series[String] = Field(nullable=False, isin=["Homo sapiens"])
    tissue_label: Series[String] = Field(nullable=True)
    tissue_id: Series[String] = Field(nullable=True)
    cell_type_label: Series[String] = Field(nullable=False)
    cell_type_id: Series[String] = Field(nullable=False)
    cell_line_label: Series[String] = Field(nullable=True)
    cell_line_id: Series[String] = Field(nullable=True)
    sex_label: Series[String] = Field(nullable=True)
    sex_id: Series[String] = Field(nullable=True, str_contains=":")
    developmental_stage_label: Series[String] = Field(nullable=True)
    developmental_stage_id: Series[String] = Field(nullable=True, str_contains=":")
    disease_label: Series[String] = Field(nullable=True)
    disease_id: Series[String] = Field(nullable=True)
    # study details
    study_title: Series[String] = Field(nullable=False)
    study_uri: Series[String] = Field(nullable=False)
    study_year: Series[Int64] = Field(nullable=False, ge=1900, le=2100)
    first_author: Series[String] = Field(nullable=False)
    last_author: Series[String] = Field(nullable=False)
    # experiment details
    experiment_title: Series[String] = Field(nullable=False)
    experiment_summary: Series[String] = Field(nullable=True)
    number_of_perturbed_targets: Series[Int64] = Field(
        nullable=False,
        ge=1,
        description="Total number of perturbed targets in the experiment.",
    )
    number_of_perturbed_samples: Series[Int64] = Field(
        nullable=False,
        ge=1,
        description="Total number of perturbed samples/cells in the experiment.",
    )  # perturbation details
    library_generation_type_id: Series[String] = Field(nullable=True)
    library_generation_type_label: Series[String] = Field(nullable=True)
    library_generation_method_id: Series[String] = Field(nullable=True)
    library_generation_method_label: Series[String] = Field(nullable=True)
    enzyme_delivery_method_id: Series[String] = Field(nullable=True)
    enzyme_delivery_method_label: Series[String] = Field(nullable=True)
    library_delivery_method_id: Series[String] = Field(nullable=True)
    library_delivery_method_label: Series[String] = Field(nullable=True)
    enzyme_integration_state_id: Series[String] = Field(nullable=True)
    enzyme_integration_state_label: Series[String] = Field(nullable=True)
    library_integration_state_id: Series[String] = Field(nullable=True)
    library_integration_state_label: Series[String] = Field(nullable=True)
    enzyme_expression_control_id: Series[String] = Field(nullable=True)
    enzyme_expression_control_label: Series[String] = Field(nullable=True)
    # library details
    library_expression_control_id: Series[String] = Field(nullable=True)
    library_expression_control_label: Series[String] = Field(nullable=True)
    library_name: Series[String] = Field(nullable=True)
    library_uri: Series[String] = Field(nullable=True)
    library_format_id: Series[String] = Field(nullable=True)
    library_format_label: Series[String] = Field(nullable=True)
    library_scope_id: Series[String] = Field(nullable=True)
    library_scope_label: Series[String] = Field(nullable=True)
    library_perturbation_type_id: Series[String] = Field(nullable=True)
    library_perturbation_type_label: Series[String] = Field(nullable=True)
    library_manufacturer: Series[String] = Field(nullable=True)
    library_lentiviral_generation: Series[String] = Field(nullable=True)
    library_grnas_per_target: Series[String] = Field(nullable=True)
    library_total_grnas: Series[Int64] = Field(nullable=True, ge=0)
    library_total_variants: Int64 = Field(nullable=True, ge=0)
    # assay details
    readout_dimensionality_id: Series[String] = Field(nullable=True)
    readout_dimensionality_label: Series[String] = Field(nullable=True)
    readout_type_id: Series[String] = Field(nullable=True)
    readout_type_label: Series[String] = Field(nullable=True)
    readout_technology_id: Series[String] = Field(nullable=True)
    readout_technology_label: Series[String] = Field(nullable=True)
    method_name_id: Series[String] = Field(nullable=True)
    method_name_label: Series[String] = Field(nullable=True)
    method_uri: Series[String] = Field(nullable=True)
    sequencing_library_kit_id: Series[String] = Field(nullable=True)
    sequencing_library_kit_label: Series[String] = Field(nullable=True)
    sequencing_platform_id: Series[String] = Field(nullable=True)
    sequencing_platform_label: Series[String] = Field(nullable=True)
    sequencing_strategy_id: Series[String] = Field(nullable=True)
    sequencing_strategy_label: Series[String] = Field(nullable=True)
    software_counts_id: Series[String] = Field(nullable=True)
    software_counts_label: Series[String] = Field(nullable=True)
    software_analysis_id: Series[String] = Field(nullable=True)
    software_analysis_label: Series[String] = Field(nullable=True)
    reference_genome_id: Series[String] = Field(nullable=True)
    reference_genome_label: Series[String] = Field(nullable=True)
    # associated datasets
    associated_datasets: Series[String] = Field(
        nullable=True,
        coerce=True,
        description="List of associated datasets with each dataset having 'dataset_accession', 'dataset_uri', 'dataset_description', 'dataset_file_name' keys.",
    )

    class Config:
        strict = True
        coerce = False
        ordered = True


# adata.var schema
class VarSchema(DataFrameModel):
    index: Index[str] = Field(
        nullable=False,
        unique=True,
        # str_startswith=("ENSG", "control"),
        # isin=gene_ont.ensembl_gene_id.values,
        check_name=True,
    )
    ensembl_gene_id: Series[str] = Field(
        nullable=True,
        str_startswith=("ENSG", "control"),
    )
    gene_symbol: Series[str] = Field(
        nullable=True,
        coerce=True,
    )

    class Config:
        strict = "filter"
        coerce = True
        ordered = True
