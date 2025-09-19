import pandas as pd
from pandera import Field, DataFrameModel
from pandera.typing import Series, Index, String, Int64, Float32
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
    dataset_id: Series[String] = Field(
        nullable=False,
        description="Unique identifier for the dataset, follows the format <firstauthor_year>",
    )
    sample_id: Series[String] = Field(
        nullable=False, description="Unique identifier for the sample."
    )
    data_modality: Series[String] = Field(
        nullable=False,
        description="Data modality of the dataset.",
        isin=["Perturb-seq", "CRISPR screen", "MAVE"],
    )
    significant: Series[String] = Field(
        nullable=True, description="Indicates whether the perturbation had a significant effect.",
        coerce=True,
        isin=["True", "False"]
    )
    significance_criteria: Series[String] = Field(
        nullable=True, description="Criteria used to determine significance, e.g., FDR < 0.05."
    )
    perturbation_name: Series[String] = Field(
        nullable=False,
        description="Name of the perturbation, often a name of the targeted gene or genomic coordinate.",
    )
    perturbed_target_coord: Series[String] = Field(
        nullable=True,
        description="Genomic coordinates of the perturbed target. Format: chr:start-end;strand",
    )
    perturbed_target_chromosome: Series[String] = Field(
        nullable=True, description="Chromosome of the perturbed target."
    )
    perturbed_target_chromosome_encoding: Int64 = Field(
        nullable=True,
        ge=0,
        description="Numeric encoding of the chromosome of the perturbed target. Required for data partitioning in BigQuery.",
    )
    perturbed_target_number: Series[Int64] = Field(
        nullable=False, ge=0, description="Number of perturbed targets in the samples."
    )
    perturbed_target_ensg: Series[String] = Field(
        nullable=True, description="Ensembl gene ID(s) of the perturbed target."
    )
    perturbed_target_symbol: Series[String] = Field(
        nullable=True, description="Gene symbol(s) of the perturbed target."
    )
    perturbed_target_biotype: Series[String] = Field(
        nullable=True, description="Biotype(s) of the perturbed target."
    )
    guide_sequence: Series[String] = Field(
        nullable=True,
        regex=r"^[ACGTN]+$",
        coerce=True,
        description="Guide RNA sequence in 5' to 3' direction consisting of A, C, G, T, N characters only.",
    )
    perturbation_type_label: Series[String] = Field(
        nullable=False,
        description="Perturbation type ontology term label of the investigated sample.",
        isin=["CRISPRn", "CRISPRi", "CRISPRa", "DMS"],
    )
    perturbation_type_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Perturbation type ontology term ID of the investigated sample.",
    )
    perturbation_type_label: Series[String] = Field(
        nullable=False,
        description="Perturbation type ontology term label of the investigated sample.",
        isin=["CRISPRn", "CRISPRi", "CRISPRa", "DMS"],
    )
    perturbation_type_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Perturbation type ontology term ID of the investigated sample.",
    )
    timepoint: Series[String] = Field(
        nullable=True,
        regex=r"^P\d+DT\d{1,2}H\d{1,2}M\d{1,2}S$",
        description="Timepoint of the investigated sample in ISO 8601 format. Example: P1DT12H30M15S",
    )
    treatment_label: Series[String] = Field(
        nullable=True,
        description="Treatment/compound ontology term label used to stimulate the investigated sample. ChEMBL compound label.",
    )
    treatment_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Treatment/compound ontology term ID used to stimulate the investigated sample. ChEMBL compound ID.",
    )
    # model system details
    model_system_label: Series[String] = Field(
        nullable=False,
        description="Model system ontology term label of the investigated sample.",
        isin=["cell_line", "primary_cell", "organoid", 'yeast'],
    )
    model_system_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Model system ontology term ID of the investigated sample.",
    )
    species: Series[String] = Field(
        nullable=False,
        description="Species name of the investigated sample.",
        isin=["Homo sapiens"],
    )
    tissue_label: Series[String] = Field(
        nullable=True,
        description="Tissue ontology term label of the investigated sample. Must be part of the UBERON ontology.",
    )
    tissue_id: Series[String] = Field(
        nullable=True,
        description="Tissue ontology term ID of the investigated sample. Must be part of the UBERON ontology.",
    )
    cell_type_label: Series[String] = Field(
        nullable=True,
        description="Cell type ontology term label of the investigated sample. Must be part of the Cell Ontology (CL).",
    )
    cell_type_id: Series[String] = Field(
        nullable=True,
        description="Cell type ontology term ID of the investigated sample. Must be part of the Cell Ontology (CL).",
    )
    cell_line_label: Series[String] = Field(
        nullable=True,
        description="Cell line ontology term label of the investigated sample. Must be part of the Cell Line Ontology (CLO).",
    )
    cell_line_id: Series[String] = Field(
        nullable=True,
        description="Cell line ontology term ID of the investigated sample. Must be part of the Cell Line Ontology (CLO).",
    )
    sex_label: Series[String] = Field(
        nullable=True,
        description="Sex ontology term label of the investigated sample.",
        isin=["female", "male", "mixed", "unknown"]
    )
    sex_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Sex ontology term ID of the investigated sample.",
    )
    developmental_stage_label: Series[String] = Field(
        nullable=True,
        description="Developmental stage ontology term label of the investigated sample.",
        isin=["embryonic", "fetal", "neonatal", "child", "adolescent", "adult", "senior adult"],
    )
    developmental_stage_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Developmental stage ontology term ID of the investigated sample.",
    )
    disease_label: Series[String] = Field(
        nullable=True,
        description="Disease ontology term label of the investigated sample. Must be part of the MONDO ontology.",
    )
    disease_id: Series[String] = Field(
        nullable=True,
        description="Disease ontology term ID of the investigated sample. Must be part of the MONDO ontology.",
    )
    # study details
    study_title: Series[String] = Field(
        nullable=False, description="Title of the study/publication."
    )
    study_uri: Series[String] = Field(
        nullable=False, description="URI/DOI of the study/publication."
    )
    study_year: Series[Int64] = Field(
        nullable=False,
        ge=1900,
        le=2100,
        description="Publication year of the study/publication.",
    )
    first_author: Series[String] = Field(
        nullable=True,
        description="Full name of the first author of the study/publication.",
    )
    last_author: Series[String] = Field(
        nullable=True,
        description="Full name of the last author of the study/publication.",
    )
    # experiment details
    experiment_title: Series[String] = Field(
        nullable=False, description="Title of the experiment."
    )
    experiment_summary: Series[String] = Field(
        nullable=True, description="Summary of the experiment."
    )
    number_of_perturbed_targets: Series[String] = Field(
        nullable=False,
        coerce=True,
        description="Total number of perturbed targets in the experiment.",
    )
    number_of_perturbed_samples: Series[String] = Field(
        nullable=True,
        coerce=True,
        description="Total number of perturbed samples/cells in the experiment.",
    )  # perturbation details
    library_generation_type_id: Series[String] = Field(
        nullable=True,
        description="Library generation type ontology term ID, defined in EFO under parent term EFO:0022867 (genetic perturbation)",
    )
    library_generation_type_label: Series[String] = Field(
        nullable=True,
        description="Library generation type ontology term label, defined in EFO under parent term EFO:0022867 (genetic perturbation)",
    )
    library_generation_method_id: Series[String] = Field(
        nullable=True,
        description="Library generation method ontology term ID, defined in EFO under parent term EFO:0022868/EFO:0022869 (Endogenous/Exogenous genetic perturbation method)",
    )
    library_generation_method_label: Series[String] = Field(
        nullable=True,
        description="Library generation method ontology term label, defined in EFO under parent term EFO:0022868/EFO:0022869 (Endogenous/Exogenous genetic perturbation method)",
    )
    enzyme_delivery_method_id: Series[String] = Field(
        nullable=True,
        description="Enzyme delivery method ontology term ID.",
    )
    enzyme_delivery_method_label: Series[String] = Field(
        nullable=True,
        description="Enzyme delivery method ontology term label.",
        isin=["lipofection", "nucleofection", "retrovirus transduction", "lentivirus transduction", "transformation"]
    )
    library_delivery_method_id: Series[String] = Field(
        nullable=True, description="Library delivery method ontology term ID."
    )
    library_delivery_method_label: Series[String] = Field(
        nullable=True,
        description="Library delivery method ontology term label.",
        isin=["lipofection", "nucleofection", "retrovirus transduction", "lentivirus transduction", "transformation"]
    )
    enzyme_integration_state_id: Series[String] = Field(
        nullable=True, description="Enzyme integration state ontology term ID."
    )
    enzyme_integration_state_label: Series[String] = Field(
        nullable=True,
        description="Enzyme integration state ontology term label.",
        isin=["random locus integration", "targeted locus integration", "native locus replacement", "non-integrative transgene expression"]
    )
    library_integration_state_id: Series[String] = Field(
        nullable=True, description="Library integration state ontology term ID."
    )
    library_integration_state_label: Series[String] = Field(
        nullable=True,
        description="Library integration state ontology term label.",
        isin=["random locus integration", "targeted locus integration", "native locus replacement", "non-integrative transgene expression"]
    )
    enzyme_expression_control_id: Series[String] = Field(
        nullable=True, description="Enzyme expression control ontology term ID."
    )
    enzyme_expression_control_label: Series[String] = Field(
        nullable=True,
        description="Enzyme expression control ontology term label.",
        isin=["constitutive transgene expression", "inducible transgene expression", "native promoter-driven transgene expression", "degradation domain-based transgene control"]
    )
    # library details
    library_expression_control_id: Series[String] = Field(
        nullable=True, description="Library expression control ontology term ID."
    )
    library_expression_control_label: Series[String] = Field(
        nullable=True,
        description="Library expression control ontology term label.",
        isin=["constitutive transgene expression", "inducible transgene expression", "native promoter-driven transgene expression", "degradation domain-based transgene control"]
    )
    library_name: Series[String] = Field(
        nullable=True,
        description="Name of the perturbation library. Example: Bassik Human CRISPR Knockout Library",
    )
    library_uri: Series[String] = Field(
        nullable=True, description="URI/accession of the perturbation library."
    )
    library_format_id: Series[String] = Field(
        nullable=True, description="Perturbation library format ontology term ID."
    )
    library_format_label: Series[String] = Field(
        nullable=True,
        description="Perturbation library format ontology term label.",
        isin=["pooled", "arrayed", "arrayed|pooled", "in vivo"],
    )
    library_scope_id: Series[String] = Field(
        nullable=True, description="Perturbation library scope ontology term ID."
    )
    library_scope_label: Series[String] = Field(
        nullable=True,
        description="Perturbation library scope ontology term label.",
        isin=["focused", "genome-wide"],
    )
    library_perturbation_type_id: Series[String] = Field(
        nullable=True, description="Ontology term ID for the library perturbation type."
    )
    library_perturbation_type_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the library perturbation type.",
        isin=["knockout", "inhibition", "activation", "base editing", "prime editing", "mutagenesis"],
    )
    library_manufacturer: Series[String] = Field(
        nullable=True,
        description="Name of the library manufacturer/vendor/origin lab. Example: Bassik",
    )
    library_lentiviral_generation: Series[String] = Field(
        nullable=True,
        description="Generation number of the lentiviral library. Example: 3",
    )
    library_grnas_per_target: Series[String] = Field(
        nullable=True, description="Number of gRNAs per target. Example: 4, 5-7"
    )
    library_total_grnas: Series[String] = Field(
        nullable=True,
        coerce=True,
        description="Total number of gRNAs in the library. Example: 20,000",
    )
    library_total_variants: Int64 = Field(
        nullable=True,
        ge=0,
        description="Only for MAVE studies; Total number of variants in the library. Example: 5,000",
    )
    # assay details
    readout_dimensionality_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the dimensionality of the readout assay.",
    )
    readout_dimensionality_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the dimensionality of the readout assay.",
        isin=["single-dimensional assay", "high-dimensional assay"],
    )
    readout_type_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the type of the readout assay.",
    )
    readout_type_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the type of the readout assay.",
        isin=["transcriptomic", "proteomic", "phenotypic"],
    )
    readout_technology_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the technology used in the readout assay.",
    )
    readout_technology_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the technology used in the readout assay.",
        isin=["single-cell rna-seq", "population growth assay", "flow cytometry"],
    )
    method_name_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the method name used in the readout assay.",
    )
    method_name_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the method name used in the readout assay.",
        isin=["Perturb-seq", "scRNA-seq", "proliferation CRISPR screen", "DMS-TileSeq", "DMS-BarSeq", "Joined and refined DMS-BarSeq and DMS-TileSeq", "Combined DMS-BarSeq and DMS-TileSeq"],
    )
    method_uri: Series[String] = Field(
        nullable=True,
        description="URI associated with the method used in the readout assay.",
    )
    sequencing_library_kit_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing library kit.",
    )
    sequencing_library_kit_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing library kit.",
        isin=[
            "10x Genomics Chromium GEM-X Single Cell 5-prime kit v3",
            "10x Genomics Single Cell 3-prime",
            "Nextera XT DNA Library Preparation Kit",
        ],
    )
    sequencing_platform_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing platform.",
    )
    sequencing_platform_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing platform.",
        isin=["Illumina NovaSeq X Plus", "Illumina HiSeq 4000", "Illumina HiSeq 2500", "Illumina HiSeq 2000", "Illumina NovaSeq 6000", "Illumina NextSeq 500"],
    )
    sequencing_strategy_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing strategy.",
    )
    sequencing_strategy_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing strategy.",
        isin=["barcode sequencing", "direct sequencing", "barcode sequencing|direct sequencing"],
    )
    software_counts_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID for the software used for generating counts.",
    )
    software_counts_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the software used for generating counts.",
        isin=["custom", "MaGeCK", "CellRanger", "Drop-seq Tools"],
    )
    software_analysis_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID for the software used for analysis.",
    )
    software_analysis_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the software used for analysis.",
        isin=["custom", "MAGeCK", "Achilles"],
    )
    score_interpretation: Series[String] = Field(
        nullable=True,
        description="Interpretation of the perturbation effect score."
    )
    reference_genome_id: Series[String] = Field(
        nullable=True, description="Ontology term ID for the reference genome."
    )
    reference_genome_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the reference genome.",
        isin=["GRCh38", "GRCh37"],
    )
    # associated datasets
    associated_datasets: Series[String] = Field(
        nullable=True,
        coerce=True,
        description="List of associated datasets with each dataset having 'dataset_accession', 'dataset_uri', 'dataset_description', 'dataset_file_name' keys.",
    )

    class Config:
        strict = True
        # coerce = False
        ordered = True


# adata.var schema
class VarSchema(DataFrameModel):
    index: Index[str] = Field(
        nullable=False,
        unique=True,
        check_name=True,
        description="Unique identifier for each gene. Usually the Ensembl gene ID, or whatever unique IDs the dataset came with"
    )
    ensembl_gene_id: Series[str] = Field(
        nullable=True,
        str_startswith=("ENSG", "control"),
        description="Ensembl gene ID"
    )
    gene_symbol: Series[str] = Field(
        nullable=True,
        coerce=True,
        description="Gene symbol"
    )

    class Config:
        strict = "filter"
        coerce = True
        ordered = True
