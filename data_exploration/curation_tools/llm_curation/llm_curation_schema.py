"""Pydantic schemas used for metadata curation."""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class CurationSchema(BaseModel):

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    dataset_id_evidence: str | None = Field(
        default=None,
    )
    dataset_id: str = Field(
        default=...,
        description="Unique identifier for the dataset, follows the format <firstauthor_year>",
    )

    data_modality_evidence: str | None = Field(
        default=None,
    )
    data_modality: Literal["Perturb-seq", "CRISPR screen", "MAVE"] = Field(
        default=...,
        description="Data modality of the dataset.",
    )

    perturbation_type_label_evidence: str | None = Field(
        default=None,
    )
    perturbation_type_label: Literal["CRISPRn", "CRISPRi", "CRISPRa", "DMS"] = Field(
        default=...,
        description="Perturbation type ontology term label of the investigated sample.",
    )

    model_system_label_evidence: str | None = Field(
        default=None,
    )
    model_system_label: Literal[
        "cell_line", "primary_cell", "organoid", "yeast", "Other"
    ] = Field(
        default=...,
        description="Model system ontology term label of the investigated sample.",
    )

    species_evidence: str | None = Field(
        default=None,
    )
    species: Literal["Homo sapiens"] = Field(
        default=...,
        description="Species name of the investigated sample.",
    )

    tissue_label_evidence: str | None = Field(
        default=None,
    )
    tissue_label: str | None = Field(
        default=...,
        description="Tissue ontology term label of the investigated sample. Must be part of the UBERON ontology.",
    )

    cell_type_label_evidence: str | None = Field(
        default=None,
    )
    cell_type_label: str | None = Field(
        default=...,
        description="Cell type ontology term label of the investigated sample. Must be part of the Cell Ontology (CL).",
    )

    cell_line_label_evidence: str | None = Field(
        default=None,
    )
    cell_line_label: str | None = Field(
        default=...,
        description="Cell line ontology term label of the investigated sample. Must be part of the Cell Line Ontology (CLO).",
    )

    sex_label_evidence: str | None = Field(
        default=None,
    )
    sex_label: Literal["female", "male", "mixed", "unknown"] | None = Field(
        default=...,
        description="Sex ontology term label of the investigated sample.",
    )

    developmental_stage_label_evidence: str | None = Field(
        default=None,
    )
    developmental_stage_label: (
        Literal[
            "embryonic",
            "fetal",
            "neonatal",
            "child",
            "adolescent",
            "adult",
            "senior adult",
        ]
        | None
    ) = Field(
        default=...,
        description="Developmental stage ontology term label of the investigated sample.",
    )

    disease_label_evidence: str | None = Field(
        default=None,
    )
    disease_label: str | None = Field(
        default=...,
        description="Disease ontology term label of the investigated sample. Must be part of the MONDO ontology.",
    )

    study_title_evidence: str | None = Field(
        default=None,
    )
    study_title: str = Field(
        default=...,
        description="Title of the study/publication.",
    )

    study_uri_evidence: str | None = Field(
        default=None,
    )
    study_uri: str = Field(
        default=...,
        description="URI/DOI of the study/publication.",
    )

    study_year_evidence: str | None = Field(
        default=None,
    )
    study_year: int = Field(
        default=...,
        description="Publication year of the study/publication.",
    )

    first_author_evidence: str | None = Field(
        default=None,
    )
    first_author: str | None = Field(
        default=...,
        description="Full name of the first author of the study/publication.",
    )

    last_author_evidence: str | None = Field(
        default=None,
    )
    last_author: str | None = Field(
        default=...,
        description="Full name of the last author of the study/publication.",
    )

    experiment_title_evidence: str | None = Field(
        default=None,
    )
    experiment_title: str = Field(
        default=...,
        description="Title of the experiment.",
    )

    experiment_summary_evidence: str | None = Field(
        default=None,
    )
    experiment_summary: str | None = Field(
        default=...,
        description="Summary of the experiment.",
    )

    library_generation_type_label_evidence: str | None = Field(
        default=None,
    )
    library_generation_type_label: (
        Literal["endogenous genetic perturbation method", "exogenous genetic perturbation method"]
        | None
    ) = Field(
        default=...,
        description="Library generation type ontology term label, defined in EFO under parent term EFO:0022867 (genetic perturbation)",
    )

    library_generation_method_label_evidence: str | None = Field(
        default=None,
    )
    library_generation_method_label: (
        Literal[
            "doped oligo synthesis",
            "error-prone PCR",
            "microarray synthesis",
            "nicking mutagenesis",
            "oligo-directed mutagenic PCR",
            "site-directed mutagenesis",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Library generation method ontology term label, defined in EFO under parent term EFO:0022868/EFO:0022869 (Endogenous/Exogenous genetic perturbation method)",
    )

    enzyme_delivery_method_label_evidence: str | None = Field(
        default=None,
    )
    enzyme_delivery_method_label: (
        Literal[
            "lipofection",
            "nucleofection",
            "retrovirus transduction",
            "lentivirus transduction",
            "transformation",
            "nanoparticle-mediated transfection",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Enzyme delivery method ontology term label.",
    )

    library_delivery_method_label_evidence: str | None = Field(
        default=None,
    )
    library_delivery_method_label: (
        Literal[
            "lipofection",
            "nucleofection",
            "retrovirus transduction",
            "lentivirus transduction",
            "transformation",
            "nanoparticle-mediated transfection",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Library delivery method ontology term label.",
    )

    enzyme_integration_state_label_evidence: str | None = Field(
        default=None,
    )
    enzyme_integration_state_label: (
        Literal[
            "random locus integration",
            "targeted locus integration",
            "native locus replacement",
            "non-integrative transgene expression",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Enzyme integration state ontology term label.",
    )

    library_integration_state_label_evidence: str | None = Field(
        default=None,
    )
    library_integration_state_label: (
        Literal[
            "random locus integration",
            "targeted locus integration",
            "native locus replacement",
            "non-integrative transgene expression",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Library integration state ontology term label.",
    )

    enzyme_expression_control_label_evidence: str | None = Field(
        default=None,
    )
    enzyme_expression_control_label: (
        Literal[
            "constitutive transgene expression",
            "inducible transgene expression",
            "native promoter-driven transgene expression",
            "degradation domain-based transgene control",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Enzyme expression control ontology term label.",
    )

    library_expression_control_label_evidence: str | None = Field(
        default=None,
    )
    library_expression_control_label: (
        Literal[
            "constitutive transgene expression",
            "inducible transgene expression",
            "native promoter-driven transgene expression",
            "degradation domain-based transgene control",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Library expression control ontology term label.",
    )

    library_name_evidence: str | None = Field(
        default=None,
    )
    library_name: str | None = Field(
        default=...,
        description="Name of the perturbation library. Example: Bassik Human CRISPR Knockout Library",
    )

    library_uri_evidence: str | None = Field(
        default=None,
    )
    library_uri: str | None = Field(
        default=...,
        description="URI/accession of the perturbation library.",
    )

    library_format_label_evidence: str | None = Field(
        default=None,
    )
    library_format_label: (
        Literal["pooled", "arrayed", "arrayed|pooled", "in vivo"] | None
    ) = Field(
        default=...,
        description="Perturbation library format ontology term label.",
    )

    library_scope_label_evidence: str | None = Field(
        default=None,
    )
    library_scope_label: Literal["focused", "genome-wide"] | None = Field(
        default=...,
        description="Perturbation library scope ontology term label.",
    )

    library_perturbation_type_label_evidence: str | None = Field(
        default=None,
    )
    library_perturbation_type_label: (
        Literal[
            "knockout",
            "inhibition",
            "activation",
            "base editing",
            "prime editing",
            "mutagenesis",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Ontology term label for the library perturbation type.",
    )

    library_manufacturer_evidence: str | None = Field(
        default=None,
    )
    library_manufacturer: str | None = Field(
        default=...,
        description="Name of the library manufacturer/vendor/origin lab. Example: Bassik",
    )

    library_lentiviral_generation_evidence: str | None = Field(
        default=None,
    )
    library_lentiviral_generation: str | None = Field(
        default=...,
        description="Generation number of the lentiviral library. Example: 3",
    )

    library_grnas_per_target_evidence: str | None = Field(
        default=None,
    )
    library_grnas_per_target: str | None = Field(
        default=...,
        description="Number of gRNAs per target. Example: 4, 5-7",
    )

    library_total_grnas_evidence: str | None = Field(
        default=None,
    )
    library_total_grnas: str | None = Field(
        default=...,
        description="Total number of gRNAs in the library. Example: 20,000",
    )

    library_total_variants_evidence: str | None = Field(
        default=None,
    )
    library_total_variants: int | None = Field(
        default=...,
        description="Only for MAVE studies; Total number of variants in the library. Example: 5,000",
        ge=0,
    )

    readout_dimensionality_label_evidence: str | None = Field(
        default=None,
    )
    readout_dimensionality_label: (
        Literal["single-dimensional assay", "high-dimensional assay"] | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the dimensionality of the readout assay.",
    )

    readout_type_label_evidence: str | None = Field(
        default=None,
    )
    readout_type_label: (
        Literal["transcriptomic", "proteomic", "phenotypic", "Other"] | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the type of the readout assay.",
    )

    readout_technology_label_evidence: str | None = Field(
        default=None,
    )
    readout_technology_label: (
        Literal[
            "single-cell rna-seq", "population growth assay", "flow cytometry", "Other"
        ]
        | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the technology used in the readout assay.",
    )

    readout_technology_label_evidence: str | None = Field(
        default=None,
    )
    readout_measurment_label_evidence: str | None = Field(
        default=None,
    )
    readout_measurment_label: (
        Literal[
            "surface protein expression",
            "cell viability",
            "gene expression",
            "protein abundance",
            "ligand binding",
            "cell proliferation",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the measurement type of the readout assay.",
    )

    method_name_label_evidence: str | None = Field(
        default=None,
    )
    method_name_label: (
        Literal[
            "Perturb-seq",
            "Perturb-CITE-seq",
            "scRNA-seq",
            "proliferation CRISPR screen",
            "DMS-TileSeq",
            "DMS-BarSeq",
            "Joined and refined DMS-BarSeq and DMS-TileSeq",
            "Combined DMS-BarSeq and DMS-TileSeq",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the method name used in the readout assay.",
    )

    method_uri_evidence: str | None = Field(
        default=None,
    )
    method_uri: str | None = Field(
        default=...,
        description="URI associated with the method used in the readout assay.",
    )

    sequencing_library_kit_label_evidence: str | None = Field(
        default=None,
    )
    sequencing_library_kit_label: (
        Literal[
            "10x Genomics Chromium GEM-X Single Cell 5-prime kit v3",
            "10x Genomics Chromium Next GEM Single Cell 5-prime HT Kit v2",
            "10x Genomics Single Cell 3-prime",
            "10x Genomics Single Cell 3-prime v2",
            "10x Genomics Single Cell 3-prime v3",
            "Nextera XT DNA Library Preparation Kit",
            "GEM-X Flex Gene Expression Human n-plex kit",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the sequencing library kit.",
    )

    sequencing_platform_label_evidence: str | None = Field(
        default=None,
    )
    sequencing_platform_label: (
        Literal[
            "Illumina NovaSeq X",
            "Illumina NovaSeq X Plus",
            "Illumina HiSeq 4000",
            "Illumina HiSeq 2500",
            "Illumina HiSeq 2000",
            "Illumina NovaSeq 6000",
            "Illumina NextSeq 500",
            "Ultima Genomics UG100",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the sequencing platform.",
    )

    sequencing_strategy_label_evidence: str | None = Field(
        default=None,
    )
    sequencing_strategy_label: (
        Literal[
            "barcode sequencing",
            "direct sequencing",
            "barcode sequencing|direct sequencing",
            "Other",
        ]
        | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the sequencing strategy.",
    )

    software_counts_label_evidence: str | None = Field(
        default=None,
    )
    software_counts_label: (
        Literal["custom", "MaGeCK", "CellRanger", "Drop-seq Tools", "Other"] | None
    ) = Field(
        default=...,
        description="Ontology term label for the software used for generating counts.",
    )

    software_analysis_label_evidence: str | None = Field(
        default=None,
    )
    software_analysis_label: (
        Literal[
            "custom", "MAGeCK", "Achilles", "TRADE", "Seurat", "MAST", "scanpy", "Other"
        ]
        | None
    ) = Field(
        default=...,
        description="Ontology term label for the software used for analysis.",
    )

    reference_genome_label_evidence: str | None = Field(
        default=None,
    )
    reference_genome_label: Literal["GRCh38", "GRCh37", "Other"] | None = Field(
        default=...,
        description="Ontology term label for the reference genome.",
    )

    associated_datasets_evidence: str | None = Field(
        default=None,
    )
    associated_datasets: str | None = Field(
        default=...,
        description="List of associated datasets with each dataset having 'dataset_accession', 'dataset_uri', 'dataset_description', 'dataset_file_name' keys.",
    )

    license_label_evidence: str | None = Field(
        default=None,
    )
    license_label: Literal[
        "CC0", "CC BY", "CC BY-SA", "CC BY-NC", "CC BY-ND", "Other"
    ] = Field(
        default=...,
        description="License type for data usage and distribution. Should be one of the terms from under SWO:0000002 (license).",
    )
