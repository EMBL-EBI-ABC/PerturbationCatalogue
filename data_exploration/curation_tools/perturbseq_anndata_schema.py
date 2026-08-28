import pandas as pd
from pandera.pandas import Field, DataFrameModel, dataframe_check
from pandera.typing import Series, Index, String, Int64, Float32
from pathlib import Path


class ObsSchema(DataFrameModel):
    dataset_id: Series[String] = Field(
        nullable=False,
        description="Unique identifier for the dataset, follows the format <firstauthor_year>",
    )
    sample_id: Series[String] = Field(
        nullable=False, coerce=True, description="Unique identifier for the sample."
    )
    cell_barcode: Series[String] = Field(
        nullable=True,
        coerce=True,
        description="Unique cell barcode.",
    )

    data_modality: Series[String] = Field(
        nullable=False,
        description="Data modality of the dataset.",
        isin=["Perturb-seq", "CRISPR screen", "MAVE"],
    )
    significant: Series[String] = Field(
        nullable=True,
        description="Indicates whether the perturbation had a significant effect.",
        coerce=True,
        isin=["True", "False"],
    )
    significance_criteria: Series[String] = Field(
        nullable=True,
        description="Criteria used to determine significance, e.g., FDR < 0.05.",
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
        coerce=True,
        description="Numeric encoding of the chromosome of the perturbed target. Required for data partitioning in BigQuery.",
    )
    perturbed_target_number: Series[Int64] = Field(
        nullable=False,
        ge=0,
        coerce=True,
        description="Number of perturbed targets in the samples.",
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
    timepoint_post_transfection: Series[String] = Field(
        nullable=True,
        regex=r"^P\d+DT\d{1,2}H\d{1,2}M\d{1,2}S$",
        description="Timepoint of the investigated sample in ISO 8601 format, starting from the time of library transfection. Example: P1DT12H30M15S",
    )
    differentiation_timepoint: Series[String] = Field(
        nullable=True,
        regex=r"^P\d+DT\d{1,2}H\d{1,2}M\d{1,2}S$",
        description="Differentiation timepoint of the investigated sample in ISO 8601 format, starting from the moment the induction of differentiation began. Example: P1DT12H30M15S",
    )
    experimental_timepoint: Series[String] = Field(
        nullable=True,
        regex=r"^P\d+DT\d{1,2}H\d{1,2}M\d{1,2}S$",
        description="Experimental timepoint of the investigated sample in ISO 8601 format. Example: P1DT12H30M15S",
    )
    treatment_type_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label describing the treatment type.",
        isin=[
            "untreated control", # NCIT:C184729
            "scrambled control oligonucleotide", # XCO:0001141
            "culture medium", # BAO:0000114
            "chemical entity", # CHEBI:24431
            "protein", # BAO:0000175
            "protein complex", # BAO:0002554
            "peptide", # BAO:0000325
            "antibody", # BAO:0000502
            "lipid", # BAO:0000171
            "PNA", # BAO:0000226
            "DNA", # BAO:0000269
            "RNA", # BAO:0000270
            "mRNA", # BAO:0000274
            "rRNA", # BAO:0000275
            "tRNA", # BAO:0000276
            "cDNA", # BAO:0000315
            "genomic DNA", # BAO:0000316
            "plasmid DNA", # BAO:0000317
            "miRNA", # BAO:0000322
            "shRNA", # BAO:0000323
            "siRNA", # BAO:0000324
            "LNA", # BAO:0000412
            "RNA aptamer", # BAO:0000496
            "riboswitch", # BAO:0000498
            "esiRNA", # BAO:0000544
        ],
    )
    treatment_type_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID for the treatment type.",
        isin=[
            "NCIT:C184729",  # untreated control
            "XCO:0001141",  # scrambled control oligonucleotide
            "BAO:0000114",  # culture medium
            "CHEBI:24431",  # chemical entity
            "BAO:0000175",  # protein
            "BAO:0002554",  # protein complex
            "BAO:0000325",  # peptide
            "BAO:0000502",  # antibody
            "BAO:0000171",  # lipid
            "BAO:0000226",  # peptide nucleic acid (PNA)
            "BAO:0000269",  # DNA
            "BAO:0000270",  # RNA
            "BAO:0000274",  # messenger RNA (mRNA)
            "BAO:0000275",  # ribosomal RNA (rRNA)
            "BAO:0000276",  # transfer RNA (tRNA)
            "BAO:0000315",  # complementary DNA (cDNA)
            "BAO:0000316",  # genomic DNA
            "BAO:0000317",  # plasmid DNA
            "BAO:0000322",  # microRNA (miRNA)
            "BAO:0000323",  # short hairpin RNA (shRNA)
            "BAO:0000324",  # small interfering RNA (siRNA)
            "BAO:0000412",  # locked nucleic acid (LNA)
            "BAO:0000496",  # RNA aptamer
            "BAO:0000498",  # riboswitch
            "BAO:0000544",  # endoribonuclease-prepared siRNA (esiRNA)
        ],
    )
    treatment_label: Series[String] = Field(
        nullable=True,
        description="Treatment/compound ontology term label used to stimulate the investigated sample. ChEMBL compound label for chemical entities. Use 'untreated control' for untreated samples where other samples were treated.",
    )
    treatment_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Treatment/compound ontology term ID used to stimulate the investigated sample. ChEMBL compound ID.",
    )
    treatment_dose: Series[Float32] = Field(
        nullable=True,
        coerce=True,
        description="Treatment/compound dose used to stimulate the investigated sample.",
    )
    treatment_unit: Series[String] = Field(
        nullable=True,
        description="Treatment/compound unit used to stimulate the investigated sample. Use 'u' for micro (e.g., 'uM' instead of 'μM').",
        isin=[
            # Concentration (molar)
            "pM",
            "nM",
            "uM",
            "mM",
            "M",
            # Mass/volume concentration
            "pg/mL",
            "ng/mL",
            "ug/mL",
            "mg/mL",
            "g/mL",
            # Mass/mass concentration
            "pg/kg",
            "ug/kg",
            "mg/kg",
            "g/kg",
            # Count-based
            "cells/uL",
            "cells/mL",
            "MOI",
            # Volume
            "uL",
            "mL",
            # Other
            "%",
            "IU/mL",
        ],
    )
    technical_replicate: Series[String] = Field(
        nullable=True, description="Technical replicate id."
    )
    biological_replicate: Series[String] = Field(
        nullable=True, description="Biological replicate id."
    )
    # model system details
    model_system_label: Series[String] = Field(
        nullable=False,
        description="Model system ontology term label of the investigated sample.",
        isin=[
            "cell_line", # CLO:0000031
            "primary_cell", # BAO:0000239
            "organoid", # NCIT:C172259
            "yeast", # NCIT:C19617
            "bacteria", # NCIT:C19167
            "bacteriophage", # NCIT:C14188
            "animal_model", # NCIT:C71164
            "cell_free_system", # mesh:D002474
            "Other",
        ],
    )
    
    model_system_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Model system ontology term ID of the investigated sample.",
        isin=[
            "CLO:0000031", # cell_line
            "BAO:0000239", # primary_cell
            "NCIT:C172259", # organoid
            "NCIT:C19617", # yeast
            "NCIT:C19167", # bacteria
            "NCIT:C14188", # bacteriophage
            "NCIT:C71164", # animal_model
            "mesh:D002474", # cell_free_system
        ],
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
        isin=[
            "female", # PATO:0000383
            "male", # PATO:0000384
            "hermaphrodite", # PATO:0001340
            "unknown"
        ],
    )
    sex_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Sex ontology term ID of the investigated sample.",
        isin=[
            "PATO:0000383", # female
            "PATO:0000384", # male
            "PATO:0001340" # hermaphrodite
        ],
    )
    developmental_stage_label: Series[String] = Field(
        nullable=True,
        description="Developmental stage ontology term label of the investigated sample.",
        isin=[
            "embryonic", # HsapDv:0000002
            "fetal", # HsapDv:0000037
            "neonatal", # HsapDv:0000262
            "child", # HsapDv:0000265 1-4 yo
            "juvenile", # HsapDv:0000271 5-14 yo
            "adult", # HsapDv:0000258 16-59 yo
            "elderly", # HsapDv:0000227 60+ yo
        ],
    )
    developmental_stage_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Developmental stage ontology term ID of the investigated sample.",
        isin=[
            "HsapDv:0000002", # embryonic
            "HsapDv:0000037", # fetal
            "HsapDv:0000262", # neonatal
            "HsapDv:0000265", # child
            "HsapDv:0000271", # juvenile
            "HsapDv:0000258", # adult
            "HsapDv:0000227", # elderly
        ],
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
    library_generation_type_label: Series[String] = Field(
        nullable=True,
        description="Library generation type ontology term label, defined in EFO under parent term EFO:0022867 (genetic perturbation)",
        isin=[
            "endogenous genetic perturbation method", # EFO:0022868
            "exogenous genetic perturbation method", # EFO:0022869
        ],
    )
    library_generation_type_id: Series[String] = Field(
        nullable=True,
        description="Library generation type ontology term ID, defined in EFO under parent term EFO:0022867 (genetic perturbation)",
        isin=[
            "EFO:0022868", # endogenous genetic perturbation method
            "EFO:0022869", # exogenous genetic perturbation method
        ],
    )
    library_generation_method_label: Series[String] = Field(
        nullable=True,
        description="Library generation method ontology term label, defined in EFO under parent term EFO:0022868/EFO:0022869 (Endogenous/Exogenous genetic perturbation method)",
        isin=[
            "doped oligo synthesis",
            "error-prone PCR",
            "microarray synthesis",
            "nicking mutagenesis",
            "oligo-directed mutagenic PCR",
            "site-directed mutagenesis",
            "silicon microarray synthesis",
            "POPCode mutagenesis",
            "insertional mutagenesis",
            "solid-phase oligonucleotide synthesis",
            "multiplexed site-directed mutagenesis",
            "microchip-based massive parallel oligo synthesis",
            "Other",
        ],
    )
    library_generation_method_id: Series[String] = Field(
        nullable=True,
        description="Library generation method ontology term ID, defined in EFO under parent term EFO:0022868/EFO:0022869 (Endogenous/Exogenous genetic perturbation method)",
    )
    enzyme_delivery_method_label: Series[String] = Field(
        nullable=True,
        description="Enzyme delivery method ontology term label.",
        isin=[
            "lipofection", # EFO:0920076
            "nucleofection", # EFO:0920075
            "electroporation", # EFO:0920074
            "adeno-associated virus transduction", # EFO:0920071
            "adenovirus transduction", # EFO:0920070
            "retrovirus transduction", # EFO:0920068
            "lentivirus transduction", # EFO:0920067
            "nanoparticle-mediated transfection", # EFO:0920077
            "transformation",
            "chemical-mediated transfection",
            "hydrodynamic injection",
            "molecular cloning",
            "influenza A virus infection",
            "Other",
        ],
    )
    enzyme_delivery_method_id: Series[String] = Field(
        nullable=True,
        description="Enzyme delivery method ontology term ID.",
    )
    library_delivery_method_label: Series[String] = Field(
        nullable=True,
        description="Library delivery method ontology term label.",
        isin=[
            "lipofection", # EFO:0920076
            "nucleofection", # EFO:0920075
            "electroporation", # EFO:0920074
            "adeno-associated virus transduction", # EFO:0920071
            "adenovirus transduction", # EFO:0920070
            "retrovirus transduction", # EFO:0920068
            "lentivirus transduction", # EFO:0920067
            "nanoparticle-mediated transfection", # EFO:0920077
            "transformation",
            "chemical-mediated transfection",
            "hydrodynamic injection",
            "molecular cloning",
            "influenza A virus infection",
            "Other",
        ],
    )
    library_delivery_method_id: Series[String] = Field(
        nullable=True, description="Library delivery method ontology term ID."
    )
    enzyme_integration_state_label: Series[String] = Field(
        nullable=True,
        description="Enzyme integration state ontology term label.",
        isin=[
            "random locus integration", # EFO:0920082
            "targeted locus integration", # EFO:0920083
            "native locus replacement", # EFO:0920084
            "non-integrative transgene expression", # EFO:0920085
            "bacteriophage genome integration",
            "Other",
        ],
    )
    enzyme_integration_state_id: Series[String] = Field(
        nullable=True, description="Enzyme integration state ontology term ID.",
        isin=[
            "EFO:0920082", # random locus integration
            "EFO:0920083", # targeted locus integration
            "EFO:0920084", # native locus replacement
            "EFO:0920085", # non-integrative transgene expression
        ],
    )
    library_integration_state_label: Series[String] = Field(
        nullable=True,
        description="Library integration state ontology term label.",
        isin=[
            "random locus integration", # EFO:0920082
            "targeted locus integration", # EFO:0920083
            "native locus replacement", # EFO:0920084
            "non-integrative transgene expression", # EFO:0920085
            "bacteriophage genome integration",
            "Other",
        ],
    )
    library_integration_state_id: Series[String] = Field(
        nullable=True, description="Library integration state ontology term ID.",
        isin=[
            "EFO:0920082", # random locus integration
            "EFO:0920083", # targeted locus integration
            "EFO:0920084", # native locus replacement
            "EFO:0920085", # non-integrative transgene expression
        ],
    )
    enzyme_expression_control_label: Series[String] = Field(
        nullable=True,
        description="Enzyme expression control ontology term label.",
        isin=[
            "constitutive transgene expression",
            "inducible transgene expression",
            "native promoter-driven transgene expression",
            "degradation domain-based transgene control",
            "transient transgene expression",
            "minimal promoter-driven transgene expression",
            "Other",
        ],
    )
    enzyme_expression_control_id: Series[String] = Field(
        nullable=True, description="Enzyme expression control ontology term ID."
    )
    # library details
    library_expression_control_label: Series[String] = Field(
        nullable=True,
        description="Library expression control ontology term label.",
        isin=[
            "constitutive transgene expression",
            "inducible transgene expression",
            "native promoter-driven transgene expression",
            "degradation domain-based transgene control",
            "transient transgene expression",
            "minimal promoter-driven transgene expression",
            "Other",
        ],
    )
    library_expression_control_id: Series[String] = Field(
        nullable=True, description="Library expression control ontology term ID."
    )
    library_name: Series[String] = Field(
        nullable=True,
        description="Name of the perturbation library. Example: Bassik Human CRISPR Knockout Library",
    )
    library_uri: Series[String] = Field(
        nullable=True, description="URI/accession of the perturbation library."
    )
    library_format_label: Series[String] = Field(
        nullable=True,
        description="Perturbation library format ontology term label.",
        isin=["pooled", "arrayed", "arrayed|pooled", "in vivo"],
    )
    library_format_id: Series[String] = Field(
        nullable=True, description="Perturbation library format ontology term ID."
    )
    library_scope_label: Series[String] = Field(
        nullable=True,
        description="Perturbation library scope ontology term label.",
        isin=["focused", "genome-wide"],
    )
    library_scope_id: Series[String] = Field(
        nullable=True, description="Perturbation library scope ontology term ID."
    )
    library_perturbation_type_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the library perturbation type.",
        isin=[
            "knockout", # EFO:0000506
            "inhibition", # INO:0000085
            "activation", # INO:0000075
            "base editing", # EFO:0022873
            "prime editing", # EFO:0022872
            "mutagenesis", # NCIT:C17376
            "Other",
        ],
    )
    library_perturbation_type_id: Series[String] = Field(
        nullable=True, description="Ontology term ID for the library perturbation type.",
        isin=[
            "EFO:0000506", # knockout
            "INO:0000085", # inhibition
            "INO:0000075", # activation
            "EFO:0022873", # base editing
            "EFO:0022872", # prime editing
            "NCIT:C17376", # mutagenesis
            "Other",
        ],
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
    readout_dimensionality_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the dimensionality of the readout assay.",
        isin=["single-dimensional assay", "high-dimensional assay"],
    )
    readout_dimensionality_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the dimensionality of the readout assay.",
    )
    readout_type_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the type of the readout assay.",
        isin=[
            "transcriptomic", # EFO:0001032
            "proteomic", # EFO:0000746
            "phenotypic", # EFO:0920062
            "Other"
        ],
    )
    readout_type_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the type of the readout assay.",
        isin=[
            "EFO:0001032", # transcriptomic
            "EFO:0000746", # proteomic
            "EFO:0920062", # phenotypic
            "Other",
        ],
    )
    readout_technology_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the technology used in the readout assay.",
        isin=[
            "single-cell rna-seq", # EFO:0008913
            "population growth assay", # EFO:0002907
            "flow cytometry", # BAO:0000005
            "high-throughput dna sequencing", # EFO:0002693
            "patch-clamp electrophysiology", # EFO:0022948
            "fluorometry", # mesh:D005470
            "Other",
        ],
    )
    readout_technology_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the technology used in the readout assay.",
        isin=[
            "EFO:0008913", # single-cell rna-seq
            "EFO:0002907", # population growth assay
            "BAO:0000005", # flow cytometry
            "EFO:0002693", # high-throughput dna sequencing
            "EFO:0022948", # patch-clamp electrophysiology
            "mesh:D005470", # fluorometry
            "Other",
        ],
    )
    readout_measurement_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the measurement type of the readout assay.",
        isin=[
            "protein abundance", # BAO:0010252
            "protein stability", # BAO:0002804
            "protein activity", # APO:0000022
            "protein ubiquitination", # GO:0016567
            "surface protein expression",
            "cell viability", # PATO:0000169
            "cell proliferation", # BAO:0002805
            "gene expression", # BAO:0002785
            "RNA splicing", # BAO:0003000
            "DNA repair", # GO:0006281
            "ligand binding", # NCIT:C178030
            "ion channel activity", # BAO:0002997
            "fluorescence", # BAO:0000363
            "Other",
        ],
    )
    readout_measurement_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the measurement type of the readout assay.",
        isin=[
            "BAO:0010252", # protein abundance
            "BAO:0002804", # protein stability
            "APO:0000022", # protein activity
            "GO:0016567", # protein ubiquitination
            "PATO:0000169", # cell viability
            "BAO:0002805", # cell proliferation
            "BAO:0002785", # gene expression
            "BAO:0003000", # RNA splicing
            "GO:0006281", # DNA repair
            "NCIT:C178030", # ligand binding
            "BAO:0002997", # ion channel activity
            "BAO:0000363", # fluorescence
            "Other",
        ],
    )
    method_name_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the method name used in the readout assay.",
        isin=[
            "Perturb-seq", # EFO:0008860
            "scRNA-seq", # EFO:0008913
            "Perturb-CITE-seq",
            "proliferation CRISPR screen",
            "yeast surface display", # MI:0115
            "yeast one-hybrid assay", # OBI:0001681
            "yeast two-hybrid assay", # BAO:0002494
            "bacterial two-hybrid assay", # OBI:0001682
            "mammalian two-hybrid assay", # BAO:0002493
            "phage display", # MI:0084
            "abundance protein fragment complementation assay", # MI:0090
            "pooled growth competition assay", # EFO:0002907
            "massively parallel reporter assay", # EFO:0008822
            "patch-clamp electrophysiology", # EFO:0022948
            "DMS-TileSeq",
            "DMS-BarSeq",
            "Joined and refined DMS-BarSeq and DMS-TileSeq",
            "Combined DMS-BarSeq and DMS-TileSeq",
            "flow cytometry-based sequencing assay",
            "MITE",
            "VAMP-seq",
            "saturation genome editing",
            "saturation prime editing",
            "polysome profiling",
            "Other",
        ],
    )
    method_name_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the method name used in the readout assay.",
        isin=[
            "EFO:0008860", # Perturb-seq
            "EFO:0008913", # scRNA-seq
            "EFO:0002907", # pooled growth competition assay
            "EFO:0008822", # massively parallel reporter assay
            "MI:0115", # yeast surface display
            "OBI:0001682", # bacterial two-hybrid assay
            "BAO:0002493", # mammalian two-hybrid assay
            "OBI:0001681", # yeast one-hybrid assay
            "MI:0084", # phage display
            "BAO:0002494", # yeast two-hybrid assay
            "EFO:0022948", # patch-clamp electrophysiology
            "MI:0090", # abundance protein fragment complementation assay
        ],
    )
    method_uri: Series[String] = Field(
        nullable=True,
        description="URI associated with the method used in the readout assay.",
    )
    sequencing_library_kit_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing library kit.",
        isin=[
            "10x Genomics Single Cell 3-prime v1", # EFO:0009901
            "10x Genomics Single Cell 3-prime v2", # EFO:0009899
            "10x Genomics Single Cell 3-prime v3", # EFO:0009922
            "10x Genomics Single Cell 3-prime v3.1", # EFO:0022980
            "10x Genomics Chromium GEM-X Flex v1", # EFO:0920088
            "10x Genomics Chromium GEM-X Single Cell 5-prime kit v3",
            "10x Genomics Chromium Next GEM Single Cell 5-prime HT Kit v2",
            "Nextera XT DNA Library Preparation Kit",
            "Parse Biosciences Evercode Whole Transcriptome Mega v1 kit",
            "TruSeq Nano DNA Library Prep Kit",
            "Ovation Ultralow Library System",
            "custom PCR library preparation",
            "Nextera DNA Library Preparation Kit",
            "PacBio SMRTbell Template Prep Kit",
            "PacBio SMRTbell Template Prep Kit v1",
            "PacBio SMRTbell Template Prep Kit v2",
            "PacBio SMRTbell Template Prep Kit v3",
            "Beckman Coulter DTCS DNA sequencing kit",
            "Other",
        ],
    )
    sequencing_library_kit_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing library kit.",
        isin=[
            "EFO:0009901", # 10x Genomics Single Cell 3-prime v1
            "EFO:0009899", # 10x Genomics Single Cell 3-prime v2
            "EFO:0009922", # 10x Genomics Single Cell 3-prime v3
            "EFO:0022980", # 10x Genomics Single Cell 3-prime v3.1
            "EFO:0920088", # 10x Genomics Chromium GEM-X Flex v1
        ],
    )
    sequencing_platform_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing platform.",
        isin=[
            "Illumina Genome Analyzer", # "EFO:0004200"
            "Illumina Genome Analyzer II", # "EFO:0004201"
            "Illumina Genome Analyzer IIx", # "EFO:0004202"
            "Illumina HiSeq 2000", # "EFO:0004203"
            "Illumina HiSeq 1000", # "EFO:0004204"
            "Illumina MiSeq", # "EFO:0004205"
            "454 GS 20 sequencer", # "EFO:0004206"
            "454 GS sequencer", # "EFO:0004431"
            "454 GS FLX sequencer", # "EFO:0004432"
            "454 GS FLX Titanium sequencer", # "EFO:0004433"
            "454 GS Junior sequencer", # "EFO:0004434"
            "AB SOLiD System", # "EFO:0004435"
            "AB SOLiD 5500xl", # "EFO:0004436"
            "AB SOLiD PI System", # "EFO:0004437"
            "AB SOLiD 4 System", # "EFO:0004438"
            "AB SOLiD System 3.0", # "EFO:0004439"
            "AB SOLiD 5500", # "EFO:0004440"
            "AB SOLiD 4hq System", # "EFO:0004441"
            "AB SOLiD System 2.0", # "EFO:0004442"
            "Illumina HiSeq 4000", # "EFO:0008563"
            "Illumina HiSeq 3000", # "EFO:0008564"
            "Illumina HiSeq 2500", # "EFO:0008565"
            "Illumina NextSeq 550", # "EFO:0008566"
            "Illumina HiSeq X", # "EFO:0008567"
            "PacBio Sequel system", # "EFO:0008630"
            "PacBio RS II", # "EFO:0008631"
            "ONT MinION", # "EFO:0008632"
            "ONT GridION X5", # "EFO:0008633"
            "ONT PromethION", # "EFO:0008634"
            "Illumina iSeq 100", # "EFO:0008635"
            "Illumina MiniSeq", # "EFO:0008636"
            "Illumina NovaSeq 6000", # "EFO:0008637"
            "Illumina NextSeq 500", # "EFO:0009173"
            "Illumina NextSeq 1000", # "EFO:0010962"
            "Illumina NextSeq 2000", # "EFO:0010963"
            "Illumina HiSeq 1500", # "EFO:0011027"
            "Illumina NovaSeq X", # "EFO:0022840"
            "Illumina NovaSeq X Plus", # "EFO:0022841"
            "Singular G4", # "EFO:0022843"
            "PacBio Sequel II system", # "EFO:0700015"
            "BGI MGISEQ-2000", # "EFO:0700018"
            "ONT PromethION 2 Solo", # "EFO:0700019"
            "Ultima UG100", # "EFO:0920005"
            "PacBio Revio", # "EFO:0920006"
            "PacBio Onso", # "EFO:0920007"
            "Element Aviti", # "EFO:0920008"
            "Illumina MiSeq i100", # "EFO:0920010"
            "MGI DNBSEQ-T7", # "EFO:0920057"
            "Ion Torrent PGM", # GENEPIO:0100136
            "Other",
        ],
    )
    sequencing_platform_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing platform.",
        isin=[
            "EFO:0004200",  # "Illumina Genome Analyzer"
            "EFO:0004201",  # "Illumina Genome Analyzer II"
            "EFO:0004202",  # "Illumina Genome Analyzer IIx"
            "EFO:0004203",  # "Illumina HiSeq 2000"
            "EFO:0004204",  # "Illumina HiSeq 1000"
            "EFO:0004205",  # "Illumina MiSeq"
            "EFO:0004206",  # "454 GS 20 sequencer"
            "EFO:0004431",  # "454 GS sequencer"
            "EFO:0004432",  # "454 GS FLX sequencer"
            "EFO:0004433",  # "454 GS FLX Titanium sequencer"
            "EFO:0004434",  # "454 GS Junior sequencer"
            "EFO:0004435",  # "AB SOLiD System"
            "EFO:0004436",  # "AB SOLiD 5500xl"
            "EFO:0004437",  # "AB SOLiD PI System"
            "EFO:0004438",  # "AB SOLiD 4 System"
            "EFO:0004439",  # "AB SOLiD System 3.0"
            "EFO:0004440",  # "AB SOLiD 5500"
            "EFO:0004441",  # "AB SOLiD 4hq System"
            "EFO:0004442",  # "AB SOLiD System 2.0"
            "EFO:0008563",  # "Illumina HiSeq 4000"
            "EFO:0008564",  # "Illumina HiSeq 3000"
            "EFO:0008565",  # "Illumina HiSeq 2500"
            "EFO:0008566",  # "Illumina NextSeq 550"
            "EFO:0008567",  # "Illumina HiSeq X"
            "EFO:0008630",  # "PacBio Sequel system"
            "EFO:0008631",  # "PacBio RS II"
            "EFO:0008632",  # "ONT MinION"
            "EFO:0008633",  # "ONT GridION X5"
            "EFO:0008634",  # "ONT PromethION"
            "EFO:0008635",  # "Illumina iSeq 100"
            "EFO:0008636",  # "Illumina MiniSeq"
            "EFO:0008637",  # "Illumina NovaSeq 6000"
            "EFO:0009173",  # "Illumina NextSeq 500"
            "EFO:0010962",  # "Illumina NextSeq 1000"
            "EFO:0010963",  # "Illumina NextSeq 2000"
            "EFO:0011027",  # "Illumina HiSeq 1500"
            "EFO:0022840",  # "Illumina NovaSeq X"
            "EFO:0022841",  # "Illumina NovaSeq X Plus"
            "EFO:0022843",  # "Singular G4"
            "EFO:0700015",  # "PacBio Sequel II system"
            "EFO:0700018",  # "BGI MGISEQ-2000"
            "EFO:0700019",  # "ONT PromethION 2 Solo"
            "EFO:0920005",  # "Ultima UG100"
            "EFO:0920006",  # "PacBio Revio"
            "EFO:0920007",  # "PacBio Onso"
            "EFO:0920008",  # "Element Aviti"
            "EFO:0920010",  # "Illumina MiSeq i100"
            "EFO:0920057",  # "MGI DNBSEQ-T7"
            "GENEPIO:0100136",  # "Ion Torrent PGM"
        ]
    )
    sequencing_strategy_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing strategy.",
        isin=[
            "barcode sequencing",
            "direct sequencing", # NCIT:C116154
            "barcode sequencing|direct sequencing",
            "Other",
        ],
    )
    sequencing_strategy_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing strategy.",
        isin=[
            "NCIT:C116154", # direct sequencing
        ],
    )
    software_counts_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the software used for generating counts.",
        isin=[
            "custom",
            "MaGeCK",
            "CellRanger",
            "Drop-seq Tools",
            "Enrich2",
            "Enrich",
            "Novoalign",
            "TileSEQ Analysis Package",
            "DiMSum",
            "dms_tools",
            "dms_tools2",
            "dms_variants",
            "CRISPResso2",
            "mapmuts",
            "ORFcall",
            "Other",
        ],
    )
    software_counts_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID for the software used for generating counts.",
    )
    software_analysis_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the software used for analysis.",
        isin=[
            "custom",
            "MAGeCK",
            "Achilles",
            "TRADE",
            "Seurat",
            "MAST",
            "scanpy",
            "Enrich2",
            "DiMSum",
            "DESeq2",
            "dms_tools",
            "dms_tools2",
            "ALDEx2",
            "multidms",
            "dmsPipeline",
            "phydms",
            "Enrich",
            "Rosetta",
            "Other",
        ],
    )
    software_analysis_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID for the software used for analysis.",
    )
    score_interpretation: Series[String] = Field(
        nullable=True, description="Interpretation of the perturbation effect score."
    )
    reference_genome_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the reference genome.",
        isin=[
            "GRCh38",
            "GRCh37",
            "cDNA reference sequence",
            "mm9",
            "S288c",
            "hg19",
            "Wuhan-Hu-1",
            "non-standard reference sequence",
            "Other",
        ],
    )
    reference_genome_id: Series[String] = Field(
        nullable=True, description="Ontology term ID for the reference genome."
    )
    # associated datasets
    associated_datasets: Series[String] = Field(
        nullable=True,
        coerce=True,
        description="List of associated datasets with each dataset having 'dataset_accession', 'dataset_uri', 'dataset_description', 'dataset_file_name' keys.",
    )
    license_label: Series[String] = Field(
        nullable=False,
        description="License type for data usage and distribution. Should be one of the terms from under SWO:0000002 (license).",
        isin=[
            "CC0", # SWO:1000049
            "CC BY", # SWO:1000050
            "CC BY-SA", # SWO:1000052
            "CC BY-NC", # SWO:1000079
            "CC BY-ND", # SWO:1000077
            "CC0 1.0", # SWO:1000049
            "CC BY 2.0", # SWO:1000050
            "CC BY-SA 2.0", # SWO:1000052
            "CC BY 4.0", # SWO:1000065
            "CC BY 2.0 UK", # SWO:1000067
            "CC BY 2.1 JP", # SWO:1000072
            "CC BY 2.5", # SWO:1000073
            "CC BY 3.0 AU", # SWO:1000074
            "CC BY 3.0", # SWO:1000075
            "CC BY 3.0 US", # SWO:1000076
            "CC BY-ND 3.0", # SWO:1000077
            "CC BY-ND 4.0", # SWO:1000078
            "CC BY-NC 3.0", # SWO:1000079
            "CC BY-NC 4.0", # SWO:1000080
            "CC BY-NC-ND 3.0", # SWO:1000081
            "CC BY-NC-ND 2.5", # SWO:1000083
            "CC BY-NC-ND 2.5 CH", # SWO:1000084
            "CC BY-NC-ND 4.0", # SWO:1000085
            "CC BY-NC-SA 2.5", # SWO:1000086
            "CC BY-NC-SA 3.0", # SWO:1000087
            "CC BY-NC-SA 3.0 US", # SWO:1000088
            "CC BY-NC-SA 2.5 IN", # SWO:1000089
            "CC BY-NC-SA 4.0", # SWO:1000090
            "CC BY-SA 2.1 JP", # SWO:1000091
            "CC BY-SA 3.0", # SWO:1000092
            "CC BY-SA 3.0 US", # SWO:1000093
            "CC BY-SA 4.0", # SWO:1000094
            "Other"
        ],
    )
    license_id: Series[String] = Field(
        nullable=True,
        description="License ontology term ID for data usage and distribution. Should be one of the terms from under SWO:0000002 (license).",
        isin=[
            "SWO:1000049", # CC0
            "SWO:1000050", # CC BY
            "SWO:1000052", # CC BY-SA
            "SWO:1000079", # CC BY-NC
            "SWO:1000077", # CC BY-ND
            "SWO:1000065", # CC BY 4.0
            "SWO:1000067", # CC BY 2.0 UK
            "SWO:1000072", # CC BY 2.1 JP
            "SWO:1000073", # CC BY 2.5
            "SWO:1000074", # CC BY 3.0 AU
            "SWO:1000075", # CC BY 3.0
            "SWO:1000076", # CC BY 3.0 US
            "SWO:1000078", # CC BY-ND 4.0
            "SWO:1000080", # CC BY-NC 4.0
            "SWO:1000081", # CC BY-NC-ND 3.0
            "SWO:1000083", # CC BY-NC-ND 2.5
            "SWO:1000084", # CC BY-NC-ND 2.5 CH
            "SWO:1000085", # CC BY-NC-ND 4.0
            "SWO:1000086", # CC BY-NC-SA 2.5
            "SWO:1000087", # CC BY-NC-SA 3.0
            "SWO:1000088", # CC BY-NC-SA 3.0 US
            "SWO:1000089", # CC BY-NC-SA 2.5 IN
            "SWO:1000090", # CC BY-NC-SA 4.0
            "SWO:1000091", # CC BY-SA 2.1 JP
            "SWO:1000092", # CC BY-SA 3.0
            "SWO:1000093", # CC BY-SA 3.0 US
            "SWO:1000094", # CC BY-SA 4.0
            
        ]
    )
    curation_agent_type: Series[String] = Field(
        nullable=False,
        description="Type of agent that curated this dataset: 'human' for manual curation, 'LLM' for automated curation by a language model.",
        isin=["human", "LLM"],
    )
    curation_agent_name: Series[String] = Field(
        nullable=False,
        description="Name or identifier of the curator. For humans: full name (e.g., 'John Doe'). For LLMs: model identifier (e.g., 'google/gemini-3.5-flash').",
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
        description="Unique identifier for each gene. Usually the Ensembl gene ID, or whatever unique IDs the dataset came with",
    )
    ensembl_gene_id: Series[str] = Field(
        nullable=True,
        str_matches=r"^(ENSG|control)",  # starts with either ENSG or control
        description="Ensembl gene ID",
    )
    gene_symbol: Series[str] = Field(
        nullable=True, coerce=True, description="Gene symbol"
    )

    class Config:
        strict = "filter"
        coerce = True
        ordered = True
