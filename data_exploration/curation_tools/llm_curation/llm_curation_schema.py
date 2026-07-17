"""Pydantic schemas used for metadata curation."""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class CurationSchema(BaseModel):

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    dataset_id: str = Field(
        default=...,
        description="Unique identifier for the dataset, follows the format <firstauthor_year>",
    )

    data_modality: Literal["Perturb-seq", "CRISPR screen", "MAVE"] = Field(
        default=...,
        description="Data modality of the dataset.",
    )

    perturbation_type_label: Literal["CRISPRn", "CRISPRi", "CRISPRa", "DMS"] = Field(
        default=...,
        description="Perturbation type ontology term label of the investigated sample.",
    )

    timepoint_post_transfection: str | None = Field(
        default=None,
        description="Timepoint of the investigated sample in ISO 8601 format, starting from the time of library transfection. Example: P1DT12H30M15S",
    )

    differentiation_timepoint: str | None = Field(
        default=None,
        description="Differentiation timepoint of the investigated sample in ISO 8601 format, starting from the moment the induction of differentiation began. Example: P1DT12H30M15S",
    )

    experimental_timepoint: str | None = Field(
        default=None,
        description="Experimental timepoint of the investigated sample in ISO 8601 format, starting from the moment the main experiment began. Example: P1DT12H30M15S",
    )

    treatment_label: str | None = Field(
        default=...,
        description="Treatment/compound ontology term label used to stimulate the investigated sample. ChEMBL compound label for chemical entities. Use 'untreated control' for untreated samples where other samples were treated.",
    )

    treatment_dose: float | None = Field(
        default=None,
        description="Treatment/compound dose used to stimulate the investigated sample.",
    )

    treatment_unit: (
        Literal[
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
        ]
        | None
    ) = Field(
        default=None,
        description="Treatment/compound unit used to stimulate the investigated sample. Use 'u' for micro (e.g., 'uM' instead of 'μM').",
    )

    model_system_label: Literal[
        "cell_line", "primary_cell", "organoid", "yeast", "bacteria", "Other"
    ] = Field(
        default=...,
        description="Model system ontology term label of the investigated sample.",
    )

    species: Literal["Homo sapiens"] = Field(
        default=...,
        description="Species name of the investigated sample.",
    )

    tissue_label: str | None = Field(
        default=...,
        description="Tissue ontology term label of the investigated sample. Must be part of the UBERON ontology.",
    )

    cell_type_label: str | None = Field(
        default=...,
        description="Cell type ontology term label of the investigated sample. Must be part of the Cell Ontology (CL).",
    )

    cell_line_label: str | None = Field(
        default=...,
        description="Cell line ontology term label of the investigated sample. Must be part of the Cell Line Ontology (CLO).",
    )

    sex_label: Literal["female", "male", "mixed", "unknown"] | None = Field(
        default=...,
        description="Sex ontology term label of the investigated sample.",
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
        description="Developmental stage ontology term label of the investigated sample. Age brackets: embryonic - upto 8th week of gestation; fetal - 8th week - 40 weeks of gestation; child (0-12 years old); adolescent (13-18 years old); adult (19-59 years old); senior adult (60+ years old).",
    )

    disease_label: str | None = Field(
        default=...,
        description="Disease ontology term label of the investigated sample. Must be part of the MONDO ontology.",
    )

    study_title: str = Field(
        default=...,
        description="Title of the study/publication.",
    )

    study_uri: str = Field(
        default=...,
        description="URI/DOI of the study/publication.",
    )

    study_year: int = Field(
        default=...,
        description="Publication year of the study/publication.",
    )

    first_author: str | None = Field(
        default=...,
        description="Full name of the first author of the study/publication.",
    )

    last_author: str | None = Field(
        default=...,
        description="Full name of the last author of the study/publication.",
    )

    experiment_title: str = Field(
        default=...,
        description="Title of the experiment.",
    )

    experiment_summary: str | None = Field(
        default=...,
        description="Summary of the experiment.",
    )

    library_generation_type_label: (
        Literal[
            "endogenous genetic perturbation method",
            "exogenous genetic perturbation method",
        ]
        | None
    ) = Field(
        default=...,
        description="Library generation type ontology term label, defined in EFO under parent term EFO:0022867 (genetic perturbation)",
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

    library_name: str | None = Field(
        default=...,
        description="Name of the perturbation library. Example: Bassik Human CRISPR Knockout Library",
    )

    library_uri: str | None = Field(
        default=...,
        description="URI/accession of the perturbation library.",
    )

    library_format_label: (
        Literal["pooled", "arrayed", "arrayed|pooled", "in vivo"] | None
    ) = Field(
        default=...,
        description="Perturbation library format ontology term label.",
    )

    library_scope_label: Literal["focused", "genome-wide"] | None = Field(
        default=...,
        description="Perturbation library scope ontology term label.",
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

    library_manufacturer: str | None = Field(
        default=...,
        description="Name of the library manufacturer/vendor/origin lab. Example: Bassik",
    )

    library_lentiviral_generation: str | None = Field(
        default=...,
        description="Generation number of the lentiviral library. Example: 3",
    )

    library_grnas_per_target: str | None = Field(
        default=...,
        description="Number of gRNAs per target. Example: 4, 5-7",
    )

    library_total_grnas: str | None = Field(
        default=...,
        description="Total number of gRNAs in the library. Example: 20,000",
    )

    library_total_variants: int | None = Field(
        default=...,
        description="Only for MAVE studies; Total number of variants in the library. Example: 5,000",
        ge=0,
    )

    readout_dimensionality_label: (
        Literal["single-dimensional assay", "high-dimensional assay"] | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the dimensionality of the readout assay.",
    )

    readout_type_label: (
        Literal["transcriptomic", "proteomic", "phenotypic", "Other"] | None
    ) = Field(
        default=...,
        description="Ontology term label associated with the type of the readout assay.",
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

    readout_measurement_label: (
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

    method_uri: str | None = Field(
        default=...,
        description="URI associated with the method used in the readout assay.",
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

    software_counts_label: (
        Literal["custom", "MaGeCK", "CellRanger", "Drop-seq Tools", "Other"] | None
    ) = Field(
        default=...,
        description="Ontology term label for the software used for generating counts.",
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

    reference_genome_label: Literal["GRCh38", "GRCh37", "Other"] | None = Field(
        default=...,
        description="Ontology term label for the reference genome.",
    )

    significance_criteria: str | None = Field(
        default=...,
        description="Criteria used to determine significance, e.g., FDR < 0.05.",
    )

    score_interpretation: str | None = Field(
        default=...,
        description="Interpretation of the perturbation effect score, e.g. negative values = depletion; positive values = enrichment, or negative values = decreased phosphatase activity; positive values = increased phosphatase activity",
    )

    associated_datasets: str | None = Field(
        default=...,
        description="List of associated datasets with each dataset having 'dataset_accession', 'dataset_uri', 'dataset_description', 'dataset_file_name' keys.",
    )

    license_label: Literal[
        "CC0", "CC BY", "CC BY-SA", "CC BY-NC", "CC BY-ND", "Other"
    ] = Field(
        default=...,
        description="License type for data usage and distribution. Should be one of the terms from under SWO:0000002 (license).",
    )

    curation_agent_type: Literal["human", "LLM"] = Field(
        default=...,
        description="Type of agent that curated this dataset: 'human' for manual curation, 'LLM' for automated curation by a language model.",
    )
    curation_agent_name: str = Field(
        default=...,
        description="Name or identifier of the curator. For humans: full name (e.g., 'John Doe'). For LLMs: model identifier (e.g., 'google/gemini-3.5-flash').",
    )


class EvidenceExtractionSchema(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    dataset_id_evidence: str | None = Field(
        default=None,
        description="Unique identifier for the dataset, follows the format <firstauthor_year>. Example: The complete sequence datasets were deposited under the accession number GEO: GSE188426.",
    )
    data_modality_evidence: str | None = Field(
        default=None,
        description="Data modality of the dataset. Example: In this work, we employed Perturb-seq, a multiplexed assay of variant effect, to query all single amino acid changes.",
    )
    perturbation_type_label_evidence: str | None = Field(
        default=None,
        description="Perturbation type of the investigated sample. Example: We utilized CRISPRi to silence target promoters across the cell pool.",
    )
    timepoint_post_transfection_evidence: str | None = Field(
        default=None,
        description="Timepoint of the investigated sample in ISO 8601 format, starting from the time of library transfection into the model system. Example: Cells were harvested 72 hours post-transfection for downstream analysis.",
    )
    differentiation_timepoint_evidence: str | None = Field(
        default=None,
        description="Differentiation timepoint of the investigated sample in ISO 8601 format, starting from the moment the induction of differentiation began. Example: On day 5 of differentiation, cells were collected to evaluate lineage markers.",
    )
    experimental_timepoint_evidence: str | None = Field(
        default=None,
        description="Experimental timepoint of the investigated sample. Example: Cells were collected at 24 hours post-treatment to assess early transcriptional responses. or Example: Samples were harvested at 48 hours post-infection to evaluate viral replication dynamics.",
    )
    treatment_label_evidence: str | None = Field(
        default=None,
        description="Treatment/compound used to stimulate the investigated sample. Use 'untreated control' for untreated samples where other samples were treated. Example: The engineered cell lines were treated with cisplatin for a duration of 48 hours.",
    )
    treatment_dose_evidence: str | None = Field(
        default=None,
        description="Treatment/compound dose (numeric value) used to stimulate the investigated sample. Example: Cells were stimulated with 10 uM of TNFalpha for 24 hours. or Example: The engineered cell lines were treated with 25 ug/mL of Amyloid-beta peptide for 12 hours.",
    )
    treatment_unit_evidence: str | None = Field(
        default=None,
        description="Treatment/compound unit used to stimulate the investigated sample. Example: Cells were stimulated with 10 uM of TNFalpha for 24 hours. or Example: The engineered cell lines were treated with 25 ug/mL of Amyloid-beta peptide for 12 hours.",
    )
    model_system_label_evidence: str | None = Field(
        default=None,
        description="Experimental model system of the investigated sample. Example: All parallel screenings were carried out in an immortalized human cell line. or Example: A humanised yeast strain was used to assess the functional impact of BRCA1 variants in a high-throughput manner.",
    )
    species_evidence: str | None = Field(
        default=None,
        description="Species name of the investigated sample. Example: Cells were derived from Homo sapiens donors according to approved protocols.",
    )
    tissue_label_evidence: str | None = Field(
        default=None,
        description="Tissue of the investigated sample. If primary cells are used, specify the tissue of origin. If cell lines are used, specify the tissue from which the cell line was originally derived. Example: Primary cells were isolated from human peripheral blood mononuclear cells.",
    )
    cell_type_label_evidence: str | None = Field(
        default=None,
        description="Cell type of the investigated sample. Example: We obtained human CD4+ T cells using negative magnetic selection.",
    )
    cell_line_label_evidence: str | None = Field(
        default=None,
        description="Cell line of the investigated sample. Only relevant for experiments where model system is a cell line. Example: All expression constructs were stably transfected into HEK293T cells.",
    )
    sex_label_evidence: str | None = Field(
        default=None,
        description="Sex of the investigated sample. Example: These primary cells were derived from a 24-year-old female donor.",
    )
    developmental_stage_label_evidence: str | None = Field(
        default=None,
        description="Developmental stage of the investigated sample. Example: The study analyzed adult human hepatocytes sourced from tissue biopsies.",
    )
    disease_label_evidence: str | None = Field(
        default=None,
        description="Disease associated with the investigated sample. Example: We profiled patient-derived models of chronic myeloid leukemia.",
    )
    study_title_evidence: str | None = Field(
        default=None,
        description="Title of the study/publication. Extract as is. Example: The functional impact of BRCA1 BRCT domain variants was analyzed using multiplexed assays.",
    )
    study_uri_evidence: str | None = Field(
        default=None,
        description="URI/DOI of the study/publication. Example: The peer-reviewed paper is available online via https://doi.org/10.1016/j.ajhg.2022.01.019.",
    )
    study_year_evidence: str | None = Field(
        default=None,
        description="Publication year of the study/publication. Example: First published online in January 2022.",
    )
    first_author_evidence: str | None = Field(
        default=None,
        description="Full name of the first author of the study/publication. Example: Authors: John A Smith, Peter B Johnson, Roberta C Williams.",
    )
    last_author_evidence: str | None = Field(
        default=None,
        description="Full name of the last author of the study/publication. Example: Authors: Jane B Doe, Alexander C Brown, Josephine D White.",
    )
    experiment_title_evidence: str | None = Field(
        default=None,
        description="Title of the specific experiment. Example: We conducted a cisplatin resistance growth assay to measure BRCA1 variant function.",
    )
    experiment_summary_evidence: str | None = Field(
        default=None,
        description="Summary of the experiment, briefly describing the key elements of the experimental setup. Example: The multiplexed reporter assay separates functionally normal and abnormal variants by selecting for cisplatin resistance.",
    )
    library_generation_type_label_evidence: str | None = Field(
        default=None,
        description="Broad library generation type. Must represent whether the library was generated endogenously (manipulation of the host organism's genome) or exogenously (introduction of genetic material from an external source). Example: Mutagenesis library was generated using error-prone PCR method. or Example: The CRISPRi library was generated using a lentiviral vector system to introduce guide RNAs into the host cells.",
    )
    library_generation_method_label_evidence: str | None = Field(
        default=None,
        description="Specific library generation method. Example: The libraries were constructed using site-directed mutagenesis libraries. or Example: The libraries were constructed using a microarray synthesis method to generate a diverse set of oligonucleotides.",
    )
    enzyme_delivery_method_label_evidence: str | None = Field(
        default=None,
        description="Enzyme delivery method. Example: Transfections were performed with 30 μL Lipofectamine and 1 mL of Opti-MEM per plate using lipofection.",
    )
    library_delivery_method_label_evidence: str | None = Field(
        default=None,
        description="Library delivery method. Example: The mutant plasmid library was introduced into the cells via lentivirus transduction.",
    )
    enzyme_integration_state_label_evidence: str | None = Field(
        default=None,
        description="Enzyme integration state into the host genome. Example: Cas9 expression was maintained from a stably integrated transgene with random locus integration in the host genome.",
    )
    library_integration_state_label_evidence: str | None = Field(
        default=None,
        description="Library integration state into the host genome. Example: The stably integrated variant library was selected with random locus integration using hygromycin B.",
    )
    enzyme_expression_control_label_evidence: str | None = Field(
        default=None,
        description="Enzyme expression control mechanism. Example: The cell line features TET-ON inducible transgene expression system of the Cas9 enzyme.",
    )
    library_expression_control_label_evidence: str | None = Field(
        default=None,
        description="Library expression control mechanism. Example: Expression of the guide library was driven by the constitutive transgene expression human U6 promoter.",
    )
    library_name_evidence: str | None = Field(
        default=None,
        description="Name of the perturbation library. Example: We synthesized the human CRISPRi v2 library targeting transcription factors. or Example: Dolcetto library was used to perform genome-wide CRISPRi screens in K562 cells.",
    )
    library_uri_evidence: str | None = Field(
        default=None,
        description="URI/accession of the perturbation library. Only relevant for libraries that have been made public through repositories such as Addgene. Example: The physical plasmids were obtained from Addgene under catalog number #1000000019.",
    )
    library_format_label_evidence: str | None = Field(
        default=None,
        description="Perturbation library format. Example: We performed a pooled survival assay to screen all variants in parallel. or Example: The arrayed library was screened in a 96-well plate format, with each well containing a unique variant.",
    )
    library_scope_label_evidence: str | None = Field(
        default=None,
        description="Perturbation library scope. Example: The focused library targets 100 essential kinases in the human genome. or Example: The genome-wide library targeting all protein-coding genes was used to perform a comprehensive CRISPR screen in THP-1 cells.",
    )
    library_perturbation_type_label_evidence: str | None = Field(
        default=None,
        description="Perturbation type of the library, specifying activation/inhibition/knockout/base editing/prime editing etc. Example: A site-saturation mutagenesis library was constructed to introduce single amino acid substitutions. or Example: The ABE8 base editor was used to generate a library of point mutations in the second exon of BRCA1 gene, resulting in a library of single-nucleotide variants.",
    )
    library_manufacturer_evidence: str | None = Field(
        default=None,
        description="Name of the library manufacturer/origin lab. Example: The library was a generous gift from the Bassik lab.",
    )
    library_lentiviral_generation_evidence: str | None = Field(
        default=None,
        description="Generation number of the lentiviral library. Example: Lentivirus was generated using a second-generation packaging system containing psPAX2 and pMD2.G.",
    )
    library_grnas_per_target_evidence: str | None = Field(
        default=None,
        description="Number of gRNAs per target. Example: The library contains an average of 10 sgRNAs per target gene.",
    )
    library_total_grnas_evidence: str | None = Field(
        default=None,
        description="Total number of gRNAs in the library. Example: A total of 180,000 unique sgRNAs were included in the synthesized pool.",
    )
    library_total_variants_evidence: str | None = Field(
        default=None,
        description="Only for MAVE studies; Total number of variants in the library. Extract from the MAVE DB metadata or from the evidence. Example: A total of 1,427 single-residue variants were assessed in the reporter assay.",
    )
    readout_dimensionality_label_evidence: str | None = Field(
        default=None,
        description="Dimensionality of the readout assay. Example: We performed single-cell transcriptomic profiling on the harvested population. or Example: We conducted a cell viability assay after drug treatment to measure the survival rate of cells.",
    )
    readout_type_label_evidence: str | None = Field(
        default=None,
        description="Type of the readout assay. Example: The survival rate of cells was measured using a cell viability phenotypic assay. or Example: We performed single-cell RNA sequencing to profile the transcriptomic changes in response to perturbations.",
    )
    readout_technology_label_evidence: str | None = Field(
        default=None,
        description="A specific technology type used in the readout assay. Example: Single cells were partitioned using the 10x Genomics Chromium controller with single-cell RNA-seq. or Example: The cell viability assay was performed using a luminescence-based readout on a plate reader.",
    )
    readout_measurement_label_evidence: str | None = Field(
        default=None,
        description="Measurement type associated with the readout assay. Example: We measured cellular abundance as a proxy for variant fitness after drug treatment targeting cell viability. or Example: We quantified surface protein expression using flow cytometry to assess the impact of perturbations on cell signaling pathways.",
    )
    method_name_label_evidence: str | None = Field(
        default=None,
        description="Method name used in the readout assay. Example: We mapped genetic interactions using Perturb-seq with direct guide capture. or Example: DMS-TileSeq was employed to assess the functional impact of BRCA1 variants by measuring their effects on protein stability and activity.",
    )
    method_uri_evidence: str | None = Field(
        default=None,
        description="URI associated with the method used in the readout assay. Example: The protocol details are registered on Zenodo under accession 10.5281/zenodo.7574261.",
    )
    sequencing_library_kit_label_evidence: str | None = Field(
        default=None,
        description="Sequencing library kit used. Example: Sequencing libraries were prepared using the 10x Genomics Chromium GEM-X Single Cell 5-prime kit v3.",
    )
    sequencing_platform_label_evidence: str | None = Field(
        default=None,
        description="Sequencing platform used. Example: The library was sequenced on an Illumina NovaSeq 6000 instrument.",
    )
    sequencing_strategy_label_evidence: str | None = Field(
        default=None,
        description="Sequencing strategy used. Example: Both the guide barcodes and cellular indexes were sequenced using barcode sequencing on the Illumina HiSeq 4000.",
    )
    software_counts_label_evidence: str | None = Field(
        default=None,
        description="Software used for generating counts. Example: Raw sequencing reads were processed using the CellRanger software package to generate count matrices.",
    )
    software_analysis_label_evidence: str | None = Field(
        default=None,
        description="Software used for analysis. Example: Differential expression analysis was performed using the Seurat packages.",
    )
    reference_genome_label_evidence: str | None = Field(
        default=None,
        description="Reference genome used. Example: All transcripts were aligned to the GRCh38 human reference genome.",
    )
    significance_criteria_evidence: str | None = Field(
        default=None,
        description="Criteria used to determine significance, e.g., FDR < 0.05. Example: A false discovery rate threshold of FDR < 0.05 was used to define significant hits.",
    )
    score_interpretation_evidence: str | None = Field(
        default=None,
        description="Interpretation of the perturbation effect score, e.g. negative values = depletion; positive values = enrichment, or negative values = decreased phosphatase activity; positive values = increased phosphatase activity. Example: Functional scores were calculated as log2 fold-change ratios, where negative values indicate depletion.",
    )
    associated_datasets_evidence: str | None = Field(
        default=None,
        description="List of associated datasets deposited in public repositories. Example: Raw sequencing reads have been deposited in the Gene Expression Omnibus under GSE188426.",
    )
    license_label_evidence: str | None = Field(
        default=None,
        description="License type for data usage and distribution. Example: All datasets are distributed under the Creative Commons Attribution CC BY International License.",
    )


class SpecificTermExtractionSchema(BaseModel):

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    dataset_id: str | None = Field(
        default=None,
        description="Unique identifier for the dataset, follows the format <firstauthor_year>",
    )

    data_modality: Literal["Perturb-seq", "CRISPR screen", "MAVE"] | None = Field(
        default=None,
        description="Data modality of the dataset.",
    )

    perturbation_type_label: Literal["CRISPRn", "CRISPRi", "CRISPRa", "DMS"] | None = (
        Field(
            default=None,
            description="Perturbation type ontology term label of the investigated sample.",
        )
    )

    timepoint_post_transfection: str | None = Field(
        default=None,
        description="Timepoint of the investigated sample in ISO 8601 format, starting from the time of library transfection. Example: P1DT12H30M15S",
    )

    differentiation_timepoint: str | None = Field(
        default=None,
        description="Differentiation timepoint of the investigated sample in ISO 8601 format, starting from the moment the induction of differentiation began. Example: P1DT12H30M15S",
    )

    experimental_timepoint: str | None = Field(
        default=None,
        description="Experimental timepoint of the investigated sample in ISO 8601 format. Example: P1DT12H30M15S",
    )

    treatment_label: str | None = Field(
        default=None,
        description="Treatment/compound ontology term label used to stimulate the investigated sample. ChEMBL compound label for chemical entities. Use 'untreated control' for untreated samples where other samples were treated.",
    )

    treatment_dose: float | None = Field(
        default=None,
        description="Treatment/compound dose used to stimulate the investigated sample.",
    )

    treatment_unit: (
        Literal[
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
        ]
        | None
    ) = Field(
        default=None,
        description="Treatment/compound unit used to stimulate the investigated sample. Use 'u' for micro (e.g., 'uM' instead of 'μM').",
    )

    model_system_label: (
        Literal["cell_line", "primary_cell", "organoid", "yeast", "bacteria", "Other"]
        | None
    ) = Field(
        default=None,
        description="Model system ontology term label of the investigated sample.",
    )

    species: Literal["Homo sapiens"] | None = Field(
        default=None,
        description="Species name of the investigated sample.",
    )

    tissue_label: str | None = Field(
        default=None,
        description="Tissue ontology term label of the investigated sample.",
    )

    cell_type_label: str | None = Field(
        default=None,
        description="Cell type ontology term label of the investigated sample.",
    )

    cell_line_label: str | None = Field(
        default=None,
        description="Cell line ontology term label of the investigated sample.",
    )

    sex_label: Literal["female", "male", "mixed", "unknown"] | None = Field(
        default=None,
        description="Sex ontology term label of the investigated sample.",
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
        default=None,
        description="Developmental stage ontology term label of the investigated sample. Age brackets: embryonic - upto 8th week of gestation; fetal - 8th week - 40 weeks of gestation; child (0-12 years old); adolescent (13-18 years old); adult (19-59 years old); senior adult (60+ years old).",
    )

    disease_label: str | None = Field(
        default=None,
        description="Disease ontology term label of the investigated sample.",
    )

    study_title: str | None = Field(
        default=None,
        description="Title of the study/publication.",
    )

    study_uri: str | None = Field(
        default=None,
        description="URI/DOI of the study/publication.",
    )

    study_year: int | None = Field(
        default=None,
        description="Publication year of the study/publication.",
    )

    first_author: str | None = Field(
        default=None,
        description="Full name of the first author of the study/publication.",
    )

    last_author: str | None = Field(
        default=None,
        description="Full name of the last author of the study/publication.",
    )

    experiment_title: str | None = Field(
        default=None,
        description="Title of the specific experiment. Extract from the MAVE DB metadata or from the evidence.",
    )

    experiment_summary: str | None = Field(
        default=None,
        description="Summary of the specific experiment. Extract from the MAVE DB metadata or from the evidence.",
    )

    library_generation_type_label: (
        Literal[
            "endogenous genetic perturbation method",
            "exogenous genetic perturbation method",
        ]
        | None
    ) = Field(
        default=None,
        description="Library generation type ontology term label. Endogenous genetic perturbation method - A genetic perturbation method that involves the manipulation of the host organism's genome. Exogenous genetic perturbation method - A genetic perturbation method that involves the introduction of foreign genetic material into a host organism, such as a library of synthetic sequences encoding variants of interest for a particular gene.",
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
        default=None,
        description="Library generation method ontology term label.",
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
        default=None,
        description="Enzyme delivery method ontology term label.",
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
        default=None,
        description="Library delivery method ontology term label.",
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
        default=None,
        description="Enzyme integration state ontology term label.",
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
        default=None,
        description="Library integration state ontology term label.",
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
        default=None,
        description="Enzyme expression control ontology term label.",
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
        default=None,
        description="Library expression control ontology term label.",
    )

    library_name: str | None = Field(
        default=None,
        description="Name of the perturbation library. Example: Bassik Human CRISPR Knockout Library",
    )

    library_uri: str | None = Field(
        default=None,
        description="URI/accession of the perturbation library.",
    )

    library_format_label: (
        Literal["pooled", "arrayed", "arrayed|pooled", "in vivo"] | None
    ) = Field(
        default=None,
        description="Perturbation library format ontology term label.",
    )

    library_scope_label: Literal["focused", "genome-wide"] | None = Field(
        default=None,
        description="Perturbation library scope ontology term label.",
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
        default=None,
        description="Ontology term label for the library perturbation type.",
    )

    library_manufacturer: str | None = Field(
        default=None,
        description="Name of the library vendor/origin lab. Example: Bassik",
    )

    library_lentiviral_generation: str | None = Field(
        default=None,
        description="Generation number of the lentiviral library. Example: 3",
    )

    library_grnas_per_target: str | None = Field(
        default=None,
        description="Number of gRNAs per target. Example: 4, 5-7",
    )

    library_total_grnas: int | None = Field(
        default=None,
        description="Total number of gRNAs in the library. Example: 20,000",
    )

    library_total_variants: int | None = Field(
        default=None,
        description="Only for MAVE studies; Total number of variants in the library. Example: 5,000",
        ge=0,
    )

    readout_dimensionality_label: (
        Literal["single-dimensional assay", "high-dimensional assay"] | None
    ) = Field(
        default=None,
        description="Ontology term label associated with the dimensionality of the readout assay.",
    )

    readout_type_label: (
        Literal["transcriptomic", "proteomic", "phenotypic", "Other"] | None
    ) = Field(
        default=None,
        description="Ontology term label associated with the type of the readout assay.",
    )

    readout_technology_label: (
        Literal[
            "single-cell rna-seq", "population growth assay", "flow cytometry", "Other"
        ]
        | None
    ) = Field(
        default=None,
        description="Ontology term label associated with the technology used in the readout assay.",
    )

    readout_measurement_label: (
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
        default=None,
        description="Ontology term label associated with the measurement type of the readout assay.",
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
        default=None,
        description="Ontology term label associated with the method name used in the readout assay.",
    )

    method_uri: str | None = Field(
        default=None,
        description="URI associated with the method used in the readout assay.",
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
        default=None,
        description="Ontology term label associated with the sequencing library kit.",
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
        default=None,
        description="Ontology term label associated with the sequencing platform.",
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
        default=None,
        description="Ontology term label associated with the sequencing strategy.",
    )

    software_counts_label: (
        Literal["custom", "MaGeCK", "CellRanger", "Drop-seq Tools", "Enrich2", "Other"] | None
    ) = Field(
        default=None,
        description="Ontology term label for the software used for generating counts.",
    )

    software_analysis_label: (
        Literal[
            "custom", "MAGeCK", "Achilles", "TRADE", "Seurat", "MAST", "scanpy", "Enrich2", "Other"
        ]
        | None
    ) = Field(
        default=None,
        description="Ontology term label for the software used for analysis.",
    )

    reference_genome_label: Literal["GRCh38", "GRCh37", "Other"] | None = Field(
        default=None,
        description="Ontology term label for the reference genome.",
    )

    significance_criteria: str | None = Field(
        default=None,
        description="Criteria used to determine significance, e.g., FDR < 0.05.",
    )

    score_interpretation: str | None = Field(
        default=None,
        description="Interpretation of the perturbation effect score, e.g. negative values = depletion; positive values = enrichment, or negative values = decreased phosphatase activity; positive values = increased phosphatase activity",
    )

    associated_datasets: str | None = Field(
        default=None,
        description="List of associated datasets with each dataset having 'dataset_accession', 'dataset_uri', 'dataset_description', 'dataset_file_name' keys.",
    )

    license_label: (
        Literal["CC0", "CC BY", "CC BY-SA", "CC BY-NC", "CC BY-ND", "Other"] | None
    ) = Field(
        default=None,
        description="License type for data usage and distribution. Should be one of the terms from under SWO:0000002 (license).",
    )

