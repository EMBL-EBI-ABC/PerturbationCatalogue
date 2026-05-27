from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


CONTROLLED_VOCAB_FIELD_MAPPING: dict[str, str] = {
	"phenotypic_assay_profiling_strategy": "Phenotypic Assay Profiling Strategy",
	"phenotypic_assay_sequencing_read_type": "Phenotypic Assay Sequencing Read Type",
	"variant_library_creation_method": "Variant Library Creation Method",
	"delivery_method": "Delivery Method",
	"endogenous_locus_library_method_mechanism": "Endogenous Locus Library Method Mechanism",
	"endogenous_locus_library_method_system": "Endogenous Locus Library Method System",
	"phenotypic_assay_dimensionality": "Phenotypic Assay Dimensionality",
	"phenotypic_assay_mechanism": "Phenotypic Assay Mechanism",
	"phenotypic_assay_model_system": "Phenotypic Assay Model System",
	"phenotypic_assay_method": "Phenotypic Assay Method",
	"molecular_mechanism_assessed": "Molecular Mechanism Assessed",
	"in_vitro_construct_library_method_mechanism": "In Vitro Construct Library Method Mechanism",
	"in_vitro_construct_library_method_system": "In Vitro Construct Library Method System",
}


PhenotypicAssayProfilingStrategy = Literal[
	"Barcode sequencing",
	"Direct sequencing",
	"Other",
	"Shotgun sequencing",
]

PhenotypicAssaySequencingReadType = Literal[
	"Multi-segment",
	"Other",
	"Single-segment (long read)",
	"Single-segment (short read)",
]

VariantLibraryCreationMethod = Literal[
	"Endogenous locus library method",
	"In vitro construct library method",
	"N/A (meta-analysis)",
	"Other",
]

DeliveryMethod = Literal[
	"Adeno-associated virus transduction",
	"Chemical-based transfection",
	"Chemical or heat shock transformation",
	"Electroporation",
	"Lentivirus transduction",
	"Nucleofection",
	"Other",
	"Retroviral transduction",
]

EndogenousLocusLibraryMethodMechanism = Literal[
	"Base editor",
	"Nuclease",
	"Other",
	"Prime editor",
]

EndogenousLocusLibraryMethodSystem = Literal[
	"AsCas12a",
	"Other",
	"RfsCas13d",
	"SaCas9",
	"SpCas9",
]

PhenotypicAssayDimensionality = Literal[
	"Combined functional data",
	"High-dimensional data",
	"Other",
	"Single-dimensional data",
]

PhenotypicAssayMechanism = Literal[
	"Dominant-negative effect",
	"Gain of function",
	"Loss of function",
	"Mixed functional effect",
	"Other",
]

PhenotypicAssayModelSystem = Literal[
	"Bacteria",
	"Bacteriophage",
	"Immortalized human cells",
	"Immortalized non-human mammalian cells",
	"Induced pluripotent stem cells from human female",
	"Induced pluripotent stem cells from human male",
	"Molecular display",
	"Murine primary cells",
	"Not applicable",
	"Organoids",
	"Other",
	"Patient derived primary cells (e.g. T-cells, adipocytes)",
	"Yeast",
]

PhenotypicAssayMethod = Literal[
	"Abundance assay",
	"Binding assay",
	"Bacterial two-hybrid assay",
	"Bulk RNA-sequencing",
	"Cell fitness",
	"Cell morphology assay",
	"Cell proliferation assay",
	"Cell proliferation assay with genetic complementation",
	"Direct protein function",
	"Electrophysiological method",
	"Flow cytometry assay",
	"Fluorescence in-situ hybridization (FISH) assay",
	"Genetic complementation",
	"Homology-directed DNA repair frequency measurement",
	"Imaging mass cytometry assay",
	"Ion channel assay",
	"Multiplexed fluorescent antibody imaging",
	"Not applicable (trained predictor)",
	"One-hybrid assay",
	"Other",
	"Polysome profiling assay",
	"Posttranslation modification assay",
	"Protein localization/trafficking assay",
	"Protein stability assay",
	"Reporter gene assay",
	"Single cell imaging",
	"Single-cell RNA sequencing assay",
	"Splicing assay",
	"Survival assessment assay",
	"Systematic evolution of ligands by exponential enrichment assay",
	"Yeast two-hybrid assay",
]

MolecularMechanismAssessed = Literal[
	"Calcium-mediated signaling",
	"Catalytic activity and Ornithine carbamoyltransferase activity",
	"Catalytic and Cysteine synthase activity",
	"Catalytic and Gluconokinase activity",
	"Catalytic and Hydroxymethylbilane synthase activity",
	"Catalytic and Lipid phosphatase activity",
	"Catalytic and Thiamine diphosphokinase activity",
	"Cellular response to cisplatin, DNA repair and Double-strand break repair via homologous recombination",
	"Core promoter sequence-specific DNA binding",
	"DNA and Mismatch repair",
	"DNA damage checkpoint signaling",
	"DNA damage response, signal transduction by p53 class mediator and MDM2/MDM4 family protein binding",
	"DNA repair and Double-strand break repair via homologous recombination",
	"Double-strand break repair, damage response, signal transduction by p53 class mediator and Response to etoposide",
	"Molecular condensate scaffold activity",
	"Molecular function",
	"Monoatomic ion transport and Potassium channel activity",
	"Monoatomic ion transport and Sodium ion transport",
	"Oxidative phosphorylation",
	"Protein abundance",
	"Protein carboxylation",
	"Protein glycosylation",
	"Protein localization/trafficking to plasma membrane",
	"Protein-protein interaction",
	"RNA splicing",
	"Regulation of protein stability",
	"Response to misfolded protein",
	"Sodium channel activity",
	"Translation efficiency",
]

InVitroConstructLibraryMethodMechanism = Literal[
	"Episomal delivery",
	"Extra-local construct insertion",
	"Landing pad integration",
	"Native locus replacement",
	"Other",
	"Plasmid (not integrated)",
	"Random locus viral integration",
	"Transfection of RNA",
]

InVitroConstructLibraryMethodSystem = Literal[
	"Cassette mutagenesis",
	"Doped oligo synthesis",
	"Error-prone PCR",
	"Microarray synthesis",
	"Nicking mutagenesis",
	"Oligo-directed mutagenic PCR",
	"Oligo pool synthesis",
	"Other",
	"POPCode mutagenesis",
	"Proprietary method",
	"Site-directed mutagenesis",
]


def _evidence_field(original_name: str):
	return Field(
		default=None,
		alias=f"{original_name} Evidence",
		description=(
			f"Short verbatim evidence quote from the source text supporting '{original_name}'. "
			"Use null if the source text does not explicitly support this field."
		),
	)


class MavedbMetadataSchema(BaseModel):
	model_config = ConfigDict(populate_by_name=True)

	phenotypic_assay_profiling_strategy: PhenotypicAssayProfilingStrategy | None = Field(
		default=None,
		alias="Phenotypic Assay Profiling Strategy",
		description="Phenotypic Assay Profiling Strategy. Allowed values: Barcode sequencing, Direct sequencing, Other, Shotgun sequencing.",
	)
	phenotypic_assay_sequencing_read_type: PhenotypicAssaySequencingReadType | None = Field(
		default=None,
		alias="Phenotypic Assay Sequencing Read Type",
		description="Phenotypic Assay Sequencing Read Type. Allowed values: Multi-segment, Other, Single-segment (long read), Single-segment (short read).",
	)
	variant_library_creation_method: VariantLibraryCreationMethod | None = Field(
		default=None,
		alias="Variant Library Creation Method",
		description="Variant Library Creation Method. Allowed values: Endogenous locus library method, In vitro construct library method, N/A (meta-analysis), Other.",
	)
	delivery_method: DeliveryMethod | None = Field(
		default=None,
		alias="Delivery Method",
		description="Delivery Method. Allowed values: Adeno-associated virus transduction, Chemical-based transfection, Chemical or heat shock transformation, Electroporation, Lentivirus transduction, Nucleofection, Other, Retroviral transduction.",
	)
	endogenous_locus_library_method_mechanism: EndogenousLocusLibraryMethodMechanism | None = Field(
		default=None,
		alias="Endogenous Locus Library Method Mechanism",
		description="Endogenous Locus Library Method Mechanism. Allowed values: Base editor, Nuclease, Other, Prime editor.",
	)
	endogenous_locus_library_method_system: EndogenousLocusLibraryMethodSystem | None = Field(
		default=None,
		alias="Endogenous Locus Library Method System",
		description="Endogenous Locus Library Method System. Allowed values: AsCas12a, Other, RfsCas13d, SaCas9, SpCas9.",
	)
	phenotypic_assay_dimensionality: PhenotypicAssayDimensionality | None = Field(
		default=None,
		alias="Phenotypic Assay Dimensionality",
		description="Phenotypic Assay Dimensionality. Allowed values: Combined functional data, High-dimensional data, Other, Single-dimensional data.",
	)
	phenotypic_assay_mechanism: PhenotypicAssayMechanism | None = Field(
		default=None,
		alias="Phenotypic Assay Mechanism",
		description="Phenotypic Assay Mechanism. Allowed values: Dominant-negative effect, Gain of function, Loss of function, Mixed functional effect, Other.",
	)
	phenotypic_assay_model_system: PhenotypicAssayModelSystem | None = Field(
		default=None,
		alias="Phenotypic Assay Model System",
		description="Phenotypic Assay Model System. Allowed values: Bacteria, Bacteriophage, Immortalized human cells, Immortalized non-human mammalian cells, Induced pluripotent stem cells from human female, Induced pluripotent stem cells from human male, Molecular display, Murine primary cells, Not applicable, Organoids, Other, Patient derived primary cells (e.g. T-cells, adipocytes), Yeast.",
	)
	phenotypic_assay_method: PhenotypicAssayMethod | None = Field(
		default=None,
		alias="Phenotypic Assay Method",
		description="Phenotypic Assay Method. Allowed values: Abundance assay, Binding assay, Bacterial two-hybrid assay, Bulk RNA-sequencing, Cell fitness, Cell morphology assay, Cell proliferation assay, Cell proliferation assay with genetic complementation, Direct protein function, Electrophysiological method, Flow cytometry assay, Fluorescence in-situ hybridization (FISH) assay, Genetic complementation, Homology-directed DNA repair frequency measurement, Imaging mass cytometry assay, Ion channel assay, Multiplexed fluorescent antibody imaging, Not applicable (trained predictor), One-hybrid assay, Other, Polysome profiling assay, Posttranslation modification assay, Protein localization/trafficking assay, Protein stability assay, Reporter gene assay, Single cell imaging, Single-cell RNA sequencing assay, Splicing assay, Survival assessment assay, Systematic evolution of ligands by exponential enrichment assay, Yeast two-hybrid assay.",
	)
	molecular_mechanism_assessed: MolecularMechanismAssessed | None = Field(
		default=None,
		alias="Molecular Mechanism Assessed",
		description="Molecular Mechanism Assessed. Allowed values: Calcium-mediated signaling, Catalytic activity and Ornithine carbamoyltransferase activity, Catalytic and Cysteine synthase activity, Catalytic and Gluconokinase activity, Catalytic and Hydroxymethylbilane synthase activity, Catalytic and Lipid phosphatase activity, Catalytic and Thiamine diphosphokinase activity, Cellular response to cisplatin, DNA repair and Double-strand break repair via homologous recombination, Core promoter sequence-specific DNA binding, DNA and Mismatch repair, DNA damage checkpoint signaling, DNA damage response, signal transduction by p53 class mediator and MDM2/MDM4 family protein binding, DNA repair and Double-strand break repair via homologous recombination, Double-strand break repair, damage response, signal transduction by p53 class mediator and Response to etoposide, Molecular condensate scaffold activity, Molecular function, Monoatomic ion transport and Potassium channel activity, Monoatomic ion transport and Sodium ion transport, Oxidative phosphorylation, Protein abundance, Protein carboxylation, Protein glycosylation, Protein localization/trafficking to plasma membrane, Protein-protein interaction, RNA splicing, Regulation of protein stability, Response to misfolded protein, Sodium channel activity, Translation efficiency.",
	)
	in_vitro_construct_library_method_mechanism: InVitroConstructLibraryMethodMechanism | None = Field(
		default=None,
		alias="In Vitro Construct Library Method Mechanism",
		description="In Vitro Construct Library Method Mechanism. Allowed values: Episomal delivery, Extra-local construct insertion, Landing pad integration, Native locus replacement, Other, Plasmid (not integrated), Random locus viral integration, Transfection of RNA.",
	)
	in_vitro_construct_library_method_system: InVitroConstructLibraryMethodSystem | None = Field(
		default=None,
		alias="In Vitro Construct Library Method System",
		description="In Vitro Construct Library Method System. Allowed values: Cassette mutagenesis, Doped oligo synthesis, Error-prone PCR, Microarray synthesis, Nicking mutagenesis, Oligo-directed mutagenic PCR, Oligo pool synthesis, Other, POPCode mutagenesis, Proprietary method, Site-directed mutagenesis.",
	)


class MavedbMetadataExtractionSchema(BaseModel):
	model_config = ConfigDict(populate_by_name=True)

	phenotypic_assay_profiling_strategy_evidence: str | None = _evidence_field(
		"Phenotypic Assay Profiling Strategy"
	)
	phenotypic_assay_profiling_strategy: PhenotypicAssayProfilingStrategy | None = Field(
		default=None,
		alias="Phenotypic Assay Profiling Strategy",
		description="Phenotypic Assay Profiling Strategy. Allowed values: Barcode sequencing, Direct sequencing, Other, Shotgun sequencing.",
	)
	phenotypic_assay_sequencing_read_type_evidence: str | None = _evidence_field(
		"Phenotypic Assay Sequencing Read Type"
	)
	phenotypic_assay_sequencing_read_type: PhenotypicAssaySequencingReadType | None = Field(
		default=None,
		alias="Phenotypic Assay Sequencing Read Type",
		description="Phenotypic Assay Sequencing Read Type. Allowed values: Multi-segment, Other, Single-segment (long read), Single-segment (short read).",
	)
	variant_library_creation_method_evidence: str | None = _evidence_field(
		"Variant Library Creation Method"
	)
	variant_library_creation_method: VariantLibraryCreationMethod | None = Field(
		default=None,
		alias="Variant Library Creation Method",
		description="Variant Library Creation Method. Allowed values: Endogenous locus library method, In vitro construct library method, N/A (meta-analysis), Other.",
	)
	delivery_method_evidence: str | None = _evidence_field("Delivery Method")
	delivery_method: DeliveryMethod | None = Field(
		default=None,
		alias="Delivery Method",
		description="Delivery Method. Allowed values: Adeno-associated virus transduction, Chemical-based transfection, Chemical or heat shock transformation, Electroporation, Lentivirus transduction, Nucleofection, Other, Retroviral transduction.",
	)
	endogenous_locus_library_method_mechanism_evidence: str | None = _evidence_field(
		"Endogenous Locus Library Method Mechanism"
	)
	endogenous_locus_library_method_mechanism: EndogenousLocusLibraryMethodMechanism | None = Field(
		default=None,
		alias="Endogenous Locus Library Method Mechanism",
		description="Endogenous Locus Library Method Mechanism. Allowed values: Base editor, Nuclease, Other, Prime editor.",
	)
	endogenous_locus_library_method_system_evidence: str | None = _evidence_field(
		"Endogenous Locus Library Method System"
	)
	endogenous_locus_library_method_system: EndogenousLocusLibraryMethodSystem | None = Field(
		default=None,
		alias="Endogenous Locus Library Method System",
		description="Endogenous Locus Library Method System. Allowed values: AsCas12a, Other, RfsCas13d, SaCas9, SpCas9.",
	)
	phenotypic_assay_dimensionality_evidence: str | None = _evidence_field(
		"Phenotypic Assay Dimensionality"
	)
	phenotypic_assay_dimensionality: PhenotypicAssayDimensionality | None = Field(
		default=None,
		alias="Phenotypic Assay Dimensionality",
		description="Phenotypic Assay Dimensionality. Allowed values: Combined functional data, High-dimensional data, Other, Single-dimensional data.",
	)
	phenotypic_assay_mechanism_evidence: str | None = _evidence_field(
		"Phenotypic Assay Mechanism"
	)
	phenotypic_assay_mechanism: PhenotypicAssayMechanism | None = Field(
		default=None,
		alias="Phenotypic Assay Mechanism",
		description="Phenotypic Assay Mechanism. Allowed values: Dominant-negative effect, Gain of function, Loss of function, Mixed functional effect, Other.",
	)
	phenotypic_assay_model_system_evidence: str | None = _evidence_field(
		"Phenotypic Assay Model System"
	)
	phenotypic_assay_model_system: PhenotypicAssayModelSystem | None = Field(
		default=None,
		alias="Phenotypic Assay Model System",
		description="Phenotypic Assay Model System. Allowed values: Bacteria, Bacteriophage, Immortalized human cells, Immortalized non-human mammalian cells, Induced pluripotent stem cells from human female, Induced pluripotent stem cells from human male, Molecular display, Murine primary cells, Not applicable, Organoids, Other, Patient derived primary cells (e.g. T-cells, adipocytes), Yeast.",
	)
	phenotypic_assay_method_evidence: str | None = _evidence_field(
		"Phenotypic Assay Method"
	)
	phenotypic_assay_method: PhenotypicAssayMethod | None = Field(
		default=None,
		alias="Phenotypic Assay Method",
		description="Phenotypic Assay Method. Allowed values: Abundance assay, Binding assay, Bacterial two-hybrid assay, Bulk RNA-sequencing, Cell fitness, Cell morphology assay, Cell proliferation assay, Cell proliferation assay with genetic complementation, Direct protein function, Electrophysiological method, Flow cytometry assay, Fluorescence in-situ hybridization (FISH) assay, Genetic complementation, Homology-directed DNA repair frequency measurement, Imaging mass cytometry assay, Ion channel assay, Multiplexed fluorescent antibody imaging, Not applicable (trained predictor), One-hybrid assay, Other, Polysome profiling assay, Posttranslation modification assay, Protein localization/trafficking assay, Protein stability assay, Reporter gene assay, Single cell imaging, Single-cell RNA sequencing assay, Splicing assay, Survival assessment assay, Systematic evolution of ligands by exponential enrichment assay, Yeast two-hybrid assay.",
	)
	molecular_mechanism_assessed_evidence: str | None = _evidence_field(
		"Molecular Mechanism Assessed"
	)
	molecular_mechanism_assessed: MolecularMechanismAssessed | None = Field(
		default=None,
		alias="Molecular Mechanism Assessed",
		description="Molecular Mechanism Assessed. Allowed values: Calcium-mediated signaling, Catalytic activity and Ornithine carbamoyltransferase activity, Catalytic and Cysteine synthase activity, Catalytic and Gluconokinase activity, Catalytic and Hydroxymethylbilane synthase activity, Catalytic and Lipid phosphatase activity, Catalytic and Thiamine diphosphokinase activity, Cellular response to cisplatin, DNA repair and Double-strand break repair via homologous recombination, Core promoter sequence-specific DNA binding, DNA and Mismatch repair, DNA damage checkpoint signaling, DNA damage response, signal transduction by p53 class mediator and MDM2/MDM4 family protein binding, DNA repair and Double-strand break repair via homologous recombination, Double-strand break repair, damage response, signal transduction by p53 class mediator and Response to etoposide, Molecular condensate scaffold activity, Molecular function, Monoatomic ion transport and Potassium channel activity, Monoatomic ion transport and Sodium ion transport, Oxidative phosphorylation, Protein abundance, Protein carboxylation, Protein glycosylation, Protein localization/trafficking to plasma membrane, Protein-protein interaction, RNA splicing, Regulation of protein stability, Response to misfolded protein, Sodium channel activity, Translation efficiency.",
	)
	in_vitro_construct_library_method_mechanism_evidence: str | None = _evidence_field(
		"In Vitro Construct Library Method Mechanism"
	)
	in_vitro_construct_library_method_mechanism: InVitroConstructLibraryMethodMechanism | None = Field(
		default=None,
		alias="In Vitro Construct Library Method Mechanism",
		description="In Vitro Construct Library Method Mechanism. Allowed values: Episomal delivery, Extra-local construct insertion, Landing pad integration, Native locus replacement, Other, Plasmid (not integrated), Random locus viral integration, Transfection of RNA.",
	)
	in_vitro_construct_library_method_system_evidence: str | None = _evidence_field(
		"In Vitro Construct Library Method System"
	)
	in_vitro_construct_library_method_system: InVitroConstructLibraryMethodSystem | None = Field(
		default=None,
		alias="In Vitro Construct Library Method System",
		description="In Vitro Construct Library Method System. Allowed values: Cassette mutagenesis, Doped oligo synthesis, Error-prone PCR, Microarray synthesis, Nicking mutagenesis, Oligo-directed mutagenic PCR, Oligo pool synthesis, Other, POPCode mutagenesis, Proprietary method, Site-directed mutagenesis.",
	)


__all__ = [
	"CONTROLLED_VOCAB_FIELD_MAPPING",
	"MavedbMetadataExtractionSchema",
	"MavedbMetadataSchema",
]
