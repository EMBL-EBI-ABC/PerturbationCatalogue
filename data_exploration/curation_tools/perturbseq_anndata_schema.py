import pandas as pd
from pandera.pandas import Field, DataFrameModel, dataframe_check
from pandera.typing import Series, Index, String, Int64, Float32
from pathlib import Path


_TREATMENT_TYPE_ID_BY_LABEL = {
    "untreated control": "NCIT:C184729",
    "scrambled control oligonucleotide": "XCO:0001141",
    "culture medium": "BAO:0000114",
    "chemical entity": "CHEBI:24431",
    "protein": "BAO:0000175",
    "protein complex": "BAO:0002554",
    "peptide": "BAO:0000325",
    "antibody": "BAO:0000502",
    "lipid": "BAO:0000171",
    "PNA": "BAO:0000226",
    "DNA": "BAO:0000269",
    "RNA": "BAO:0000270",
    "mRNA": "BAO:0000274",
    "rRNA": "BAO:0000275",
    "tRNA": "BAO:0000276",
    "cDNA": "BAO:0000315",
    "genomic DNA": "BAO:0000316",
    "plasmid DNA": "BAO:0000317",
    "miRNA": "BAO:0000322",
    "shRNA": "BAO:0000323",
    "siRNA": "BAO:0000324",
    "LNA": "BAO:0000412",
    "RNA aptamer": "BAO:0000496",
    "riboswitch": "BAO:0000498",
    "esiRNA": "BAO:0000544",
}

_TREATMENT_TYPE_LABELS = frozenset(_TREATMENT_TYPE_ID_BY_LABEL)
_TREATMENT_TYPE_IDS = frozenset(_TREATMENT_TYPE_ID_BY_LABEL.values())

_DATA_MODALITIES = frozenset({"Perturb-seq", "CRISPR screen", "MAVE"})
_SIGNIFICANCE_VALUES = frozenset({"True", "False"})
_CURATION_AGENT_TYPES = frozenset({"human", "LLM"})
_PERTURBATION_TYPE_LABELS = frozenset({"CRISPRn", "CRISPRi", "CRISPRa", "DMS"})

_MODEL_SYSTEM_ID_BY_LABEL = {
    "cell_line": "CLO:0000031",
    "primary_cell": "BAO:0000239",
    "organoid": "NCIT:C172259",
    "yeast": "NCIT:C19617",
    "bacteria": "NCIT:C19167",
    "bacteriophage": "NCIT:C14188",
    "animal_model": "NCIT:C71164",
    "cell_free_system": "mesh:D002474",
}
_MODEL_SYSTEM_LABELS = frozenset(_MODEL_SYSTEM_ID_BY_LABEL)
_MODEL_SYSTEM_IDS = frozenset(_MODEL_SYSTEM_ID_BY_LABEL.values())

_SPECIES = frozenset({"Homo sapiens"})

_SEX_ID_BY_LABEL = {
    "female": "PATO:0000383",
    "male": "PATO:0000384",
    "hermaphrodite": "PATO:0001340",
}
_SEX_LABELS = frozenset((*_SEX_ID_BY_LABEL, "unknown"))
_SEX_IDS = frozenset(_SEX_ID_BY_LABEL.values())

_DEVELOPMENTAL_STAGE_ID_BY_LABEL = {
    "embryonic": "HsapDv:0000002",
    "fetal": "HsapDv:0000037",
    "neonatal": "HsapDv:0000262",
    "child": "HsapDv:0000265",
    "juvenile": "HsapDv:0000271",
    "adult": "HsapDv:0000258",
    "elderly": "HsapDv:0000227",
}
_DEVELOPMENTAL_STAGE_LABELS = frozenset(_DEVELOPMENTAL_STAGE_ID_BY_LABEL)
_DEVELOPMENTAL_STAGE_IDS = frozenset(_DEVELOPMENTAL_STAGE_ID_BY_LABEL.values())

_LIBRARY_GENERATION_TYPE_ID_BY_LABEL = {
    "endogenous genetic perturbation method": "EFO:0022868",
    "exogenous genetic perturbation method": "EFO:0022869",
}
_LIBRARY_GENERATION_TYPE_LABELS = frozenset(_LIBRARY_GENERATION_TYPE_ID_BY_LABEL)
_LIBRARY_GENERATION_TYPE_IDS = frozenset(_LIBRARY_GENERATION_TYPE_ID_BY_LABEL.values())

_LIBRARY_GENERATION_METHOD_ID_BY_LABEL = {
    "doped oligo synthesis": "EFO:0022900",
    "error-prone PCR": "EFO:0022901",
    "microarray synthesis": "EFO:0022902",
    "silicon microarray synthesis": "EFO:0022902",
    "nicking mutagenesis": "EFO:0022903",
    "oligo-directed mutagenic PCR": "EFO:0022904",
    "site-directed mutagenesis": "EFO:0022905",
    "POPCode mutagenesis": "EFO:0022905",
    "multiplexed site-directed mutagenesis": "EFO:0022905",
    "insertional mutagenesis": "NCIT:C17377",
}
_LIBRARY_GENERATION_METHOD_LABELS = frozenset(
    (
        *_LIBRARY_GENERATION_METHOD_ID_BY_LABEL,
        "solid-phase oligonucleotide synthesis",
        "microchip-based massive parallel oligo synthesis",
        "mutagenesis by integrated tiles",
    )
)
_LIBRARY_GENERATION_METHOD_IDS = frozenset(
    _LIBRARY_GENERATION_METHOD_ID_BY_LABEL.values()
)

_ENZYME_DELIVERY_METHOD_LABELS = frozenset(
    {
        "lipofection",
        "nucleofection",
        "electroporation",
        "adeno-associated virus transduction",
        "adenovirus transduction",
        "retrovirus transduction",
        "lentivirus transduction",
        "nanoparticle-mediated transfection",
        "molecular cloning",
        "transformation",
        "chemical-mediated transfection",
        "hydrodynamic injection",
        "influenza A virus infection",
    }
)
_LIBRARY_DELIVERY_METHOD_LABELS = frozenset(
    {
        "lipofection",
        "nucleofection",
        "electroporation",
        "adeno-associated virus transduction",
        "adenovirus transduction",
        "retrovirus transduction",
        "lentivirus transduction",
        "nanoparticle-mediated transfection",
        "transformation",
        "chemical-mediated transfection",
        "hydrodynamic injection",
        "molecular cloning",
        "influenza A virus infection",
    }
)

_INTEGRATION_STATE_ID_BY_LABEL = {
    "random locus integration": "EFO:0920082",
    "targeted locus integration": "EFO:0920083",
    "native locus replacement": "EFO:0920084",
    "non-integrative transgene expression": "EFO:0920085",
}
_INTEGRATION_STATE_LABELS = frozenset(
    (*_INTEGRATION_STATE_ID_BY_LABEL, "bacteriophage genome integration")
)
_INTEGRATION_STATE_IDS = frozenset(_INTEGRATION_STATE_ID_BY_LABEL.values())

_ENZYME_EXPRESSION_CONTROL_LABELS = frozenset(
    {
        "constitutive transgene expression",
        "inducible transgene expression",
        "native promoter-driven transgene expression",
        "degradation domain-based transgene control",
        "transient transgene expression",
        "minimal promoter-driven transgene expression",
    }
)
_LIBRARY_EXPRESSION_CONTROL_LABELS = _ENZYME_EXPRESSION_CONTROL_LABELS

_LIBRARY_FORMAT_LABELS = frozenset(
    {"pooled", "arrayed", "arrayed|pooled", "in vivo"}
)
_LIBRARY_SCOPE_LABELS = frozenset({"focused", "genome-wide"})

_LIBRARY_PERTURBATION_TYPE_ID_BY_LABEL = {
    "knockout": "EFO:0000506",
    "inhibition": "INO:0000085",
    "activation": "INO:0000075",
    "base editing": "EFO:0022873",
    "prime editing": "EFO:0022872",
    "mutagenesis": "NCIT:C17376",
}
_LIBRARY_PERTURBATION_TYPE_LABELS = frozenset(
    _LIBRARY_PERTURBATION_TYPE_ID_BY_LABEL
)
_LIBRARY_PERTURBATION_TYPE_IDS = frozenset(
    _LIBRARY_PERTURBATION_TYPE_ID_BY_LABEL.values()
)

_READOUT_DIMENSIONALITY_LABELS = frozenset(
    {"single-dimensional assay", "high-dimensional assay"}
)

_READOUT_TYPE_ID_BY_LABEL = {
    "transcriptomic": "EFO:0001032",
    "proteomic": "EFO:0000746",
    "phenotypic": "EFO:0920062",
}
_READOUT_TYPE_LABELS = frozenset(_READOUT_TYPE_ID_BY_LABEL)
_READOUT_TYPE_IDS = frozenset(_READOUT_TYPE_ID_BY_LABEL.values())

_READOUT_TECHNOLOGY_ID_BY_LABEL = {
    "single-cell rna-seq": "EFO:0008913",
    "population growth assay": "EFO:0002907",
    "flow cytometry": "BAO:0000005",
    "high-throughput dna sequencing": "EFO:0002693",
    "patch-clamp electrophysiology": "EFO:0022948",
    "fluorometry": "mesh:D005470",
}
_READOUT_TECHNOLOGY_LABELS = frozenset(_READOUT_TECHNOLOGY_ID_BY_LABEL)
_READOUT_TECHNOLOGY_IDS = frozenset(_READOUT_TECHNOLOGY_ID_BY_LABEL.values())

_READOUT_MEASUREMENT_ID_BY_LABEL = {
    "protein abundance": "BAO:0010252",
    "protein stability": "BAO:0002804",
    "protein activity": "APO:0000022",
    "protein ubiquitination": "GO:0016567",
    "cell viability": "PATO:0000169",
    "cell proliferation": "BAO:0002805",
    "gene expression": "BAO:0002785",
    "RNA splicing": "BAO:0003000",
    "DNA repair": "GO:0006281",
    "ligand binding": "NCIT:C178030",
    "ion channel activity": "BAO:0002997",
    "fluorescence": "BAO:0000363",
}
_READOUT_MEASUREMENT_LABELS = frozenset(
    (
        *_READOUT_MEASUREMENT_ID_BY_LABEL,
        "surface protein expression",
        "viral growth",
    )
)
_READOUT_MEASUREMENT_IDS = frozenset(_READOUT_MEASUREMENT_ID_BY_LABEL.values())

_METHOD_NAME_ID_BY_LABEL = {
    "Perturb-seq": "EFO:0008860",
    "scRNA-seq": "EFO:0008913",
    "pooled growth competition assay": "EFO:0002907",
    "massively parallel reporter assay": "EFO:0008822",
    "yeast surface display": "MI:0115",
    "bacterial two-hybrid assay": "OBI:0001682",
    "mammalian two-hybrid assay": "BAO:0002493",
    "yeast one-hybrid assay": "OBI:0001681",
    "phage display": "MI:0084",
    "mRNA display": "MI:0073",
    "yeast two-hybrid assay": "BAO:0002494",
    "patch-clamp electrophysiology": "EFO:0022948",
    "abundance protein fragment complementation assay": "MI:0090",
    "computational meta-analysis": "NCIT:C17886",
}
_METHOD_NAME_LABELS = frozenset(
    (
        *_METHOD_NAME_ID_BY_LABEL,
        "Perturb-CITE-seq",
        "proliferation CRISPR screen",
        "DMS-TileSeq",
        "DMS-BarSeq",
        "Joined and refined DMS-BarSeq and DMS-TileSeq",
        "Combined DMS-BarSeq and DMS-TileSeq",
        "flow cytometry-based sequencing assay",
        "CRISPR mutagenesis screen",
        "Saturation-Selection-Sequencing assay",
        "fluorescence-based homology-directed repair assay",
        "gap repair assay",
        "homology-directed repair assay",
        "phage-assisted continuous selection",
        "pooled deep mutational scanning",
        "protein folding sensor assay",
        "saturation genome editing",
        "saturation prime editing",
        "saturation base editing",
        "MITE",
        "VAMP-seq",
        "polysome profiling",
    )
)
_METHOD_NAME_IDS = frozenset(_METHOD_NAME_ID_BY_LABEL.values())

_SEQUENCING_LIBRARY_KIT_ID_BY_LABEL = {
    "10x Genomics Single Cell 3-prime v1": "EFO:0009901",
    "10x Genomics Single Cell 3-prime v2": "EFO:0009899",
    "10x Genomics Single Cell 3-prime v3": "EFO:0009922",
    "10x Genomics Single Cell 3-prime v3.1": "EFO:0022980",
    "10x Genomics Chromium GEM-X Flex v1": "EFO:0920088",
}
_SEQUENCING_LIBRARY_KIT_LABELS = frozenset(
    (
        *_SEQUENCING_LIBRARY_KIT_ID_BY_LABEL,
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
    )
)
_SEQUENCING_LIBRARY_KIT_IDS = frozenset(
    _SEQUENCING_LIBRARY_KIT_ID_BY_LABEL.values()
)

_SEQUENCING_PLATFORM_ID_BY_LABEL = {
    "Illumina Genome Analyzer": "EFO:0004200",
    "Illumina Genome Analyzer II": "EFO:0004201",
    "Illumina Genome Analyzer IIx": "EFO:0004202",
    "Illumina HiSeq 2000": "EFO:0004203",
    "Illumina HiSeq 1000": "EFO:0004204",
    "Illumina MiSeq": "EFO:0004205",
    "454 GS 20 sequencer": "EFO:0004206",
    "454 GS sequencer": "EFO:0004431",
    "454 GS FLX sequencer": "EFO:0004432",
    "454 GS FLX Titanium sequencer": "EFO:0004433",
    "454 GS Junior sequencer": "EFO:0004434",
    "AB SOLiD System": "EFO:0004435",
    "AB SOLiD 5500xl": "EFO:0004436",
    "AB SOLiD PI System": "EFO:0004437",
    "AB SOLiD 4 System": "EFO:0004438",
    "AB SOLiD System 3.0": "EFO:0004439",
    "AB SOLiD 5500": "EFO:0004440",
    "AB SOLiD 4hq System": "EFO:0004441",
    "AB SOLiD System 2.0": "EFO:0004442",
    "Illumina HiSeq 4000": "EFO:0008563",
    "Illumina HiSeq 3000": "EFO:0008564",
    "Illumina HiSeq 2500": "EFO:0008565",
    "Illumina NextSeq 550": "EFO:0008566",
    "Illumina HiSeq X": "EFO:0008567",
    "PacBio Sequel system": "EFO:0008630",
    "PacBio RS II": "EFO:0008631",
    "ONT MinION": "EFO:0008632",
    "ONT GridION X5": "EFO:0008633",
    "ONT PromethION": "EFO:0008634",
    "Illumina iSeq 100": "EFO:0008635",
    "Illumina MiniSeq": "EFO:0008636",
    "Illumina NovaSeq 6000": "EFO:0008637",
    "Illumina NextSeq 500": "EFO:0009173",
    "Illumina NextSeq 1000": "EFO:0010962",
    "Illumina NextSeq 2000": "EFO:0010963",
    "Illumina HiSeq 1500": "EFO:0011027",
    "Illumina NovaSeq X": "EFO:0022840",
    "Illumina NovaSeq X Plus": "EFO:0022841",
    "Singular G4": "EFO:0022843",
    "PacBio Sequel II system": "EFO:0700015",
    "BGI MGISEQ-2000": "EFO:0700018",
    "ONT PromethION 2 Solo": "EFO:0700019",
    "Ultima UG100": "EFO:0920005",
    "PacBio Revio": "EFO:0920006",
    "PacBio Onso": "EFO:0920007",
    "Element Aviti": "EFO:0920008",
    "Illumina MiSeq i100": "EFO:0920010",
    "MGI DNBSEQ-T7": "EFO:0920057",
    "Roche 454 GS FLX": "EFO:0004432",
    "Ion Torrent PGM": "GENEPIO:0100136",
    "Ultima Genomics UG100": "EFO:0920005",
}
_SEQUENCING_PLATFORM_LABELS = frozenset(
    (
        *_SEQUENCING_PLATFORM_ID_BY_LABEL,
        "Roche 454 GS FLX+",
        "Illumina NextSeq (model unspecified)",
        "Illumina sequencer (model unspecified)",
        "Illumina HiSeq (model unspecified)",
        "PacBio sequencer (model unspecified)",
    )
)
_SEQUENCING_PLATFORM_IDS = frozenset(_SEQUENCING_PLATFORM_ID_BY_LABEL.values())

_SEQUENCING_STRATEGY_ID_BY_LABEL = {"direct sequencing": "NCIT:C116154"}
_SEQUENCING_STRATEGY_LABELS = frozenset(
    {
        "barcode sequencing",
        "direct sequencing",
        "barcode sequencing|direct sequencing",
    }
)
_SEQUENCING_STRATEGY_IDS = frozenset(_SEQUENCING_STRATEGY_ID_BY_LABEL.values())

_SOFTWARE_COUNTS_LABELS = frozenset(
    {
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
        "ABSSeq",
        "satmut_utils",
        "bcftools",
        "TagDust2",
        "Subassembly",
        "pysamstats",
        "Jellyfish",
        "Tagdust2",
    }
)
_SOFTWARE_ANALYSIS_LABELS = frozenset(
    {
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
        "dms_variants",
        "mapmuts",
        "maveLLR",
        "tileseq_package",
        "tileseqMave",
        "TileseqMave",
        "ABSSeq",
        "Cluster",
        "ORFcall",
        "samtools",
    }
)
_REFERENCE_GENOME_LABELS = frozenset(
    {
        "GRCh38",
        "GRCh37",
        "cDNA reference sequence",
        "mm9",
        "S288c",
        "hg19",
        "Wuhan-Hu-1",
        "non-standard reference sequence",
    }
)

_LICENSE_ID_BY_LABEL = {
    "CC0": "SWO:1000049",
    "CC BY": "SWO:1000050",
    "CC BY-SA": "SWO:1000052",
    "CC BY-NC": "SWO:1000079",
    "CC BY-ND": "SWO:1000077",
    "CC0 1.0": "SWO:1000049",
    "CC BY 2.0": "SWO:1000050",
    "CC BY-SA 2.0": "SWO:1000052",
    "CC BY 4.0": "SWO:1000065",
    "CC BY 2.0 UK": "SWO:1000067",
    "CC BY 2.1 JP": "SWO:1000072",
    "CC BY 2.5": "SWO:1000073",
    "CC BY 3.0 AU": "SWO:1000074",
    "CC BY 3.0": "SWO:1000075",
    "CC BY 3.0 US": "SWO:1000076",
    "CC BY-ND 3.0": "SWO:1000077",
    "CC BY-ND 4.0": "SWO:1000078",
    "CC BY-NC 3.0": "SWO:1000079",
    "CC BY-NC 4.0": "SWO:1000080",
    "CC BY-NC-ND 3.0": "SWO:1000081",
    "CC BY-NC-ND 2.5": "SWO:1000083",
    "CC BY-NC-ND 2.5 CH": "SWO:1000084",
    "CC BY-NC-ND 4.0": "SWO:1000085",
    "CC BY-NC-SA 2.5": "SWO:1000086",
    "CC BY-NC-SA 3.0": "SWO:1000087",
    "CC BY-NC-SA 3.0 US": "SWO:1000088",
    "CC BY-NC-SA 2.5 IN": "SWO:1000089",
    "CC BY-NC-SA 4.0": "SWO:1000090",
    "CC BY-SA 2.1 JP": "SWO:1000091",
    "CC BY-SA 3.0": "SWO:1000092",
    "CC BY-SA 3.0 US": "SWO:1000093",
    "CC BY-SA 4.0": "SWO:1000094",
}
_LICENSE_LABELS = frozenset(_LICENSE_ID_BY_LABEL)
_LICENSE_IDS = frozenset(_LICENSE_ID_BY_LABEL.values())

_TREATMENT_UNITS = frozenset(
    {
        "pM",
        "nM",
        "uM",
        "mM",
        "M",
        "pg/mL",
        "ng/mL",
        "ug/mL",
        "mg/mL",
        "g/mL",
        "pg/kg",
        "ug/kg",
        "mg/kg",
        "g/kg",
        "cells/uL",
        "cells/mL",
        "MOI",
        "uL",
        "mL",
        "%",
        "IU/mL",
    }
)

_TREATMENT_FIELDS = (
    "treatment_type_label",
    "treatment_type_id",
    "treatment_label",
    "treatment_id",
    "treatment_dose",
    "treatment_unit",
)


def _split_pipe_value(value: object) -> list[str] | None:
    if pd.isna(value):
        return None
    return [token.strip() for token in str(value).split("|")]


def _value_has_allowed_tokens(value: object, allowed_values: frozenset[str]) -> bool:
    tokens = _split_pipe_value(value)
    return tokens is None or all(token in allowed_values for token in tokens)


def _series_has_allowed_tokens(
    values: pd.Series, allowed_values: frozenset[str]
) -> pd.Series:
    string_values = values.astype("string")
    is_pipe_delimited = string_values.str.contains("|", regex=False, na=False)
    result = values.isna() | (
        ~is_pipe_delimited & string_values.isin(allowed_values)
    )

    if is_pipe_delimited.any():
        result.loc[is_pipe_delimited] = string_values.loc[
            is_pipe_delimited
        ].map(lambda value: _value_has_allowed_tokens(value, allowed_values))

    return result


def _token_count(value: object) -> int | None:
    tokens = _split_pipe_value(value)
    return None if tokens is None else len(tokens)


def _treatment_type_values_correspond(label_value: object, id_value: object) -> bool:
    labels = _split_pipe_value(label_value)
    ids = _split_pipe_value(id_value)

    if labels is None or ids is None:
        return labels is None and ids is None

    return len(labels) == len(ids) and all(
        _TREATMENT_TYPE_ID_BY_LABEL.get(label) == treatment_id
        for label, treatment_id in zip(labels, ids)
    )


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
        isin=_DATA_MODALITIES,
    )
    significant: Series[String] = Field(
        nullable=True,
        description="Indicates whether the perturbation had a significant effect.",
        coerce=True,
        isin=_SIGNIFICANCE_VALUES,
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
        isin=_PERTURBATION_TYPE_LABELS,
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
    )
    treatment_type_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID for the treatment type.",
    )
    treatment_label: Series[String] = Field(
        nullable=True,
        description="Treatment/compound ontology term label used to stimulate the investigated sample. ChEMBL compound label for chemical entities. Use 'untreated control' for untreated samples where other samples were treated.",
    )
    treatment_id: Series[String] = Field(
        nullable=True,
        description="Treatment/compound ontology term ID used to stimulate the investigated sample. ChEMBL compound ID.",
    )
    treatment_dose: Series[String] = Field(
        nullable=True,
        coerce=True,
        ignore_na=True,
        description="Treatment/compound dose used to stimulate the investigated sample.",
    )
    treatment_unit: Series[String] = Field(
        nullable=True,
        description="Treatment/compound unit used to stimulate the investigated sample. Use 'u' for micro (e.g., 'uM' instead of 'μM').",
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
        isin=_MODEL_SYSTEM_LABELS,
    )
    
    model_system_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Model system ontology term ID of the investigated sample.",
        isin=_MODEL_SYSTEM_IDS,
    )
    
    species: Series[String] = Field(
        nullable=False,
        description="Species name of the investigated sample.",
        isin=_SPECIES,
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
        isin=_SEX_LABELS,
    )
    sex_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Sex ontology term ID of the investigated sample.",
        isin=_SEX_IDS,
    )
    developmental_stage_label: Series[String] = Field(
        nullable=True,
        description="Developmental stage ontology term label of the investigated sample.",
        isin=_DEVELOPMENTAL_STAGE_LABELS,
    )
    developmental_stage_id: Series[String] = Field(
        nullable=True,
        str_contains=":",
        description="Developmental stage ontology term ID of the investigated sample.",
        isin=_DEVELOPMENTAL_STAGE_IDS,
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
        isin=_LIBRARY_GENERATION_TYPE_LABELS,
    )
    library_generation_type_id: Series[String] = Field(
        nullable=True,
        description="Library generation type ontology term ID, defined in EFO under parent term EFO:0022867 (genetic perturbation)",
        isin=_LIBRARY_GENERATION_TYPE_IDS,
    )
    library_generation_method_label: Series[String] = Field(
        nullable=True,
        description="Library generation method ontology term label, defined in EFO under parent term EFO:0022868/EFO:0022869 (Endogenous/Exogenous genetic perturbation method)",
        isin=_LIBRARY_GENERATION_METHOD_LABELS,
    )
    library_generation_method_id: Series[String] = Field(
        nullable=True,
        description="Library generation method ontology term ID, defined in EFO under parent term EFO:0022868/EFO:0022869 (Endogenous/Exogenous genetic perturbation method)",
        isin=_LIBRARY_GENERATION_METHOD_IDS,
    )
    enzyme_delivery_method_label: Series[String] = Field(
        nullable=True,
        description="Enzyme delivery method ontology term label.",
        isin=_ENZYME_DELIVERY_METHOD_LABELS,
    )
    enzyme_delivery_method_id: Series[String] = Field(
        nullable=True,
        description="Enzyme delivery method ontology term ID.",
    )
    library_delivery_method_label: Series[String] = Field(
        nullable=True,
        description="Library delivery method ontology term label.",
        isin=_LIBRARY_DELIVERY_METHOD_LABELS,
    )
    library_delivery_method_id: Series[String] = Field(
        nullable=True, description="Library delivery method ontology term ID."
    )
    enzyme_integration_state_label: Series[String] = Field(
        nullable=True,
        description="Enzyme integration state ontology term label.",
        isin=_INTEGRATION_STATE_LABELS,
    )
    enzyme_integration_state_id: Series[String] = Field(
        nullable=True, description="Enzyme integration state ontology term ID.",
        isin=_INTEGRATION_STATE_IDS,
    )
    library_integration_state_label: Series[String] = Field(
        nullable=True,
        description="Library integration state ontology term label.",
        isin=_INTEGRATION_STATE_LABELS,
    )
    library_integration_state_id: Series[String] = Field(
        nullable=True, description="Library integration state ontology term ID.",
        isin=_INTEGRATION_STATE_IDS,
    )
    enzyme_expression_control_label: Series[String] = Field(
        nullable=True,
        description="Enzyme expression control ontology term label.",
        isin=_ENZYME_EXPRESSION_CONTROL_LABELS,
    )
    enzyme_expression_control_id: Series[String] = Field(
        nullable=True, description="Enzyme expression control ontology term ID."
    )
    # library details
    library_expression_control_label: Series[String] = Field(
        nullable=True,
        description="Library expression control ontology term label.",
        isin=_LIBRARY_EXPRESSION_CONTROL_LABELS,
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
        isin=_LIBRARY_FORMAT_LABELS,
    )
    library_format_id: Series[String] = Field(
        nullable=True, description="Perturbation library format ontology term ID."
    )
    library_scope_label: Series[String] = Field(
        nullable=True,
        description="Perturbation library scope ontology term label.",
        isin=_LIBRARY_SCOPE_LABELS,
    )
    library_scope_id: Series[String] = Field(
        nullable=True, description="Perturbation library scope ontology term ID."
    )
    library_perturbation_type_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the library perturbation type.",
        isin=_LIBRARY_PERTURBATION_TYPE_LABELS,
    )
    library_perturbation_type_id: Series[String] = Field(
        nullable=True, description="Ontology term ID for the library perturbation type.",
        isin=_LIBRARY_PERTURBATION_TYPE_IDS,
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
        isin=_READOUT_DIMENSIONALITY_LABELS,
    )
    readout_dimensionality_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the dimensionality of the readout assay.",
    )
    readout_type_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the type of the readout assay.",
        isin=_READOUT_TYPE_LABELS,
    )
    readout_type_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the type of the readout assay.",
        isin=_READOUT_TYPE_IDS,
    )
    readout_technology_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the technology used in the readout assay.",
        isin=_READOUT_TECHNOLOGY_LABELS,
    )
    readout_technology_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the technology used in the readout assay.",
        isin=_READOUT_TECHNOLOGY_IDS,
    )
    readout_measurement_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the measurement type of the readout assay.",
        isin=_READOUT_MEASUREMENT_LABELS,
    )
    readout_measurement_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the measurement type of the readout assay.",
        isin=_READOUT_MEASUREMENT_IDS,
    )
    method_name_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the method name used in the readout assay.",
        isin=_METHOD_NAME_LABELS,
    )
    method_name_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the method name used in the readout assay.",
        isin=_METHOD_NAME_IDS,
    )
    method_uri: Series[String] = Field(
        nullable=True,
        description="URI associated with the method used in the readout assay.",
    )
    sequencing_library_kit_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing library kit.",
        isin=_SEQUENCING_LIBRARY_KIT_LABELS,
    )
    sequencing_library_kit_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing library kit.",
        isin=_SEQUENCING_LIBRARY_KIT_IDS,
    )
    sequencing_platform_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing platform.",
        isin=_SEQUENCING_PLATFORM_LABELS,
    )
    sequencing_platform_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing platform.",
        isin=_SEQUENCING_PLATFORM_IDS
    )
    sequencing_strategy_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label associated with the sequencing strategy.",
        isin=_SEQUENCING_STRATEGY_LABELS,
    )
    sequencing_strategy_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID associated with the sequencing strategy.",
        isin=_SEQUENCING_STRATEGY_IDS,
    )
    software_counts_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the software used for generating counts.",
        isin=_SOFTWARE_COUNTS_LABELS,
    )
    software_counts_id: Series[String] = Field(
        nullable=True,
        description="Ontology term ID for the software used for generating counts.",
    )
    software_analysis_label: Series[String] = Field(
        nullable=True,
        description="Ontology term label for the software used for analysis.",
        isin=_SOFTWARE_ANALYSIS_LABELS,
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
        isin=_REFERENCE_GENOME_LABELS,
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
        isin=_LICENSE_LABELS,
    )
    license_id: Series[String] = Field(
        nullable=True,
        description="License ontology term ID for data usage and distribution. Should be one of the terms from under SWO:0000002 (license).",
        isin=_LICENSE_IDS
    )
    curation_agent_type: Series[String] = Field(
        nullable=False,
        description="Type of agent that curated this dataset: 'human' for manual curation, 'LLM' for automated curation by a language model.",
        isin=_CURATION_AGENT_TYPES,
    )
    curation_agent_name: Series[String] = Field(
        nullable=False,
        description="Name or identifier of the curator. For humans: full name (e.g., 'John Doe'). For LLMs: model identifier (e.g., 'google/gemini-3.5-flash').",
    )

    # Checks
    @dataframe_check(
        ignore_na=False,
        error="cell_barcode is required for Perturb-seq rows.",
    )
    def perturbseq_requires_cell_barcode(cls, df: pd.DataFrame) -> pd.Series:
        is_perturbseq = df["data_modality"].eq("Perturb-seq")
        return ~is_perturbseq | df["cell_barcode"].notna()

    @dataframe_check(
        error="Each treatment_type_label must be an allowed term, including within pipe-delimited values.",
    )
    def treatment_type_labels_are_valid(cls, df: pd.DataFrame) -> pd.Series:
        return _series_has_allowed_tokens(
            df["treatment_type_label"], _TREATMENT_TYPE_LABELS
        )

    @dataframe_check(
        error="Each treatment_type_id must be an allowed term, including within pipe-delimited values.",
    )
    def treatment_type_ids_are_valid(cls, df: pd.DataFrame) -> pd.Series:
        return _series_has_allowed_tokens(
            df["treatment_type_id"], _TREATMENT_TYPE_IDS
        )

    @dataframe_check(
        error="Each treatment_type_label must correspond to its treatment_type_id, including within pipe-delimited values.",
    )
    def treatment_type_labels_and_ids_correspond(cls, df: pd.DataFrame) -> pd.Series:
        return df[["treatment_type_label", "treatment_type_id"]].apply(
            lambda row: _treatment_type_values_correspond(
                row["treatment_type_label"], row["treatment_type_id"]
            ),
            axis=1,
        )

    @dataframe_check(
        error="Each treatment_unit must be an allowed unit, including within pipe-delimited values.",
    )
    def treatment_units_are_valid(cls, df: pd.DataFrame) -> pd.Series:
        return _series_has_allowed_tokens(
            df["treatment_unit"], _TREATMENT_UNITS
        )

    @dataframe_check(
        ignore_na=False,
        error="treatment_dose must be a float or a pipe-delimited list of floats.",
    )
    def treatment_dose_values_are_floats(cls, df: pd.DataFrame) -> pd.Series:
        def is_valid_float_list(value: object) -> bool:
            tokens = _split_pipe_value(value)
            if tokens is None:
                return True
            for token in tokens:
                try:
                    float(token)
                except (TypeError, ValueError):
                    return False
            return True

        return df["treatment_dose"].map(is_valid_float_list)

    @dataframe_check(
        error="treatment_type_label and treatment_type_id must either both be present or both be absent.",
    )
    def treatment_type_label_and_id_are_paired(cls, df: pd.DataFrame) -> pd.Series:
        return df["treatment_type_label"].notna() == df["treatment_type_id"].notna()

    @dataframe_check(
        error="treatment_type_label and treatment_label must either both be present or both be absent.",
    )
    def treatment_type_label_and_treatment_label_are_paired(
        cls, df: pd.DataFrame
    ) -> pd.Series:
        return df["treatment_type_label"].notna() == df["treatment_label"].notna()

    @dataframe_check(
        error="Pipe-delimited treatment fields must contain the same number of values.",
    )
    def treatment_fields_are_aligned(cls, df: pd.DataFrame) -> pd.Series:
        string_values = df[list(_TREATMENT_FIELDS)].astype("string")
        has_pipe_delimited_values = string_values.apply(
            lambda column: column.str.contains("|", regex=False, na=False)
        ).any(axis=1)
        result = pd.Series(True, index=df.index)

        if has_pipe_delimited_values.any():
            token_counts = string_values.loc[has_pipe_delimited_values].map(
                _token_count
            )
            result.loc[has_pipe_delimited_values] = token_counts.notna().all(
                axis=1
            ) & token_counts.nunique(axis=1).eq(1)

        return result

    @dataframe_check(
        error="Chemical entity treatments require a treatment_label.",
    )
    def chemical_entity_requires_treatment_id(cls, df: pd.DataFrame) -> pd.Series:
        treatment_type_labels = df["treatment_type_label"].astype("string")
        treatment_labels = df["treatment_label"].astype("string")
        is_pipe_delimited = treatment_type_labels.str.contains(
            "|", regex=False, na=False
        )
        is_chemical_entity = treatment_type_labels.eq("chemical entity").fillna(False)
        treatment_label_is_present = treatment_labels.notna()
        result = ~is_chemical_entity | treatment_label_is_present

        if not is_pipe_delimited.any():
            return result

        def row_is_valid(row: pd.Series) -> bool:
            treatment_type_labels = _split_pipe_value(row["treatment_type_label"])
            if treatment_type_labels is None:
                return True

            chemical_entity_positions = [
                index
                for index, label in enumerate(treatment_type_labels)
                if label == "chemical entity"
            ]
            if not chemical_entity_positions:
                return True

            treatment_labels = _split_pipe_value(row["treatment_label"])
            if treatment_labels is None or len(treatment_labels) != len(
                treatment_type_labels
            ):
                return False

            return all(treatment_labels[index] for index in chemical_entity_positions)

        result.loc[is_pipe_delimited] = df.loc[
            is_pipe_delimited, ["treatment_type_label", "treatment_label"]
        ].apply(row_is_valid, axis=1)
        return result

    @dataframe_check(
        error="treatment_dose and treatment_unit must either both be present or both be absent.",
    )
    def treatment_dose_and_unit_are_paired(cls, df: pd.DataFrame) -> pd.Series:
        return df["treatment_dose"].notna() == df["treatment_unit"].notna()

    @dataframe_check(
        ignore_na=False,
        error="perturbed_target_coord is required when perturbed_target_biotype is enhancer.",
    )
    def enhancer_requires_target_coord(cls, df: pd.DataFrame) -> pd.Series:
        is_enhancer = df["perturbed_target_biotype"].eq("enhancer")
        return ~is_enhancer | df["perturbed_target_coord"].notna()

    @dataframe_check(
        error="If model system is cell_line, then cell_line_label must be present.",
    )
    def cell_line_requires_metadata(cls, df: pd.DataFrame) -> pd.Series:
        is_cell_line = df["model_system_label"].eq("cell_line")
        return ~is_cell_line | df["cell_line_label"].notna()

    @dataframe_check(
        error="If model system is anything other than cell_line, then cell_line_label and cell_line_id must be absent.",
    )
    def non_cell_line_excludes_metadata(cls, df: pd.DataFrame) -> pd.Series:
        is_cell_line = df["model_system_label"].eq("cell_line")
        cell_line_metadata_absent = df["cell_line_label"].isna() & df[
            "cell_line_id"
        ].isna()
        return is_cell_line | cell_line_metadata_absent

    @dataframe_check(
        error="Cell line label must be present only when model system is cell line.",
    )
    def cell_line_label_requires_cell_line_model(
        cls, df: pd.DataFrame
    ) -> pd.Series:
        return ~df["cell_line_label"].notna() | df["model_system_label"].eq(
            "cell_line"
        )

    @dataframe_check(
        error="Cell line ID must be present only when model system is cell line.",
    )
    def cell_line_id_requires_cell_line_model(cls, df: pd.DataFrame) -> pd.Series:
        return ~df["cell_line_id"].notna() | df["model_system_label"].eq(
            "cell_line"
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
