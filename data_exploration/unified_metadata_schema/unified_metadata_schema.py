import json
from pydantic import (
    BaseModel,
    model_validator,
    HttpUrl,
    Field,
)
from typing import Optional, List, Literal
from fastapi import FastAPI, HTTPException
from enum import Enum
from datetime import datetime, timedelta


# Load the options from the JSON file
with open("enums.json", "r") as f:
    options = json.load(f)

# Load the CRISPR libraries from the JSON file
with open("crispr_libraries.json", "r") as f:
    crispr_libraries = json.load(f)


# Function to make enums
def make_enum(name, options_dict=options):
    values = options_dict[name]
    return Enum(name, {i.replace(" ", "_").lower(): i for i in values})


# Define the models
class Author(BaseModel):
    first_name: str = Field(..., description="First name of the author", example="John")
    last_name: str = Field(..., description="Last name of the author", example="Doe")


class StudyDetails(BaseModel):
    title: str = Field(..., description="Title of the study/publication")
    study_uri: Optional[HttpUrl] = Field(
        None, description="URI/link of the study/publication"
    )
    year: int = Field(
        ge=1950,
        le=datetime.now().year,
        description="Year of the study/publication",
        example=2024,
    )
    first_author: Author = Field(
        ..., description="First author of the study/publication"
    )
    last_author: Author = Field(..., description="Last author of the study/publication")


class SampleQuantity(BaseModel):
    sample_quantity_value: float = Field(
        ..., description="Sample quantity value", example=1.0
    )
    sample_quantity_unit: make_enum("SampleQuantityUnit") = Field(
        ..., description="Sample quantity unit", example="gram"
    )


class Treatment(BaseModel):
    term_id: str = Field(
        ...,
        description="Term ID in CURIE format of the treatment defined in ChEMBL",
        example="CHEBI:6904",
    )
    term_label: str = Field(
        ...,
        description="Label of the treatment defined in ChEMBL",
        example="metoprolol",
    )


class ExperimentDetails(BaseModel):
    title: str = Field(..., description="Title of the experiment")
    summary: str = Field(..., description="Short summary of the experiment")
    treatments: Optional[List[Treatment]] = Field(
        None, description="List of treatments used in the experiment, defined in ChEMBL"
    )
    timepoints: Optional[List[timedelta]] = Field(
        None,
        description="List of timepoints captured in the experiment. Must be in the format of a timedelta object.",
        example="P1DT6H",
    )
    replicates: make_enum("Replicates") = Field(
        ...,
        description="Types of replicates used the experiment: Biological, Technical or Biological and Technical",
        example="Biological",
    )
    number_of_samples: int = Field(
        ..., ge=1, description="Number of samples in the experiment", example=3
    )
    number_of_perturbed_cells: int = Field(
        ...,
        ge=1,
        description="Number of perturbed cells profiled in the experiment",
        example=200000,
    )
    number_of_perturbed_genes: int = Field(
        ...,
        ge=1,
        description="How many genes have been perturbed in the experiment",
        example=4,
    )


class Library(BaseModel):
    library_name: str = Field(
        ...,
        description="Name of the library used in the experiment",
        example="Bassik Human CRISPR Knockout Library",
    )
    accession: Optional[str] = Field(
        None,
        description="Accession number of the library used in the experiment (if available)",
    )
    library_format: make_enum("LibraryFormat") = Field(
        ...,
        description="Format of the library used in the experiment (Pooled or Arrayed)",
        example="Pooled",
    )
    library_scope: make_enum("LibraryScope") = Field(
        ...,
        description="Scope of the library used in the experiment (Genome-wide or Focused)",
        example="Genome-wide",
    )
    perturbation_type: make_enum("PerturbationType") = Field(
        ...,
        description="Type of perturbation used in the experiment (e.g. Knockout, Inhibition, Base editing)",
        example="Knockout",
    )
    manufacturer: str = Field(
        ...,
        description="The lab or commercial manufacturer of the library used in the experiment",
        example="Bassik",
    )
    lentiviral_generation: str = Field(
        ...,
        description="Lentiviral generation of the library used in the experiment",
        example="3",
    )
    grnas_per_gene: Optional[str] = Field(
        None,
        description="Number of gRNAs targeting each gene in the library. Can be a string representing a range, e.g. 3-5.",
        example=5,
    )
    total_grnas: Optional[str] = Field(
        None,
        description="Total number of gRNAs in the library. Can be a string (e.g. 'Varies')",
        example="10000",
    )
    total_genes: Optional[int] = Field(
        None,
        ge=1,
        description="Total number of genes targeted by the library",
        example=18000,
    )
    total_variants: Optional[int] = Field(
        None,
        ge=1,
        description="For SGE experiments, total number of variants in the library",
        example=500,
    )

    # Update the library parameters with the pre-defined values from the "crispr_libraries" dictionary
    @model_validator(mode="before")
    @classmethod
    def validate_enum_library(cls, values):
        lib_name = values["library_name"]
        if lib_name in crispr_libraries["libraries"].keys():
            # If the library name is in the crispr_libraries, update the values with the pre-defined values
            lib = crispr_libraries["libraries"][lib_name]
            values.update(lib)
            print(
                "Library name is in the crispr_libraries, updated values with pre-defined values"
            )
        return values

    # If perturbation_type is "Saturation mutagenesis", then total_variants is required
    @model_validator(mode="before")
    @classmethod
    def validate_total_variants(cls, values):
        if (
            values.get("perturbation_type") == "Saturation mutagenesis"
            and values.get("total_variants") is None
        ):
            raise ValueError("Total variants is required for saturation mutagenesis")
        if (
            values.get("perturbation_type") != "Saturation mutagenesis"
            and values.get("total_variants") is not None
        ):
            raise ValueError(
                "Total variants is not required for this perturbation type"
            )
        if (
            values.get("perturbation_type") == "Saturation mutagenesis"
            and values.get("library_scope") == "Genome-wide"
        ):
            raise ValueError(
                "Saturation mutagenesis perturbation type is not allowed for Genome-wide libraries"
            )
        return values


class PerturbationDetails(BaseModel):
    library_generation_type: Literal["Endogenous", "Exogenous"] = Field(
        ...,
        description="Whether the genetic modifications are primarily driven by intracellular mechanisms (Endogenous) or by external in vitro laboratory techniques (Exogenous)",
        example="Exogenous",
    )
    library_generation_method: str = Field(
        ...,
        description="Method of library generation (can be an enzyme, such as 'SpCas9', or a methodology, such as 'doped oligo synthesis'). Must be an EFO term in CURIE format under 'Genetic perturbation' EFO:0022867 parent",
        example="EFO:0022876",
    )
    enzyme_delivery_method: Optional[make_enum("DeliveryMethod")] = Field(
        None,
        description="Delivery method of the enzyme used in the experiment (e.g. Lentiviral transduction, Electroporation, Lipofection)",
        example="Electroporation",
    )
    library_delivery_method: make_enum("DeliveryMethod") = Field(
        ...,
        description="Delivery method of the library used in the experiment (e.g. Lentiviral transduction, Electroporation, Lipofection)",
        example="Lentiviral transduction",
    )
    enzyme_integration_state: Optional[make_enum("IntegrationState")] = Field(
        None,
        description="Integration state of the enzyme used in the experiment (e.g. Episomal expression, Random locus integration)",
        example="Random locus integration",
    )
    library_integration_state: make_enum("IntegrationState") = Field(
        ...,
        description="Integration state of the library used in the experiment (e.g. Episomal expression, Random locus integration)",
        example="Episomal expression",
    )
    enzyme_expression_control: Optional[make_enum("ExpressionControl")] = Field(
        None,
        description="Expression control of the enzyme used in the experiment (e.g. Native promoter-driven expression, Inducible)",
        example="Inducible",
    )
    library_expression_control: make_enum("ExpressionControl") = Field(
        ...,
        description="Expression control of the library used in the experiment (e.g. Native promoter-driven expression, Inducible)",
        example="Native promoter-driven expression",
    )
    library: Library

    @model_validator(mode="before")
    @classmethod
    def validate_enzyme_delivery_method(cls, values):
        enzyme_delivery_method = values.get("enzyme_delivery_method")
        library_generation_type = values.get("library_generation_type")
        if library_generation_type == "Exogenous" and enzyme_delivery_method != None:
            raise ValueError(
                "Enzyme delivery method is not required for exogenous library generation type"
            )
        if library_generation_type == "Endogenous" and enzyme_delivery_method == None:
            raise ValueError(
                "Enzyme delivery method is required for endogenous library generation type"
            )
        return values

    @model_validator(mode="before")
    @classmethod
    def validate_enzyme_expression_control(cls, values):
        enzyme_expression_control = values.get("enzyme_expression_control")
        library_generation_type = values.get("library_generation_type")
        if library_generation_type == "Exogenous" and enzyme_expression_control != None:
            raise ValueError(
                "Enzyme expression control is not required for exogenous library generation type"
            )
        if (
            library_generation_type == "Endogenous"
            and enzyme_expression_control == None
        ):
            raise ValueError(
                "Enzyme expression control is required for endogenous library generation type"
            )
        return values


class AssayDetails(BaseModel):
    readout_dimensionality: make_enum("ReadoutDimensionality") = Field(
        ...,
        description="Dimensionality of the readout (i.e Single-dimensional assay, High-dimensional assay)",
        example="Single-dimensional assay",
    )
    readout_type: make_enum("ReadoutType") = Field(
        ...,
        description="Type of the readout (e.g. Transcriptomic, Phenotypic)",
        example="Transcriptomic",
    )
    readout_technology: make_enum("ReadoutTechnology") = Field(
        ...,
        description="Technology used for the readout (e.g. Population growth assay, Flow cytometry)",
        example="Population growth assay",
    )
    method_name: make_enum("MethodName") = Field(
        ...,
        description="Name of the method used for the readout (e.g. scRNA-seq, Proliferation CRISPR screen)",
        example="Proliferation CRISPR screen",
    )
    method_uri: Optional[HttpUrl] = None
    sequencing_library_kit: make_enum("SequencingLibraryKits") = Field(
        ...,
        description="Sequencing library kit used in the experiment (e.g. 10x Genomics Single Cell 3-prime v1)",
        example="10x Genomics Single Cell 3-prime v1",
    )
    sequencing_platform: make_enum("SequencingPlatform") = Field(
        ...,
        description="Sequencing platform used in the experiment (e.g. Illumina NovaSeq 6000, ONT MinION)",
        example="Illumina NovaSeq 6000",
    )
    sequencing_strategy: make_enum("SequencingStrategy") = Field(
        ...,
        description="Sequencing strategy used in the experiment (e.g. Shotgun sequencing, Barcode sequencing)",
        example="Shotgun sequencing",
    )
    software_counts: List[make_enum("SoftwareCounts")] = Field(
        ...,
        description="Software used for counting the reads (e.g. CellRanger, MAGeCK)",
        example="CellRanger",
    )
    software_analysis: List[make_enum("SoftwareAnalysis")] = Field(
        ...,
        description="Software used for analyzing the data (e.g. Seurat, MAGeCK)",
        example="Seurat",
    )
    reference_genome: make_enum("ReferenceGenome") = Field(
        ...,
        description="Reference genome used in the experiment (e.g. GRCh38, T2T-CHM13)",
        example="GRCh38",
    )


class ModelSystemDetails(BaseModel):
    model_system: make_enum("ModelSystem") = Field(
        ...,
        description="Model system used in the experiment (e.g. Cell line, Primary cell, Tissue sample, Whole organism)",
        example="Cell line",
    )
    species: Literal["Homo Sapiens"] = Field(
        ..., description="Species used in the experiment", example="Homo sapiens"
    )
    tissue: Optional[List[str]] = Field(
        None,
        description="Specific biological tissue the sample is derived from. Must be a term ID in CURIE format from UBERON anatomical entity; child of UBERON:0001062.",
        example="UBERON:0002098",
    )
    cell_type: Optional[List[str]] = Field(
        None,
        description="Cell type/types profiled in the experiment. Must be a term ID in CURIE format from Cell Ontology; child of CL:0000000.",
        example="CL:0000235",
    )
    cell_line: Optional[str] = Field(
        None,
        description="Cell line name used in the experiment. Must be a term ID in CURIE format from Cell Line Ontology 'cultured cell' CL:0000010 parent",
        example="CLO:0009348",
    )
    sex: make_enum("Sex") = Field(
        ..., description="Model system organism's sex", example="Unknown"
    )
    developmental_stage: Optional[List[make_enum("DevelopmentalStages")]] = Field(
        None, description="Developmental stage/age of the model system", example="Adult"
    )
    passage_number: Optional[int] = Field(
        None,
        ge=1,
        description="Passage number of cultured cells (if known)",
        example=23,
    )
    sample_quantity: SampleQuantity = Field(
        ..., description="Sample quantity used in the experiment"
    )

    @model_validator(mode="before")
    @classmethod
    def validate_enum_model_system(cls, values):
        model_system = values["model_system"]

        if model_system == "Tissue sample":
            if values.get("tissue") is None:
                raise ValueError("Tissue is required for tissue sample model system")
            if values.get("cell_type") is not None:
                raise ValueError(
                    "Cell type is not required for tissue sample model system"
                )
            if values.get("cell_line") is not None:
                raise ValueError(
                    "Cell line name is not required for tissue sample model system"
                )
            if values.get("passage_number") is not None:
                raise ValueError(
                    "Passage number is not required for tissue sample model system"
                )
            if values.get("developmental_stage") is None:
                raise ValueError(
                    "Developmental stage is required for tissue sample model system"
                )

        if model_system == "Cell line":
            if values.get("cell_type") is None:
                raise ValueError("Cell type is required for cell line model system")
            if values.get("passage_number") is None:
                raise ValueError(
                    "Passage number is required for cell line model system"
                )
            if values.get("cell_line") is None:
                raise ValueError(
                    "Cell line name is required for cell line model system"
                )
            if values.get("tissue") is not None:
                raise ValueError("Tissue is not required for cell line model system")
            # if values.get("developmental_stage") is not None:
            #     raise ValueError("Developmental stage is not required for cell line model system")

        if model_system == "Primary cell":
            if values.get("cell_type") is None:
                raise ValueError("Cell type is required for primary cell model system")
            if values.get("passage_number") is None:
                raise ValueError(
                    "Passage number is required for primary cell model system"
                )
            if values.get("cell_line") is not None:
                raise ValueError(
                    "Cell line name is not required for primary cell model system"
                )
            if values.get("tissue") is not None:
                raise ValueError("Tissue is not required for primary cell model system")
            if values.get("developmental_stage") is None:
                raise ValueError(
                    "Developmental stage is required for primary cell model system"
                )

        if model_system == "Whole organism":
            if values.get("cell_line") is not None:
                raise ValueError(
                    "Cell line name is not required for whole organism model system"
                )
            if values.get("tissue") is not None:
                raise ValueError(
                    "Tissue is not required for whole organism model system"
                )
            if values.get("cell_type") is not None:
                raise ValueError(
                    "Cell type is not required for whole organism model system"
                )
            if values.get("passage_number") is not None:
                raise ValueError(
                    "Passage number is not required for whole organism model system"
                )
            if values.get("developmental_stage") is None:
                raise ValueError(
                    "Developmental stage is required for whole organism model system"
                )

        return values


class Disease(BaseModel):
    term_id: str = Field(
        ...,
        description="Term ID in CURIE format of the phenotype defined in EFO under parent term EFO:0000408 (disease)",
        example="MONDO:0004975",
    )
    term_label: str = Field(
        ...,
        description="Label of the phenotype defined in EFO under parent term EFO:0000408 (disease)",
        example="Alzheimer disease",
    )


class AssociatedDatasets(BaseModel):
    dataset_accession: str = Field(
        ...,
        description="Accession number of the dataset (e.g. GEO, ArrayExpress)",
        example="GSE123456",
    )
    dataset_uri: Optional[HttpUrl] = None
    dataset_description: str = Field(
        ...,
        description="Short description of the dataset",
        example="RNA-seq data from a CRISPR screen",
    )
    dataset_file_name: str = Field(
        ..., description="File name of the dataset", example="GSE123456_counts.txt"
    )


# Assemble the main Experiment model
class Experiment(BaseModel):
    study: StudyDetails
    experiment: ExperimentDetails
    perturbation: PerturbationDetails
    assay: AssayDetails
    model_system: ModelSystemDetails
    associated_diseases: List[Disease]
    associated_datasets: List[AssociatedDatasets]


# In-memory database for testing
pertcat_db: list[Experiment] = []

# Save the schema to a file
# schema = Experiment.model_json_schema()
# with open("unified_metadata_schema.json", "w") as f:
#     json.dump(schema, f, indent=4)


# FastAPI instance
app = FastAPI()


@app.post("/metadata/", response_model=Experiment)
def add_metadata(metadata: Experiment):
    if metadata in pertcat_db:
        raise HTTPException(status_code=400, detail="Study metadata already exists")
    pertcat_db.append(metadata)
    return metadata


@app.get("/metadata/", response_model=List[Experiment])
def read_metadata():
    return pertcat_db


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="127.0.0.1", port=8000)
