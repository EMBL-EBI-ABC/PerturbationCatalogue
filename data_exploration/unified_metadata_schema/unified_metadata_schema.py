import json
import os
from pydantic import (
    BaseModel,
    model_validator,
    field_serializer,
    HttpUrl,
    Field,
)
from pydantic_core import Url
from typing import Optional, List, Literal
from fastapi import FastAPI, HTTPException
from datetime import datetime


# Define the models

# TermOptional and TermRequired are used to define terms that are linked with ontology
# TermRequired is used when the term is already present in the ontology, hence both term_id and term_label can be defined
class TermOptional(BaseModel):
    term_id: Optional[str] = Field(
        None,
        description="Ontology term ID in CURIE format",
        pattern=r"^[a-zA-Z0-9_]+:[0-9]+$"
    )
    term_label: str = Field(
        ...,
        description="Ontology term label"
    )


# TermOptional is used when the term is not yet present in the ontology, hence only term_label can be defined
# Once the necessary terms are added to the ontology, the schema should be updated to use TermRequired instead
class TermRequired(BaseModel):
    term_id: str = Field(
        ...,
        description="Ontology term ID in CURIE format",
        pattern=r"^[a-zA-Z0-9_]+:[0-9]+$"
    )
    term_label: str = Field(
        ...,
        description="Ontology term label"
    )


# Author first and last name
class Author(BaseModel):
    first_name: str = Field(..., description="First name of the author", example="John")
    last_name: str = Field(..., description="Last name of the author", example="Doe")


class StudyDetails(BaseModel):
    title: str = Field(..., description="Title of the study/publication")
    study_uri: Optional[HttpUrl|str] = Field(
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
    
    @field_serializer('study_uri')
    def convert_uri_to_string(self, val):
        if isinstance(val, Url):
            return str(val)
        return val


class SampleQuantity(BaseModel):
    sample_quantity_value: float = Field(
        ..., description="Sample quantity value", example=1.0
    )
    sample_quantity_unit: Literal[
        "gram", "liter", "unit", "colony-forming unit", "cells"
    ] = Field(..., description="Sample quantity unit", example="gram")


class ExperimentDetails(BaseModel):
    title: str = Field(..., description="Title of the experiment")
    summary: str = Field(..., description="Short summary of the experiment")
    treatments: Optional[List[TermRequired]] = Field(
        None,
        description="List of treatments used in the experiment, defined in ChEBI under 'chemical entity' CHEBI:24431 parent",
        example="CHEBI:28262 - dimethyl sulfoxide",
    )
    timepoints: Optional[List[str]] = Field(
        None,
        description="List of timepoints captured in the experiment. Must be in the ISO 8601 format",
        example="P1DT6H30M0S",
        # pattern=r"^P\d+DT\d{1,2}H\d{1,2}M\d{1,2}S$"
    )
    replicates: str = Field(
        ...,
        description="Types of replicates used the experiment: biological, technical or biological and technical, none",
        example="biological",
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
    perturbation_type: Optional[List[TermOptional]] = Field(
        ...,
        escription="Type of perturbation used in the experiment (e.g. CRISPRko, CRISPRi, CRISPRa, Mutagenesis)",
        example="CRISPRko",
    )
    perturbed_target_category: List[str] = Field(
        ...,
        description="Biotype of the perturbed target defined by ENSEMBL (e.g. protein coding, regulatory, etc.)",
        example="protein coding",
    )
    number_of_perturbed_targets: int = Field(
        ...,
        ge=1,
        description="How many targets (genes or variants) have been perturbed in the experiment",
        example=4,
    )
    perturbed_targets: List[str] = Field(
        ...,
        description="List of ENSEMBL IDs for perturbed targets (genes or variants) in the experiment",
        example=["ENSG00000141510", "ENSG00000146648"],
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
    library_format: TermOptional = Field(
        ...,
        description="Format of the library used in the experiment (pooled or arrayed)",
        example="pooled",
    )
    library_scope: TermOptional = Field(
        ...,
        description="Scope of the library used in the experiment (genome-wide or focused)",
        example="genome-wide",
    )
    library_perturbation_type: Optional[List[TermOptional]] = Field(
        ...,
        description="Type of perturbation used in the experiment (e.g. knockout, inhibition, activation, base editing, prime editing)",
        example="knockout",
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


    # If library_perturbation_type is "saturation mutagenesis", then total_variants is required
    @model_validator(mode="before")
    @classmethod
    def validate_total_variants(cls, values):
        
        library_perturbation_type = values.get("library_perturbation_type")[0]['term_label']
        total_variants = values.get("total_variants")
        library_scope = values.get("library_scope")['term_label']
        
        if (
            library_perturbation_type == "saturation mutagenesis"
            and total_variants is None
        ):
            raise ValueError("Total variants is required for saturation mutagenesis")
        if (
            library_perturbation_type != "saturation mutagenesis"
            and total_variants is not None
        ):
            raise ValueError(
                "Total variants is not required for this perturbation type"
            )
        if (
            library_perturbation_type == "saturation mutagenesis"
            and library_scope == "genome-wide"
        ):
            raise ValueError(
                "Saturation mutagenesis perturbation type is not allowed for Genome-wide libraries"
            )
        return values


class PerturbationDetails(BaseModel):
    library_generation_type: TermRequired = Field(
        ...,
        description="Term ID in CURIE format and term label of the library generation type defined in EFO under parent term EFO:0022867 (genetic perturbation)",
        example="EFO:0022868 - endogenous",
    )
    library_generation_method: TermRequired = Field(
        ...,
        description="Library generation method used in the experiment",
        example="EFO:0022895 - dCas9-KRAB",
    )
    enzyme_delivery_method: TermOptional = Field(
        None,
        description="Enzyme delivery method used in the experiment",
        example="retroviral transduction",
    )
    library_delivery_method: TermOptional = Field(
        ...,
        description="Library delivery method used in the experiment",
        example="lentiviral transduction",
    )
    enzyme_integration_state: TermOptional = Field(
        None,
        description="Integration status of the enzyme used in the experiment",
        example="random locus integration",
    )
    library_integration_state: TermOptional = Field(
        ...,
        description="Integration status of the library sequences used in the experiment",
        example="episomal expression",
    )
    enzyme_expression_control: TermOptional = Field(
        None,
        description="Expression control of the enzyme used in the experiment",
        example="constitutive expression",
    )
    library_expression_control: TermOptional = Field(
        ...,
        description="Expression control of the library used in the experiment",
        example="inducible expression",
    )
    library: Library

    @model_validator(mode="before")
    @classmethod
    def validate_enzyme_delivery_method(cls, values):
        
        enzyme_delivery_method = values.get("enzyme_delivery_method")['term_label']
        library_generation_type = values.get("library_generation_type")['term_label']
        
        if library_generation_type == "exogenous" and enzyme_delivery_method != None:
            raise ValueError(
                "Enzyme delivery method is not required for exogenous library generation type"
            )
        if library_generation_type == "endogenous" and enzyme_delivery_method == None:
            raise ValueError(
                "Enzyme delivery method is required for endogenous library generation type"
            )
        return values

    @model_validator(mode="before")
    @classmethod
    def validate_enzyme_expression_control(cls, values):
        
        enzyme_expression_control = values.get("enzyme_expression_control")['term_label']
        library_generation_type = values.get("library_generation_type")['term_label']
        
        if library_generation_type == "exogenous" and enzyme_expression_control != None:
            raise ValueError(
                "Enzyme expression control is not required for exogenous library generation type"
            )
        if (
            library_generation_type == "endogenous"
            and enzyme_expression_control == None
        ):
            raise ValueError(
                "Enzyme expression control is required for endogenous library generation type"
            )
        return values


class AssayDetails(BaseModel):
    readout_dimensionality: TermOptional = Field(
        ...,
        description="Dimensionality of the readout (i.e single-dimensional assay, high-dimensional assay)",
        example="single-dimensional assay",
    )
    readout_type: TermOptional = Field(
        ...,
        description="Type of the readout (e.g. transcriptomic, phenotypic)",
        example="transcriptomic",
    )
    readout_technology: TermOptional = Field(
        ...,
        description="Technology used for the readout (e.g. population growth assay, flow cytometry)",
        example="population growth assay",
    )
    method_name: TermOptional = Field(
        ...,
        description="Name of the method used for the readout (e.g. scRNA-seq, Proliferation CRISPR screen)",
        example="Proliferation CRISPR screen",
    )
    method_uri: Optional[HttpUrl|str] = None
    sequencing_library_kit: TermOptional = Field(
        ...,
        description="Sequencing library kit used in the experiment (e.g. 10x Genomics Single Cell 3-prime v1)",
        example="10x Genomics Single Cell 3-prime v1",
    )
    sequencing_platform: TermOptional = Field(
        ...,
        description="Sequencing platform used in the experiment (e.g. Illumina NovaSeq 6000, ONT MinION)",
        example="Illumina NovaSeq 6000",
    )
    sequencing_strategy: TermOptional = Field(
        ...,
        description="Sequencing strategy used in the experiment (e.g. shotgun sequencing, barcode sequencing)",
        example="shotgun sequencing",
    )
    software_counts: TermOptional = Field(
        ...,
        description="Software used for counting the reads (e.g. CellRanger, MAGeCK)",
        example="CellRanger",
    )
    software_analysis: TermOptional = Field(
        ...,
        description="Software used for analyzing the data (e.g. Seurat, MAGeCK)",
        example="Seurat",
    )
    reference_genome: TermOptional = Field(
        ...,
        description="Reference genome used in the experiment (e.g. GRCh38, GRCh37, T2T-CHM13)",
        example="GRCh38",
    )
    
    @field_serializer('method_uri')
    def convert_uri_to_string(self, val):
        if isinstance(val, Url):
            return str(val)
        return val


class ModelSystemDetails(BaseModel):
    model_system: Optional[List[TermOptional]] = Field(
        ...,
        description="Model system used in the experiment. Must be a term ID in CURIE format from EFO 'disease' EFO:0000408 parent",
    )
    species: Literal["Homo sapiens"] = Field(
        ...,
        description="Species used in the experiment; must be 'Homo sapiens'",
        example="Homo sapiens",
    )
    tissue: Optional[List[TermRequired]] = Field(
        None,
        description="Specific biological tissue the sample is derived from. Must be a term ID in CURIE format from UBERON 'anatomical entity' UBERON:0001062 parent",
        example="UBERON:0002098",
    )
    cell_type: Optional[List[TermRequired]] = Field(
        None,
        description="Cell type/types profiled in the experiment. Must be a term ID in CURIE format from Cell Ontology 'cell' CL:0000000 parent",
        example="CL:0000235",
    )
    cell_line: Optional[List[TermRequired]] = Field(
        None,
        description="Cell line name used in the experiment. Must be a term ID in CURIE format from Cell Line Ontology 'cultured cell' CL:0000010 parent",
        example="CLO:0009348",
    )
    sex: Optional[List[TermOptional]] = Field(
        None, description="Model system organism's sex", example="female"
    )
    developmental_stage: Optional[List[TermOptional]] = Field(
        None, description="Developmental stage/age of the model system", example="adult"
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


class AssociatedDatasets(BaseModel):
    dataset_accession: str = Field(
        ...,
        description="Accession number of the dataset (e.g. GEO, ArrayExpress)",
        example="GSE123456",
    )
    dataset_uri: Optional[HttpUrl|str] = None
    dataset_description: str = Field(
        ...,
        description="Short description of the dataset",
        example="RNA-seq data from a CRISPR screen",
    )
    dataset_file_name: str = Field(
        ..., description="File name of the dataset", example="GSE123456_counts.txt"
    )
    
    @field_serializer('dataset_uri')
    def convert_uri_to_string(self, val):
        if isinstance(val, Url):
            return str(val)
        return val


# Assemble the main Experiment model
class Experiment(BaseModel):
    study: StudyDetails = Field(..., description="Details of the study/publication")
    experiment: ExperimentDetails = Field(..., description="Details of the experiment")
    perturbation: PerturbationDetails = Field(
        ..., description="Details of the perturbation"
    )
    assay: AssayDetails = Field(..., description="Details of the assay")
    model_system: ModelSystemDetails = Field(
        ..., description="Details of the model system"
    )
    associated_diseases: Optional[List[TermRequired]] = Field(
        None,
        description="Term ID in CURIE format of the phenotype defined in EFO under parent term EFO:0000408 (disease)",
        examples=["MONDO:000497 - Alzheimer disease"],
    )
    associated_datasets: List[AssociatedDatasets] = Field(
        ...,
        description="List of datasets associated with the experiment",
    )


# In-memory database for testing
pertcat_db: list[Experiment] = []


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
