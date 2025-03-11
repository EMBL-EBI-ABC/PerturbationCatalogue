import json
from pydantic import BaseModel, model_validator, HttpUrl, Field
from typing import Optional, List
from fastapi import FastAPI, HTTPException
from enum import Enum

# function to make enums
def make_enum(name, values):
    return Enum(name, {i.replace(" ", "_").lower(): i for i in values})

# load the options from the JSON file
with open("enums.json", "r") as f:
    options = json.load(f)

# create the enums
TimepointUnit = make_enum("TimepointUnit", options["timepoint_unit"])
Replicates = make_enum("Replicates", options["replicates"])
LibraryGenerationType = make_enum(
    "LibraryGenerationType", options["library_generation_type"]
)
LibraryGenerationMethod = make_enum(
    "LibraryGenerationMethod", options["library_generation_method"]
)
DeliveryMethod = make_enum("DeliveryMethod", options["delivery_method"])
IntegrationState = make_enum("IntegrationState", options["integration_state"])
ExpressionControl = make_enum("ExpressionControl", options["expression_control"])
LibraryName = make_enum("LibraryName", ["Other"] + list(options["libraries"]))
LibraryScope = make_enum("LibraryScope", options["library_scope"])
LibraryFormat = make_enum("LibraryFormat", options["library_format"])
PerturbationType = make_enum("PerturbationType", options["perturbation_type"])
ReadoutDimensionality = make_enum(
    "ReadoutDimensionality", options["readout_dimensionality"]
)
ReadoutType = make_enum("ReadoutType", options["readout_type"])
ReadoutTechnology = make_enum("ReadoutTechnology", options["readout_technology"])
MethodName = make_enum("MethodName", options["method_name"])
SequencingLibraryKits = make_enum("SequencingLibraryKits", options["sequencing_library_kits"])
SequencingPlatform = make_enum("SequencingPlatform", options["sequencing_platform"])
SequencingStrategy = make_enum("SequencingStrategy", options["sequencing_strategy"])
SoftwareCounts = make_enum("SoftwareCounts", options["software_counts"])
SoftwareAnalysis = make_enum("SoftwareAnalysis", options["software_analysis"])
ReferenceGenome = make_enum("ReferenceGenome", options["reference_genome"])
ModelSystem = make_enum("ModelSystem", options["model_system"])
Species = make_enum("Species", options["species"])
Tissue = make_enum("Tissue", options["tissue"])
CellType = make_enum("CellType", options["cell_type"])
ModelName = make_enum("ModelName", options["model_name"])
Sex = make_enum("Sex", options["sex"])
DevelopmentalStages = make_enum("DevelopmentalStages", options["developmental_stages"])
SampleQuantityUnit = make_enum("SampleQuantityUnit", options["sample_quantity_unit"])

# Define the models
class Author(BaseModel):
    first_name: str
    last_name: str

class StudyDetails(BaseModel):
    title: str
    uri: Optional[HttpUrl] = None
    year: int = Field(ge=1000, le=3000)
    first_author: Author
    last_author: Author

class Timepoint(BaseModel):
    timepoint_value: float
    timepoint_unit: TimepointUnit

class Treatment(BaseModel):
    term_id: str = Field(..., description="Term ID of the treatment defined in ChEMBL", example="CHEBI:6904")
    term_label: str = Field(..., description="Label of the treatment defined in ChEMBL", example="metoprolol")

class ExperimentDetails(BaseModel):
    title: str
    summary: str
    treatments: Optional[List[Treatment]] = None
    timepoints: Optional[List[Timepoint]] = None
    replicates: Replicates
    number_of_samples: int = Field(ge=1)
    number_of_perturbed_cells: int = Field(ge=1)
    number_of_perturbed_genes: int = Field(ge=1)
    coverage: float = Field(ge=0)

class Library(BaseModel):
    library_name: LibraryName
    accession: Optional[str] = None
    library_format: LibraryFormat
    library_scope: LibraryScope
    perturbation_type: PerturbationType
    manufacturer: Optional[str] = None
    lentiviral_generation: Optional[str] = None
    grnas_per_gene: Optional[str] = None
    total_grnas: Optional[int] = Field(None, ge=1)
    total_genes: Optional[int] = Field(None, ge=1)
    total_variants: Optional[int] = Field(None, ge=1)

    # Validate the Enum library. If the library name is not "Other" (i.e. one of the Enum library values),
    # update the library parameters with the pre-defined values from the "options" dictionary
    @model_validator(mode="before")
    @classmethod
    def validate_enum_library(cls, values):
        lib_name = values["library_name"]
        if lib_name != "Other":
            lib = options["libraries"][lib_name]
            values.update(lib)
            print("Library parameters have been updated")
        return values

class PerturbationDetails(BaseModel):
    library_generation_type: LibraryGenerationType
    library_generation_method: LibraryGenerationMethod
    enzyme_delivery_method: DeliveryMethod
    library_delivery_method: DeliveryMethod
    enzyme_integration_state: IntegrationState
    library_integration_state: IntegrationState
    enzyme_expression_control: ExpressionControl
    library_expression_control: ExpressionControl
    library: Library

class AssayDetails(BaseModel):
    readout_dimensionality: ReadoutDimensionality
    readout_type: ReadoutType
    readout_technology: ReadoutTechnology
    method_name: MethodName
    method_uri: Optional[HttpUrl] = None
    sequencing_library_kit: SequencingLibraryKits
    sequencing_platform: SequencingPlatform
    sequencing_strategy: SequencingStrategy
    software_counts: SoftwareCounts
    software_analysis: SoftwareAnalysis
    reference_genome: ReferenceGenome

class ModelSystemDetails(BaseModel):
    model_system: ModelSystem
    species: Species
    tissue: Optional[Tissue] = None
    cell_type: Optional[CellType] = None
    model_name: Optional[ModelName] = None
    sex: Sex
    developmental_stage: DevelopmentalStages
    passage_number: Optional[int] = Field(None, ge=1)
    sample_quantity: float = Field(gt=0)
    sample_quantity_unit: SampleQuantityUnit
    
    @model_validator(mode="before")
    @classmethod
    def validate_enum_model_system(cls, values):
        model_system = values["model_system"]
        
        if model_system == "Tissue sample":
            if values.get("tissue") is None:
                raise ValueError("Tissue is required for tissue sample model system")
            if values.get("cell_type") is not None:
                raise ValueError("Cell type is not required for tissue sample model system")
            if values.get("model_name") is not None:
                raise ValueError("Model name is not required for tissue sample model system")
            if values.get("passage_number") is not None:
                raise ValueError("Passage number is not required for tissue sample model system")
            
        if model_system == "Cell line":
            if values.get("cell_type") is None:
                raise ValueError("Cell type is required for cell line model system")
            if values.get("passage_number") is None:
                raise ValueError("Passage number is required for cell line model system")
            if values.get("model_name") is None:
                raise ValueError("Model name is required for cell line model system")
            if values.get("tissue") is not None:
                raise ValueError("Tissue is not required for cell line model system")
            
        if model_system == "Primary cell":
            if values.get("cell_type") is None:
                raise ValueError("Cell type is required for primary cell model system")
            if values.get("passage_number") is None:
                raise ValueError("Passage number is required for primary cell model system")
            if values.get("model_name") is not None:
                raise ValueError("Model name is not required for primary cell model system")
            if values.get("tissue") is not None:
                raise ValueError("Tissue is not required for primary cell model system")
        
        if model_system == "Whole organism":
            if values.get("model_name") is not None:
                raise ValueError("Model name is not required for whole organism model system")
            if values.get("tissue") is not None:
                raise ValueError("Tissue is not required for whole organism model system")
            if values.get("cell_type") is not None:
                raise ValueError("Cell type is not required for whole organism model system")
            if values.get("passage_number") is not None:
                raise ValueError("Passage number is not required for whole organism model system")
            
        return values
            
class Phenotype(BaseModel):
    term_id: str = Field(..., description="Term ID of the phenotype defined in EFO under parent term EFO:0000408 (disease)", example="MONDO:0004975")
    term_label: str = Field(..., description="Label of the phenotype defined in EFO under parent term EFO:0000408 (disease)", example="Alzheimer disease")

class AssociatedDatasets(BaseModel):
    dataset_accession: str
    dataset_uri: Optional[HttpUrl] = None
    dataset_description: str
    dataset_file_name: str

# Assemble the main Experiment model
class Experiment(BaseModel):
    study: StudyDetails
    experiment: ExperimentDetails
    perturbation: PerturbationDetails
    assay: AssayDetails
    model_system: ModelSystemDetails
    associated_phenotypes: List[Phenotype]
    associated_datasets: List[AssociatedDatasets]


# In-memory database for testing
pertcat_db: list[Experiment] = []

# save the schema to a file
schema = Experiment.model_json_schema()
with open("unified_metadata_schema.json", "w") as f:
    json.dump(schema, f, indent=4)


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
