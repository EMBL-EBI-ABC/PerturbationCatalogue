import pandas as pd
import pandera as pa
from pandera import Field, DataFrameModel
from pandera.typing import Series, Index
from typing import Optional
from pathlib import Path

# Get the absolute path to the current module file
module_path = Path(__file__).resolve()

# Navigate up to the project root and then to 'ontologies'
ont_dir = module_path.parent.parent.parent / "ontologies"

gene_ont = pd.read_parquet(ont_dir / "genes.parquet").drop_duplicates()
ctype_ont = pd.read_parquet(ont_dir / "cell_types.parquet").drop_duplicates()
cline_ont = pd.read_parquet(ont_dir / "cell_lines.parquet").drop_duplicates()
tis_ont = pd.read_parquet(ont_dir / "tissues.parquet").drop_duplicates()
dis_ont = pd.read_parquet(ont_dir / "diseases.parquet").drop_duplicates()


# adata.obs schema
class ObsSchema(DataFrameModel):
    perturbation_name: Series[str] = Field(nullable=False)
    perturbed_target_coord: Series[str] = Field(nullable=False)
    perturbed_target_number: Series[int] = Field(nullable=False, ge=0)
    perturbed_target_ensg: Series[str] = Field(nullable=True)
    perturbed_target_symbol: Optional[Series[str]] = Field(nullable=True)
    perturbed_target_biotype: Optional[Series[str]] = Field(nullable=True)
    perturbation_type_label: Series[str] = Field(nullable=False)
    perturbation_type_id: Series[str] = Field(nullable=True, str_contains=":")
    timepoint: Optional[Series[str]] = Field(
        nullable=True, regex=r"^P\d+DT\d{1,2}H\d{1,2}M\d{1,2}S$"
    )
    treatment_label: Optional[Series[str]] = Field(nullable=True)
    treatment_id: Optional[Series[str]] = Field(nullable=True, str_contains=":")
    model_system_label: Series[str] = Field(nullable=False)
    model_system_id: Series[str] = Field(nullable=True, str_contains=":")
    species: Series[str] = Field(nullable=False, isin=["Homo sapiens"])
    tissue_label: Optional[Series[str]] = Field(nullable=True)
    tissue_id: Optional[Series[str]] = Field(nullable=True)
    cell_type_label: Series[str] = Field(nullable=False)
    cell_type_id: Series[str] = Field(nullable=False)
    cell_line_label: Optional[Series[str]] = Field(nullable=True)
    cell_line_id: Optional[Series[str]] = Field(nullable=True)
    sex_label: Optional[Series[str]] = Field(nullable=True)
    sex_id: Optional[Series[str]] = Field(nullable=True, str_contains=":")
    developmental_stage_label: Optional[Series[str]] = Field(nullable=True)
    developmental_stage_id: Optional[Series[str]] = Field(
        nullable=True, str_contains=":"
    )
    disease_label: Optional[Series[str]] = Field(nullable=True)
    disease_id: Optional[Series[str]] = Field(nullable=True)

    class Config:
        strict = True
        coerce = True
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
