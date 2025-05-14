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
    perturbed_target_number: Series[int] = Field(nullable=False, ge=1)
    perturbed_target_ensg: Series[str] = Field(nullable=False)
    perturbed_target_symbol: Optional[Series[str]] = Field(nullable=True)
    perturbed_target_category: Optional[Series[str]] = Field(nullable=True)
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
    tissue_label: Optional[Series[str]] = Field(nullable=True, isin=tis_ont.name.values)
    tissue_id: Optional[Series[str]] = Field(
        nullable=True,
        isin=tis_ont.ontology_id.values,
    )
    cell_type_label: Series[str] = Field(nullable=False, isin=ctype_ont.name.values)
    cell_type_id: Series[str] = Field(
        nullable=False,
        isin=ctype_ont.ontology_id.values,
    )
    cell_line_label: Optional[Series[str]] = Field(
        nullable=True, isin=cline_ont.name.values
    )
    cell_line_id: Optional[Series[str]] = Field(
        nullable=True,
        isin=cline_ont.ontology_id.values,
    )
    sex_label: Optional[Series[str]] = Field(nullable=True)
    sex_id: Optional[Series[str]] = Field(nullable=True, str_contains=":")
    developmental_stage_label: Optional[Series[str]] = Field(nullable=True)
    developmental_stage_id: Optional[Series[str]] = Field(
        nullable=True, str_contains=":"
    )
    disease_term_label: Optional[Series[str]] = Field(
        nullable=True, isin=dis_ont.name.values
    )
    disease_term_id: Optional[Series[str]] = Field(
        nullable=True,
        isin=dis_ont.ontology_id.values,
    )

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
        # isin=gene_ont.ensembl_gene_id.values
    )
    gene_symbol: Series[str] = Field(
        nullable=True,
        coerce=True,
        # isin=gene_ont.symbol.values
    )
    original_gene_symbol: Optional[Series[str]] = Field(nullable=False)
    original_ensembl_gene_id: Optional[Series[str]] = Field(nullable=False)

    @pa.dataframe_check
    def validate_gene_identifier_columns(cls, df):
        """
        Validates the presence of gene identifier columns in the DataFrame.
        Exactly one of 'original_ensembl_gene_id' or 'original_gene_symbol' must be present.
        """
        columns = set(df.columns)
        has_ensembl = "original_ensembl_gene_id" in columns
        has_symbol = "original_gene_symbol" in columns
        return (has_ensembl != has_symbol)  # True if exactly one is present

    class Config:
        strict = "filter"
        coerce = True
        ordered = True
