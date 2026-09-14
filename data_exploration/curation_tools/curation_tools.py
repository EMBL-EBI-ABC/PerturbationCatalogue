import datetime
import glob
import os
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import polars as pl
import pyarrow.parquet as pq
import json
from pprint import pprint
import requests
import logging

from pydantic import ValidationError
from typing import Literal
import pandera.pandas as pa
from pandera.typing import Series, Int64, String
from tqdm import tqdm
from thefuzz import process

from libchebipy import search
import scanpy as sc
import anndata as ad

import ibis
import re
from google.cloud import bigquery
from IPython.display import display  # type: ignore

from curation_tools.perturbseq_anndata_schema import ObsSchema, VarSchema

# Module-level logger
logger = logging.getLogger(__name__)


def _coerce_metadata_for_polars(
    metadata_df: pd.DataFrame, polars_schema: dict
) -> pd.DataFrame:
    """Normalize metadata columns before applying the Polars schema.

    Metadata can arrive from CSV or Excel files as strings, including integer
    values formatted as decimal strings such as ``"20724.0"``. Polars cannot
    cast those strings directly to ``Int64``, so numeric schema columns must be
    normalized in pandas first.
    """

    normalized_df = metadata_df.copy()

    for column_name, polars_dtype in polars_schema.items():
        if column_name not in normalized_df.columns:
            continue

        series = normalized_df[column_name].replace("None", pd.NA)

        if polars_dtype == pl.Int64:
            numeric = pd.to_numeric(series, errors="coerce")
            invalid = series.notna() & (numeric.isna() | ~np.isfinite(numeric))
            if invalid.any():
                invalid_values = series[invalid].head().tolist()
                raise ValueError(
                    f"Column {column_name!r} contains non-integer values: "
                    f"{invalid_values}"
                )

            non_integer = numeric.notna() & ((numeric % 1) != 0)
            if non_integer.any():
                invalid_values = series[non_integer].head().tolist()
                raise ValueError(
                    f"Column {column_name!r} contains non-integer values: "
                    f"{invalid_values}"
                )

            normalized_df[column_name] = numeric.astype("Int64")
        elif polars_dtype == pl.Float64:
            numeric = pd.to_numeric(series, errors="coerce")
            invalid = series.notna() & numeric.isna()
            if invalid.any():
                invalid_values = series[invalid].head().tolist()
                raise ValueError(
                    f"Column {column_name!r} contains non-numeric values: "
                    f"{invalid_values}"
                )

            normalized_df[column_name] = numeric.astype("Float64")
        elif polars_dtype == pl.Boolean:
            boolean = (
                series.astype("string")
                .str.strip()
                .str.lower()
                .map({"true": True, "false": False})
            )
            invalid = series.notna() & boolean.isna()
            if invalid.any():
                invalid_values = series[invalid].head().tolist()
                raise ValueError(
                    f"Column {column_name!r} contains non-boolean values: "
                    f"{invalid_values}"
                )

            normalized_df[column_name] = boolean.astype("boolean")
        elif polars_dtype == pl.String:
            normalized_df[column_name] = series.astype("string")

    return normalized_df


# function to add a new synonym to the ontology
def add_synonym(
    ontology_type=Literal["genes", "cell_types", "cell_lines", "tissues", "diseases"],
    ref_column=str,
    syn_column=str,
    syn_map=dict,
    save=True,
):
    """
    Add a new synonym to the specified ontology term.

    Parameters
    ----------
    ontology_type : str
        The name of the ontology type (e.g.,"genes", "cell_types", "cell_lines", "tissues", "diseases").
    ref_column : str
        The name of the column in the ontology DataFrame to use for matching terms.
    syn_column : str
        The name of the column in the ontology DataFrame to add synonyms to.
    syn_map : dict
        A dictionary mapping the ontology term to the new synonyms.
    save : bool
        Whether to save the updated ontology DataFrame to a parquet file. Default is True.
    """

    if ontology_type not in [
        "genes",
        "cell_types",
        "cell_lines",
        "tissues",
        "diseases",
    ]:
        raise ValueError(
            "ontology_type must be one of 'genes', 'cell_types', 'cell_lines', 'tissues', 'diseases'"
        )

    # Get the path to the ontologies directory relative to this file
    ONTOLOGIES_DIR = Path(__file__).parent / "ontologies"

    ont = pd.read_parquet(ONTOLOGIES_DIR / f"{ontology_type}.parquet").drop_duplicates()

    if ref_column not in ont.columns:
        raise ValueError(
            f"Column `{ref_column}` not found in `{ontology_type}` ontology"
        )
    if syn_column not in ont.columns:
        raise ValueError(
            f"Column `{syn_column}` not found in `{ontology_type}` ontology"
        )

    # map the synonyms to the ontology terms
    for term, synonyms in syn_map.items():
        if term not in ont[ref_column].values:
            raise ValueError(f"Term `{term}` not found in `{ontology_type}` ontology")
        if not isinstance(synonyms, list):
            raise ValueError("`syn_map` values must be a list of synonyms")

        # add the synonyms to the ontology
        for synonym in synonyms:
            if synonym not in ont[syn_column].values:
                ont.loc[ont[ref_column] == term, syn_column] += f"|{synonym}"

    # Display the updated terms
    print(f"Updated terms in `{ontology_type}` ontology:")
    display(ont.loc[ont[ref_column].isin(syn_map.keys()), :])

    # Save the updated ontology
    if save:
        ont.to_parquet(ONTOLOGIES_DIR / f"{ontology_type}.parquet", index=False)


class CuratedDataset:

    # Get the path to the ontologies directory relative to this file
    ONTOLOGIES_DIR = Path(__file__).parent / "ontologies"

    opentargets_gene_reference_path = (
        ONTOLOGIES_DIR / "opentargets_gene_identifier_reference_26_03.parquet"
    )
    ctype_ont = pd.read_parquet(ONTOLOGIES_DIR / "cell_types.parquet").drop_duplicates()
    cline_ont = pd.read_parquet(ONTOLOGIES_DIR / "cell_lines.parquet").drop_duplicates()
    tis_ont = pd.read_parquet(ONTOLOGIES_DIR / "tissues.parquet").drop_duplicates()
    dis_ont = pd.read_parquet(ONTOLOGIES_DIR / "diseases.parquet").drop_duplicates()

    def __init__(
        self,
        obs_schema,
        var_schema,
        data_source_link=None,
        noncurated_path=None,
        curated_path=None,
    ):
        """Initialize the CuratedDataset class.
        Parameters
        ----------
        obs_schema : ObsSchema
            The schema for the obs data.
        var_schema : VarSchema
            The schema for the var data.
        data_source_link : str, optional
            The link to the data source. The default is None.
        noncurated_path : str, optional
            The path to the non-curated .h5ad data. The default is None.
        curated_path : str, optional
            The path to the curated .h5ad data. The default is None.
        """

        self.obs_schema = obs_schema
        self.var_schema = var_schema

        self.data_source_link = data_source_link
        self.noncurated_path = noncurated_path
        self.curated_path = curated_path

        if self.noncurated_path and not self.curated_path:
            self.curated_path = noncurated_path.replace(
                "non_curated", "curated"
            ).replace(".h5ad", "_curated.h5ad")
        elif self.curated_path:
            self.noncurated_path = curated_path.replace(
                "curated", "non_curated"
            ).replace("_curated.h5ad", ".h5ad")

        self.curated_parquet_data_path = self.curated_path.replace(
            ".h5ad", "_data.parquet"
        ).replace("h5ad", "parquet")
        self.curated_parquet_metadata_path = self.curated_path.replace(
            ".h5ad", "_metadata.parquet"
        ).replace("h5ad", "parquet")

        # Initialise adata object
        self.adata = None

        # Initialise a dataset ID
        if self.noncurated_path:
            self.dataset_id = self.noncurated_path.split("/")[-1].replace(".h5ad", "")
        elif self.curated_path:
            self.dataset_id = self.curated_path.split("/")[-1].replace(
                "_curated.h5ad", ""
            )
        else:
            raise ValueError(
                "Either noncurated_path or curated_path must be provided to initialize the dataset ID."
            )

    def show_obs(self, obs_columns=None):
        """
        Display the observation data.
        Parameters
        ----------
        obs_columns : list, optional
            A list of columns to display. The default is None, which displays all columns.
        """
        if obs_columns is None:
            obs_columns = self.adata.obs.columns
        elif not isinstance(obs_columns, list):
            obs_columns = [obs_columns]
        print("Observation data:")
        self.print_data(self.adata.obs[obs_columns])

    def show_var(self, var_columns=None):
        """
        Display the variable data.
        Parameters
        ----------
        var_columns : list, optional
            A list of columns to display. The default is None, which displays all columns.
        """
        if var_columns is None:
            var_columns = self.adata.var.columns
        elif not isinstance(var_columns, list):
            var_columns = [var_columns]
        print("Variable data:")
        self.print_data(self.adata.var[var_columns])

    def show_unique(self, slot=Literal["var", "obs"], column=None):
        """
        Display the unique values in a column of the specified slot of the adata object.
        Parameters
        ----------
        slot : str
            The slot to display unique values from. Can be either "obs" or "var".
        column : str, optional
            The name of the column to display unique values from. If None, displays all columns.
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')
        df = getattr(self.adata, slot)
        if column is None:
            raise ValueError("`column` must be specified")
        else:
            if column not in df.columns:
                raise ValueError(f"Column {column} not found in adata.{slot}")
            else:
                print(f"Unique values in adata.{slot}.{column}: {len(set(df[column]))}")
                self.print_data(set(df[column]))

    def download_data(self):
        """
        Download the data from the specified source.
        """
        if not os.path.exists(self.noncurated_path):
            print(
                f"Downloading data from {self.data_source_link} to {self.noncurated_path}"
            )
            os.makedirs(os.path.dirname(self.noncurated_path), exist_ok=True)
            os.system(f"wget {self.data_source_link} -O {self.noncurated_path}")
        else:
            print(f"File {self.noncurated_path} already exists. Skipping download.")

    def load_data(self, curated=False):
        """
        Load the adata from a curated or non-curated path.
        """
        if curated:
            if os.path.exists(self.curated_path):
                print(f"Loading data from {self.curated_path}")
                self.adata = sc.read_h5ad(self.curated_path)
            else:
                raise ValueError(
                    f"File {self.curated_path} does not exist. Check the path."
                )
        else:
            if os.path.exists(self.noncurated_path):
                print(f"Loading data from {self.noncurated_path}")
                self.adata = sc.read_h5ad(self.noncurated_path)
            else:
                raise ValueError(
                    f"File {self.noncurated_path} does not exist. Run download_data() first."
                )
            for col in self.adata.obs.columns:
                if self.adata.obs[col].dtype == "object":
                    self.adata.obs[col] = self.adata.obs[col].astype(
                        "str"
                    )  # .astype("category")

    def save_curated_data_h5ad(self):
        """Save the curated data to an .h5ad file."""

        ad.settings.allow_write_nullable_strings = (
            True  # Allow nullable strings in adata
        )

        adata = self.adata

        if adata is None:
            raise ValueError("adata is not loaded. Please load the data first.")

        # check if the base directory exists, if not, create it
        if not os.path.exists(os.path.dirname(self.curated_path)):
            os.makedirs(os.path.dirname(self.curated_path))
        # Replace None with np.nan in adata.obs and adata.var
        adata.obs = adata.obs.fillna(value=np.nan)
        adata.var = adata.var.fillna(value=np.nan)

        adata.write_h5ad(self.curated_path)
        print(f"✅ Curated h5ad data saved to {self.curated_path}")

    def polars_schema_from_pandera_model(self):
        """
        Extract a Polars schema dictionary from a Pandera Polars DataFrameModel.
        The dictionary maps column names to Polars data types.
        Uses isinstance() to check column dtype properly.
        """
        pandera_model = self.obs_schema
        schema_dict = {}
        for col_name, column in pandera_model.to_schema().columns.items():
            pandera_dtype = column.dtype

            if isinstance(pandera_dtype, String):
                polars_dtype = pl.String
            elif isinstance(pandera_dtype, Int64):
                polars_dtype = pl.Int64
            elif isinstance(pandera_dtype, float):
                polars_dtype = pl.Float64
            elif isinstance(pandera_dtype, bool):
                polars_dtype = pl.Boolean
            else:
                # Default to Utf8 for unknown or unhandled dtypes
                polars_dtype = pl.Utf8

            schema_dict[col_name] = polars_dtype

        return schema_dict

    def save_curated_data_parquet(
        self, split_metadata=False, save_metadata_only=False, overwrite=False
    ):
        """Save the curated data to a parquet file ready for BigQuery ingestion.

        Parameters
        ----------
        split_metadata : bool
            Whether to split the data and metadata into two separate files (default is False).
        save_metadata_only : bool
            Whether to save only the metadata and skip saving the data (default is False).
        overwrite : bool
            Whether to overwrite existing Parquet files. Defaults to False.
        """

        adata = self.adata

        if adata is None:
            raise ValueError("adata is not loaded. Please load the data first.")

        # check if the base directory exists, if not, create it
        if not os.path.exists(os.path.dirname(self.curated_path)):
            os.makedirs(os.path.dirname(self.curated_path))

        polars_schema = self.polars_schema_from_pandera_model()

        # Normalize metadata according to the schema without stringifying
        # numeric columns before converting to Polars.
        full_metadata_df = _coerce_metadata_for_polars(adata.obs, polars_schema)

        metadata_columns = full_metadata_df.columns.to_list()
        id_columns = metadata_columns[0:2]

        # Process features (e.g. genes or scores) in chunks
        feature_colnames = adata.var_names.tolist()

        if not split_metadata:
            parquet_path = self.curated_path.replace(
                ".h5ad", "_unified.parquet"
            ).replace("h5ad", "parquet")
            if os.path.exists(parquet_path) and not overwrite:
                raise FileExistsError(
                    f"File {parquet_path} already exists. Skipping write."
                )
            if not os.path.exists(os.path.dirname(parquet_path)):
                os.makedirs(os.path.dirname(parquet_path))
            else:
                X_df = adata.to_df()

                full_data_df = pd.concat(
                    [
                        full_metadata_df.reset_index(drop=True),
                        X_df.reset_index(drop=True),
                    ],
                    axis=1,
                    ignore_index=False,
                )
                full_data_df = pl.from_pandas(
                    full_data_df, schema_overrides=polars_schema
                )

                # convert to pyarrow for efficient streaming to parquet
                full_data_df = full_data_df.to_arrow()
                # TODO: if this works, try also using polars native parquet writer pl.write_parquet
                writer = pq.ParquetWriter(parquet_path, full_data_df.schema)
                print(
                    f"Created ParquetWriter and started writing full data to {parquet_path}"
                )
                writer.write_table(full_data_df)
                writer.close()
                print(f"✅ Unified data saved to {parquet_path}")

        # if split_metadata is True, save metadata and data separately
        else:
            if not os.path.exists(os.path.dirname(self.curated_parquet_data_path)):
                os.makedirs(os.path.dirname(self.curated_parquet_data_path))
            if not os.path.exists(os.path.dirname(self.curated_parquet_metadata_path)):
                os.makedirs(os.path.dirname(self.curated_parquet_metadata_path))

            if (
                not overwrite
                and (
                    os.path.exists(self.curated_parquet_data_path)
                    or os.path.exists(self.curated_parquet_metadata_path)
                )
            ):
                print(
                    f"Files {self.curated_parquet_data_path} or {self.curated_parquet_metadata_path} already exist. Skipping write."
                )
                return
            # if save_metadata_only is True, save only the metadata and skip saving the data
            if save_metadata_only:
                # Write metadata to parquet
                full_metadata_df.to_parquet(self.curated_parquet_metadata_path, index=False)
                print(f"✅ Metadata saved to {self.curated_parquet_metadata_path}")
                return

            else:
                # Write metadata to parquet
                full_metadata_df.to_parquet(self.curated_parquet_metadata_path, index=False)
                print(f"✅ Metadata saved to {self.curated_parquet_metadata_path}")

                print("Processing data...")
                X_df = adata.to_df()

                full_data_df = pd.concat(
                    [
                        full_metadata_df.reset_index(drop=True),
                        X_df.reset_index(drop=True),
                    ],
                    axis=1,
                    ignore_index=False,
                )
                full_data_df = pl.from_pandas(
                    full_data_df, schema_overrides=polars_schema
                )

                # Select only the ID and feature columns for the data file
                data_subset_df = full_data_df.select(id_columns + feature_colnames)

                data_subset_df = data_subset_df.unpivot(
                    on=feature_colnames,
                    index=id_columns,
                    variable_name="score_name",
                    value_name="score_value",
                )
                data_subset_df = data_subset_df.with_columns(
                    pl.col("score_value").cast(pl.Float64)
                )

                # save data_subset_df to parquet
                print(f"Saving data to {self.curated_parquet_data_path}...")
                data_subset_df.write_parquet(self.curated_parquet_data_path)
                print(f"✅ Data saved to {self.curated_parquet_data_path}")

    def upload_parquet_to_bq(
        self,
        project_id: str,
        bq_dataset_id: Literal["perturb_seq", "crispr", "mavedb"],
        bq_table_name: Literal["data", "metadata"],
        key_columns,
        parquet_path=None,
        verbose=True,
    ):
        """Upload a curated parquet file to BigQuery.

        Parameters
        ----------
        bq_table_name : Literal["data", "metadata"]
            The name of the BigQuery table to upload to: {project_id}.{bq_dataset_id}.{bq_table_name}
        key_columns : list[str]
            Columns used to merge staging rows into the destination table.
        parquet_path : str, optional
            Explicit parquet file path. If omitted, bq_table_name must be provided.
        bq_dataset_id : Literal["perturb_seq", "crispr", "mavedb"]
            BigQuery dataset ID.
        project_id : str, optional
            BigQuery project ID. Defaults to the BQ_PROJECT environment variable.
        verbose : bool
            Whether to print upload progress.
        """
        if parquet_path is None:
            if bq_table_name == "metadata":
                parquet_path = self.curated_parquet_metadata_path
            elif bq_table_name == "data":
                parquet_path = self.curated_parquet_data_path
            else:
                raise ValueError(
                    "parquet_path must be provided unless bq_table_name is "
                    "'metadata' or 'data'."
                )

        parquet_path = Path(parquet_path)
        if not parquet_path.exists():
            raise FileNotFoundError(f"Parquet file not found: {parquet_path}")

        # BigQuery MERGE needs stable keys to decide which rows to update or insert.
        if not key_columns:
            raise ValueError("key_columns must contain at least one column.")

        _upload_parquet_to_bq(
            parquet_path=parquet_path,
            project_id=project_id,
            bq_dataset_id=bq_dataset_id,
            bq_table_name=bq_table_name,
            key_columns=key_columns,
            verbose=verbose,
        )

    def chromosome_encoding(self, chromosome_col="perturbed_target_chromosome"):
        """
        Encode the chromosome column (default='perturbed_target_chromosome') in the adata.obs DataFrame.
        Parameters
        ----------
        chromosome_col : str
            The name of the column containing chromosome information.
        """
        if chromosome_col not in self.adata.obs.columns:
            raise ValueError(f"Column {chromosome_col} not found in adata.obs")

        # Create a mapping for chromosome encoding
        chromosome_encoding_dict = {
            **{str(i): i for i in range(1, 23)},
            "X": 23,
            "Y": 24,
            "MT": 25
        }

        # Apply the mapping to the chromosome column
        self.adata.obs["perturbed_target_chromosome_encoding"] = [
            chromosome_encoding_dict[x] if x in chromosome_encoding_dict else 0
            for x in self.adata.obs[chromosome_col]
        ]

        print(
            f"Chromosome encoding applied to {chromosome_col} in adata.obs and stored as 'perturbed_target_chromosome_encoding'."
        )

    def create_columns(self, col_dict, slot=Literal["var", "obs"], overwrite=False):
        """
        Create new columns in the specified slot of the adata object based on the provided dictionary.
        Parameters
        ----------
        col_dict : dict
            A dictionary containing the column names and their values.
        slot : str
            The slot to create columns in. Can be either "obs" or "var".
        overwrite : bool
            Whether to overwrite existing columns. If False, it raises an error if any column already exists.
            Default is False.
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')

        df = getattr(self.adata, slot).copy()

        # Check if the columns already exist
        column_names = set(col_dict.keys())
        existing_columns = set(df.columns)
        if any([col in existing_columns for col in column_names]) and not overwrite:
            existing_columns = column_names.intersection(existing_columns)
            raise ValueError(
                f"Columns {existing_columns} already exist in adata.{slot}. Review the column names or set overwrite=True to replace them."
            )

        if df.empty:
            raise ValueError(f"adata.{slot} is empty")
        # Create the new columns
        for column_name, column_value in col_dict.items():
            df[column_name] = column_value
            print(f"Column {column_name} added to adata.{slot}")

        setattr(self.adata, slot, df)

    def rename_columns(self, name_dict, slot=Literal["var", "obs"]):
        """
        Rename the columns of the specified slot of the adata object based on the provided dictionary.
        Parameters
        ----------
        name_dict : dict
            A dictionary containing the old and new column names.
        slot : str
            The slot to rename columns in. Can be either "obs" or "var".
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')

        df = getattr(self.adata, slot)

        # Check if the columns to be renamed exist in the DataFrame
        old_names = set(name_dict.keys())
        df_columns = set(df.columns)
        if not all([col in df_columns for col in old_names]):
            missing_cols = old_names - df_columns
            raise ValueError(f"Columns {missing_cols} not found in adata.{slot}")
        if df.empty:
            raise ValueError(f"adata.{slot} is empty")
        # Rename the columns
        df = df.rename(columns=name_dict)
        print(f"Renamed columns in adata.{slot}: {name_dict}")

        setattr(self.adata, slot, df)

    def replace_entries(self, slot=Literal["var", "obs"], column=None, map_dict=None):
        """
        Replace entries in a column of the named slot of the adata object. Note that values are replaced in the defined order.
        Parameters
        ----------
        slot : str
            The slot to replace entries in. Can be either "obs" or "var".
        column : str
            The name of the column to replace entries in.
        map_dict : dict
            A dictionary mapping old values to new values.
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')
        df = getattr(self.adata, slot)
        if column not in df.columns:
            raise ValueError(f"Column {column} not found in adata.{slot}")
        if df[column].empty:
            raise ValueError(f"Column {column} is empty in adata.{slot}")
        if map_dict is None:
            raise ValueError("map_dict must be provided")
        if not isinstance(map_dict, dict):
            raise ValueError("map_dict must be a dictionary")

        for old_val, new_val in map_dict.items():
            if df[column].str.upper().str.contains(old_val.upper()).any():
                df[column] = (
                    df[column]
                    .str.upper()
                    .str.replace(old_val.upper(), new_val, regex=True)
                )
                print(
                    f"Replaced '{old_val}' with '{new_val}' in column {column} of adata.{slot}"
                )
            else:
                raise ValueError(
                    f"Column {column} has no entries matching {old_val} in adata.{slot}. Check the map_dict."
                )

        setattr(self.adata, slot, df)

    def map_values_from_column(self, ref_col, target_col, map_dict):
        """
        Replace values in target_col based on corresponding values in ref_col using ref_value and target_value.
        """

        df = self.adata.obs

        if ref_col not in df.columns:
            raise ValueError(f"Column {ref_col} not found in adata.obs")
        if df[ref_col].empty:
            raise ValueError(f"Column {ref_col} is empty in adata.obs")

        if target_col not in df.columns:
            df[target_col] = np.nan  # Create target_col if it doesn't exist
            print(f"Column {target_col} created in adata.obs")

        # Ensure target_col is a string type
        df[target_col] = df[target_col].astype(str)

        for ref_value, target_value in map_dict.items():
            if ref_value not in df[ref_col].values:
                print(
                    f"Value {ref_value} not found in column {ref_col} of adata.obs. Skipping this entry."
                )

            df.loc[df[ref_col] == ref_value, target_col] = target_value
            print(
                f"Mapped value {ref_value} in column {ref_col} to {target_value} in column {target_col} of adata.obs"
            )

        # Update the adata.obs with the modified DataFrame
        setattr(self.adata, "obs", df)

    def remove_entries(self, slot=Literal["var", "obs"], column=None, to_remove=None):
        """
        Remove entries in a column of the named slot of the adata object.
        Parameters
        ----------
        slot : str
            The slot to remove entries from. Can be either "obs" or "var".
        column : str
            The name of the column to remove entries from.
        to_remove : str
            The value to remove. Must be a regex-like string.
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')
        df = getattr(self.adata, slot)
        if column not in df.columns:
            raise ValueError(f"Column {column} not found in adata.{slot}")
        if df[column].empty:
            raise ValueError(f"Column {column} is empty in adata.{slot}")

        # remove the entries from adata
        if df[column].str.contains(to_remove).any():
            entries_to_remove = df[column].str.contains(to_remove, regex=True)
            self.adata = self.adata[~entries_to_remove]
        else:
            raise ValueError(
                f"Column {column} has no entries matching {to_remove} in adata.{slot}"
            )

        print(
            f"Removed {sum(entries_to_remove)} entries {to_remove} from column {column} of adata.{slot}"
        )

    def remove_na(self, slot=Literal["var", "obs"], column=None):
        """
        Remove NA entries in a column of the named slot of the adata object.
        Parameters
        ----------
        slot : str
            The slot to remove NA entries from. Can be either "obs" or "var".
        column : str
            The name of the column to remove NA entries from.
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')
        df = getattr(self.adata, slot)
        if column not in df.columns:
            raise ValueError(f"Column {column} not found in adata.{slot}")
        if df[column].empty:
            raise ValueError(f"Column {column} is empty in adata.{slot}")

        # remove the NA entries from adata
        if df[column].isna().any():
            na_entries = df[column].isna()
            self.adata = self.adata[~na_entries]
            print(
                f"Removed {sum(na_entries)} NA entries from column {column} of adata.{slot}"
            )
        else:
            print(f"Column {column} has no NA entries in adata.{slot}")

    @staticmethod
    def remove_version_from_genes(df, column, sep="."):
        """
        Remove version numbers from gene symbols or ENSG IDs in a column of adata.var or adata.obs.

        Args:
            df: DataFrame with gene symbols.
            column: Name of the column containing gene symbols/ENSG IDs.
            sep: Separator used between the gene symbols/ENSG IDs and the version (default is ".").
        """

        if column not in df.columns:
            raise ValueError(f"Column {column} not found in the df")
        if df[column].empty:
            raise ValueError(f"Column {column} is empty in the df")

        df[column] = df[column].str.split(sep).str[0]
        
        print(f"Removed version numbers from {column}")

        return df

    def count_entries(
        self,
        slot=Literal["var", "obs"],
        input_column=None,
        count_column_name=None,
        sep="|",
    ):
        """
        Count the number of entries (e.g. number of perturbations in a cell) in a column of the named slot of the adata object.
        Parameters
        ----------
        slot : str
            The slot to count entries in. Can be either "obs" or "var".
        input_column : str
            The name of the column to count entries in.
        count_column_name : str
            The name of the column to store the count of entries.
        sep : str
            The separator used to split the entries in the column. The default is '|'.
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')

        df = getattr(self.adata, slot)
        if input_column not in df.columns:
            raise ValueError(f"Column {input_column} not found in adata.{slot}")
        if df[input_column].empty:
            raise ValueError(f"Column {input_column} is empty in adata.{slot}")
        if count_column_name is None:
            raise ValueError("count_column_name must be provided")
        if count_column_name in df.columns:
            raise ValueError(
                f"Column {count_column_name} already exists in adata.{slot}"
            )

        # Count unique entries in the column
        df[count_column_name] = [
            len(set(x.split(sep))) if x is not None else 1 for x in df[input_column]
        ]

        # if the entry contains "untreated", set the count to 0
        df.loc[
            df[input_column].str.contains("untreated", na=False), count_column_name
        ] = 0

        setattr(self.adata, slot, df)
        print(
            f"Counted entries in column {input_column} of adata.{slot} and stored in {count_column_name}"
        )

    @staticmethod
    def get_chebi_compound(compound_name):
        """
        Search for a compound in ChEBI and return its standardized name and ChEBI ID.
        Parameters:
            compound_name (str): The name of the compound to search for.
        Returns:
            dict: A dictionary containing the response from ChEBI.
        """
        
        r = requests.get(f"https://www.ebi.ac.uk/chebi/backend/api/public/es_search/?term={compound_name}&page=1&size=1")
        
        if r.ok:
            return r.json().get('results')[0].get('_source')
        else:
            print(f"Error: {r.status_code} - {r.text}")
            return None

    def standardize_compounds(self, column=None, overwrite=False):
        """
        Standardize compound names in a DataFrame column using ChEBI.

        Parameters:
            column (str): The name of the column containing compound names to be standardized.
            overwrite (bool): Whether to overwrite existing 'treatment_label' and 'treatment_id' columns. Default is False.
        """

        df = self.adata.obs

        if column is None:
            raise ValueError(
                "Column name must be provided for standardizing compounds."
            )
        if column not in df.columns:
            raise ValueError(f"Column {column} not found in adata.obs")
        if df[column].empty:
            raise ValueError(f"Column {column} is empty in adata.obs")

        # get the unique compound names from the specified column
        compound_names = df[column].dropna().unique()

        # Initialize a list to store the search results
        search_results_df = pd.DataFrame(
            columns=["original_name", "treatment_label", "treatment_id"]
        )
        mapped_compounds = []
        unmapped_compounds = []

        # Iterate over each compound name
        for compound_name in compound_names:
            # Search for the compound in ChEBI
            chebi_results = self.get_chebi_compound(compound_name)
            if chebi_results is None:
                print(f"No results found for compound '{compound_name}'")
                unmapped_compounds.append(compound_name)
                continue
            else:
                # Merge the search results with the original DataFrame
                chebi_results_df = pd.DataFrame(chebi_results)[['name', 'chebi_accession']]
                chebi_results_df['original_name'] = compound_name
                chebi_results_df = chebi_results_df.rename(
                    columns={
                        "name": "treatment_label",
                        "chebi_accession": "treatment_id",
                    }
                )
                search_results_df = pd.concat(
                    [search_results_df, chebi_results_df], ignore_index=True
                )
                mapped_compounds.append(compound_name)

        search_results_df = search_results_df.drop_duplicates()
        # Add the mapped results to the DataFrame
        if search_results_df.empty:
            print(f"None of the compounds in column {column} found in ChEBI.")
            print(f"Unmapped compounds: {unmapped_compounds}")
        else:
            df = df.merge(
                search_results_df,
                how="left",
                left_on=column,
                right_on="original_name",
            ).drop(columns=["original_name"])

            setattr(self.adata, "obs", df)
            print(f"Successfully mapped {len(mapped_compounds)}/{len(compound_names)} compounds: {mapped_compounds}")
            display(search_results_df)
            if unmapped_compounds:
                print(f"Failed to map compounds: {unmapped_compounds}")

    def standardize_genes(
        self,
        slot=Literal["var", "obs"],
        input_column=None,
        remove_version=False,
        version_sep=".",
        multiple_entries=False,
        multiple_entries_sep=None,
        keep_unmapped=False,
    ):
        """
        Standardize gene symbols or ENSG in a DataFrame column using Open Targets gene reference table combined with Ensembl outdated ID mapping.
        Args:
            slot: Which AnnData attribute to use: "var" or "obs".
            input_column: Column name containing gene symbols/ENSG IDs
            remove_version: Boolean indicating whether to remove version numbers from gene symbols/ENSG IDs (default is False)
            version_sep: Separator used between the gene symbols/ENSG IDs and the version (default is ".")
            multiple_entries: Boolean indicating whether to handle multiple entries. Default is False.
            multiple_entries_sep: Separator used between multiple entries (default is None).
            keep_unmapped: Boolean indicating whether to keep unmapped terms as-is.
                ENSG-like terms are kept in the Ensembl ID column; all other terms
                are kept in the gene symbol column.
        Returns:
            DataFrame with standardized gene symbols and ENSG IDs
        """

        df = getattr(self.adata, slot)

        # Check if the column exists in the DataFrame
        if input_column not in df.columns:
            raise ValueError(f"Column {input_column} not found in DataFrame")
        # Check if the column is empty
        if df[input_column].empty:
            raise ValueError(f"Column {input_column} is empty")

        # reset index to avoid duplicate gene symbols
        df.index.name = 'original_index'
        df = df.reset_index()

        # initialize the converted DataFrame
        conv_df = df[[input_column, 'original_index']].copy()
        conv_df['positional_index'] = range(len(conv_df))

        if multiple_entries:
            if multiple_entries_sep is None:
                raise ValueError("multiple_entries_sep must be provided if multiple_entries is True")
            conv_df[input_column] = conv_df[input_column].str.split(multiple_entries_sep)
            conv_df = conv_df.explode(input_column)

        # Remove version numbers from gene symbols/ENSG IDs
        if remove_version:
            conv_df = self.remove_version_from_genes(df=conv_df, column=input_column, sep=version_sep)

        reference = self.load_opentargets_gene_reference()
        conv_df["normalized_input_identifier"] = conv_df[input_column].map(
            self.normalize_gene_identifier
        )

        conv_df = conv_df.merge(
            reference,
            how="left",
            on="normalized_input_identifier",
        )

        control_terms = [
            "control_nontargeting",
            "control_gsh",
            "control_genedesert",
            "control_intergenic",
            "control_positive",
            "control_guideonly",
            "control_casonly",
        ]
        control_mask = conv_df[input_column].isin(control_terms)
        for column in [
            "ensembl_gene_id",
            "gene_symbol",
            "biotype",
            "gene_coord",
            "chromosome_name",
        ]:
            conv_df.loc[control_mask, column] = conv_df.loc[control_mask, input_column]

        # Handle unmapped terms: ENSG-like terms are kept in the Ensembl ID column; all other terms are kept in the gene symbol column.
        if keep_unmapped:
            missing_mask = conv_df["ensembl_gene_id"].isna()
            ensg_mask = conv_df["normalized_input_identifier"].str.match(
                r"^ENSG[0-9]+", na=False
            )
            conv_df.loc[missing_mask & ensg_mask, "ensembl_gene_id"] = conv_df.loc[
                missing_mask & ensg_mask, input_column
            ]
            conv_df.loc[missing_mask & ~ensg_mask, "gene_symbol"] = conv_df.loc[
                missing_mask & ~ensg_mask, input_column
            ]

        mapped_count = conv_df["ensembl_gene_id"].dropna().nunique()
        input_count = conv_df[input_column].dropna().nunique()
        input_genes = conv_df[input_column].dropna().unique()
        print(
            f"{'-'*50}\nSuccessfully mapped {mapped_count} out of {input_count} genes.\nInput genes: {input_genes}\n{'-'*50}"
        )

        if multiple_entries:
            # collapse the DataFrame
            conv_df = self.collapse_df(conv_df, unique_val_column='positional_index')
            conv_df = conv_df.set_index('positional_index')

        # ensure the length of the converted DataFrame is the same as the original DataFrame
        if len(conv_df) != len(df):
            raise ValueError(
                f"Length of converted DataFrame ({len(conv_df)}) does not match "
                f"length of original DataFrame ({len(df)})"
            )

        # rename the columns depending on the slot
        if slot == "obs":
            new_colnames_map = {
                "ensembl_gene_id": "perturbed_target_ensg",
                "gene_symbol": "perturbed_target_symbol",
                "biotype": "perturbed_target_biotype",
                "gene_coord": "perturbed_target_coord",
                "chromosome_name": "perturbed_target_chromosome",
            }
        elif slot == "var":
            new_colnames_map = {
                "ensembl_gene_id": "ensembl_gene_id",
                "gene_symbol": "gene_symbol",
            }
        
        # if somehow the new column names already exist in the DataFrame, rename them to avoid conflicts
        for new_col in new_colnames_map.values():
            if new_col in conv_df.columns:
                conv_df = conv_df.rename(columns={new_col: f"original_{new_col}"})
                print(f"Renamed existing column {new_col} to original_{new_col} to avoid conflicts.")

        conv_df = conv_df.rename(columns=new_colnames_map)
        conv_df = conv_df.replace("None", None)
        # keep only relevant columns
        conv_df = conv_df[list(new_colnames_map.values()) + ['original_index']]

        # drop overlapping columns in the original df to avoid conflicts when merging, but keep the "original_index" column
        out_df = df[list(set(df.columns) - set(conv_df.columns))]

        # merge the converted DataFrame to the original DataFrame
        out_df = out_df.merge(conv_df, "left", left_index=True, right_index=True)

        out_df.index.name = 'index'

        # if keep_unmapped is False, remove unmapped genes from the DataFrame
        # this subset is done on adata.obs or adata.var, so that the adata object is updated as a whole
        mapped_output_column = new_colnames_map["ensembl_gene_id"]
        if not keep_unmapped:
            mapped_mask = out_df[mapped_output_column].notna()
            removed_count = len(out_df) - mapped_mask.sum()
            if removed_count:
                print(
                    f"Removing {removed_count} unmapped genes from adata.{slot}."
                    f"Unmapped genes: {out_df.loc[~mapped_mask, input_column].unique()}"
                )

            out_df = out_df.loc[mapped_mask].copy()
            if slot == "obs":
                self.adata = self.adata[mapped_mask.to_numpy(), :].copy()
            elif slot == "var":
                self.adata = self.adata[:, mapped_mask.to_numpy()].copy()

        setattr(self.adata, slot, out_df)

    @staticmethod
    def normalize_gene_identifier(identifier):
        """Normalize gene identifiers to match the Open Targets reference table."""
        if identifier is None or pd.isna(identifier):
            return None

        normalized = re.sub(r"\s+", " ", str(identifier).strip())
        if not normalized:
            return None

        if re.match(r"^ENSG[0-9]+(?:\.[0-9]+)?$", normalized, re.IGNORECASE):
            return normalized.split(".", maxsplit=1)[0].upper()

        return normalized.upper()

    @classmethod
    def load_opentargets_gene_reference(cls, reference_path=None):
        """Load the Open Targets gene standardization reference."""
        if reference_path is None:
            reference_path = cls.opentargets_gene_reference_path

        return pd.read_parquet(
            reference_path,
            columns=[
                "normalized_input_identifier",
                "ensembl_gene_id",
                "gene_symbol",
                "biotype",
                "gene_coord",
                "chromosome_name",
            ],
        )

    def standardize_ontology(
        self,
        input_column=None,
        column_type=Literal["term_name", "term_id"],
        ontology_type=Literal["cell_type", "cell_line", "tissue", "disease"],
        overwrite=False,
    ):
        """
        Standardize ontology terms in a DataFrame column using the provided ontology type.
        Args:
            input_column: Column name containing ontology terms
            column_type: Type of the input column, either 'term_name' or 'term_id'
            ontology_type: Type of the ontology, either 'cell_type', 'cell_line', 'tissue', or 'disease'
            overwrite: Boolean indicating whether to overwrite existing columns (default is False)
        """

        df = self.adata.obs

        if input_column not in self.adata.obs.columns:
            raise ValueError(f"Column {input_column} not found in adata.obs")
        if self.adata.obs[input_column].empty:
            raise ValueError(f"Column {input_column} is empty in adata.obs")
        if column_type not in ["term_name", "term_id"]:
            raise ValueError("column_type must be either 'term_name' or 'term_id'")
        if ontology_type not in ["cell_type", "cell_line", "tissue", "disease"]:
            raise ValueError(
                "ontology_type must be one of 'cell_type', 'cell_line', 'tissue', or 'disease'"
            )

        # Select the appropriate ontology DataFrame based on the ontology type
        if ontology_type == "cell_type":
            ont_df = self.ctype_ont
            output_column_names = {"label": "cell_type_label", "id": "cell_type_id"}
        elif ontology_type == "cell_line":
            ont_df = self.cline_ont
            output_column_names = {"label": "cell_line_label", "id": "cell_line_id"}
        elif ontology_type == "tissue":
            ont_df = self.tis_ont
            output_column_names = {"label": "tissue_label", "id": "tissue_id"}
        elif ontology_type == "disease":
            ont_df = self.dis_ont
            output_column_names = {"label": "disease_label", "id": "disease_id"}

        # if inpjt column has all None values, skip the mapping
        if df[input_column].isnull().all():
            print(
                f"Column {input_column} contains only None values. Skipping ontology mapping."
            )
            return

        # get the original column values for mapping
        conv_df = df[[input_column]].drop_duplicates().reset_index(drop=True).copy()

        # rename the input column to avoid naming conflicts
        conv_df = conv_df.rename(columns={input_column: "input_column"})

        # convert the input column to lowercase for case-insensitive matching
        conv_df["input_column_lower"] = conv_df["input_column"].str.lower()

        # replace underscores with colons in the input column for matching
        if "_" in conv_df["input_column_lower"][0]:
            conv_df["input_column_lower"] = conv_df["input_column_lower"].str.replace(
                "_", ":", regex=False
            )

        if column_type == "term_id":
            # create the lower `ontology_id` column for case-insensitive matching
            ont_df["ontology_id_lower"] = ont_df["ontology_id"].str.lower()

            # map the ontology IDs to the input column
            mapping_df = conv_df.merge(
                ont_df,
                how="left",
                left_on="input_column_lower",
                right_on="ontology_id_lower",
                indicator=True,
            )

        elif column_type == "term_name":
            # create lower `name`` and `synonym`` columns for case-insensitive matching
            ont_df["name_lower"] = ont_df["name"].str.lower()
            ont_df["synonyms_lower"] = ont_df["synonyms"].str.lower()

            # create a pluralised version of the term names
            ont_df["name_lower_plural"] = ont_df["name_lower"] + "s"
            ont_df["synonyms_lower_plural"] = ont_df["synonyms_lower"] + "s"

            # concatenate the dataframes with the names, synonyms and pluralised forms into a single column for matching
            ont_df_concat = (
                pd.concat(
                    [
                        ont_df[["name_lower", "ontology_id", "name"]].assign(
                            matching_type="name"
                        ),
                        ont_df[["name_lower_plural", "ontology_id", "name"]]
                        .rename(columns={"name_lower_plural": "name_lower"})
                        .assign(matching_type="pluralised name"),
                        ont_df[["synonyms_lower", "ontology_id", "name"]]
                        .rename(columns={"synonyms_lower": "name_lower"})
                        .assign(matching_type="synonym"),
                        ont_df[["synonyms_lower_plural", "ontology_id", "name"]]
                        .rename(columns={"synonyms_lower_plural": "name_lower"})
                        .assign(matching_type="pluralised synonym"),
                    ],
                    ignore_index=True,
                )
                .dropna(subset="name_lower")
                .drop_duplicates(subset="name_lower")
            )

            # explode the name_lower column to handle synonyms in a single cell
            ont_df_concat["name_lower"] = ont_df_concat["name_lower"].str.split("|")
            ont_df_concat = ont_df_concat.explode("name_lower").drop_duplicates(
                subset="name_lower"
            )

            # map the term names to the ontology term names
            mapping_df = conv_df.merge(
                ont_df_concat,
                how="left",
                left_on="input_column_lower",
                right_on="name_lower",
                indicator=True,
            )

        else:
            raise ValueError("column_type must be either 'term_name' or 'term_id'")

        # Filter the mapping DataFrame to get the mapped and unmapped terms
        mapped_df = (
            mapping_df[mapping_df["_merge"] == "both"]
            .drop(columns="_merge")
            .drop_duplicates()
        )

        unmapped_df = (
            mapping_df[mapping_df["_merge"] == "left_only"]
            .drop(columns="_merge")
            .drop_duplicates()
        )

        if not mapped_df.empty:
            print(
                f"Mapped {len(mapped_df)} {ontology_type} ontology terms from `{input_column}` column to ontology terms"
            )
            self.print_data(mapped_df.iloc[:, :4])
        else:
            print(
                f"Warning: No {ontology_type} ontology terms could be mapped from `{input_column}` column to ontology terms. Map the terms manually or check the input column for errors."
            )
            return

        if not unmapped_df.empty:
            print(
                f"{len(unmapped_df)} {ontology_type} ontology terms from `{input_column}` column could not be mapped to ontology terms"
            )
            self.print_data(unmapped_df.iloc[:, :4])

        # Check if the output columns already exist in the DataFrame
        for col in output_column_names.values():
            if col in df.columns:
                if not overwrite:
                    raise ValueError(
                        f"Column {col} already exists in adata.obs. Set overwrite=True to replace it."
                    )
                else:
                    print(f"Overwriting column {col} in adata.obs")
                    temp_col_name = f"temp_{col}"
                    df[temp_col_name] = df[col]
                    df = df.drop(columns=col)
                    input_column = temp_col_name

        # Merge the mapped DataFrame with the original DataFrame
        df = df.merge(
            mapped_df[["input_column", "name", "ontology_id"]],
            how="left",
            left_on=input_column,
            right_on="input_column",
        )
        # Rename the columns to match the output schema
        df = df.rename(
            columns={
                "name": output_column_names["label"],
                "ontology_id": output_column_names["id"],
            }
        ).drop(columns="input_column")

        setattr(self.adata, "obs", df)

    def match_schema_columns(
        self,
        slot=Literal["var", "obs"],
    ):
        """
        Match the columns of the specified slot of the adata object to the provided schema.
        Parameters
        ----------
        slot : str
            The slot to match columns in. Can be either "obs" or "var".
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')

        df = getattr(self.adata, slot)
        if df.empty:
            raise ValueError(f"adata.{slot} is empty")

        schema = self.obs_schema if slot == "obs" else self.var_schema

        schema_columns = schema.to_schema().columns.keys()

        df = df[schema_columns]

        setattr(self.adata, slot, df)

        print(f"Matched columns of adata.{slot} to the {slot+'_schema'}.")

    def validate_data(self, slot=Literal["var", "obs"], verbose=True):
        """
        Validate the data in the specified slot of the adata object against the schema.
        Parameters
        ----------
        slot : str
            The slot to validate. Can be either "obs" or "var".
        verbose : bool
            Whether to print the validation results. Defaults to True.
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')

        df = getattr(self.adata, slot)
        if df.empty:
            raise ValueError(f"adata.{slot} is empty")

        if slot == "obs":
            schema = self.obs_schema
            dtype_map = self.get_schema_dtype_map(schema_cls=schema)
            # Only cast columns present in the dataframe
            dtype_map = {
                col: dtype for col, dtype in dtype_map.items() if col in df.columns
            }
            if dtype_map:
                logger.debug(
                    "Applying dtype casting on adata.%s for columns: %s",
                    slot,
                    list(dtype_map.keys()),
                )
                df = df.astype(dtype_map)
        else:
            schema = self.var_schema

        try:
            validated_obs = schema.validate(df, lazy=True)

            setattr(self.adata, slot, validated_obs)

            logger.info("adata.%s is valid according to the %s_schema.", slot, slot)
            if verbose:
                # Log a concise preview to avoid huge logs
                try:
                    logger.debug(
                        "Validated adata.%s preview (shape=%s):\n%s",
                        slot,
                        validated_obs.shape,
                        validated_obs.head(5).to_string(),
                    )
                except Exception:
                    logger.debug("Validated adata.%s (shape=%s)", slot, validated_obs.shape)
                # Keep notebook-friendly display for interactive use, if available
                try:
                    display(validated_obs)
                except Exception:
                    # If display is unavailable, we've already logged a preview
                    pass

        except pa.errors.SchemaErrors as e:
            try:
                msg = json.dumps(e.message, indent=2)
            except Exception:
                msg = str(e)
            logger.error("Validation errors for adata.%s: %s", slot, msg)

    def _get_vals(self, column):
        """
        Get the unique values of a column and convert them to a list.
        Args:
            column: str
                The name of the column to get the values from.
        """

        if self.adata is None:
            raise ValueError("adata is not loaded. Please load the data first.")

        df = self.adata.obs.copy()
        if column not in df.columns:
            raise ValueError(f"{column} is not a column in adata.obs")
        unique_values = df[column].dropna().unique()
        if len(unique_values) == 0:
            return None
        else:
            return [str(x) for x in unique_values.tolist()]

    def _get_dict_vals(self, term_id, term_label):
        """
        Get the values from adata obs and convert them to a list with dictionaries.
        Args:
            term_id: str
                The name of the term ID column.
            term_label: str
                The name of the term label column.
        """
        if self.adata is None:
            raise ValueError("adata is not loaded. Please load the data first.")

        df = self.adata.obs.copy()

        if term_id not in df.columns:
            raise ValueError(f"{term_id} is not a column in adata.obs")
        if term_label not in df.columns:
            raise ValueError(f"{term_label} is not a column in adata.obs")

        # get the values from the adata object
        df = df[[term_id, term_label]].drop_duplicates()
        # rename the columns
        df = df.rename(columns={term_id: "term_id", term_label: "term_label"})
        # convert the values to a list of dictionaries
        dict_vals = []
        for index, row in df.iterrows():
            dict_vals.append(
                {"term_id": row["term_id"], "term_label": row["term_label"]}
            )

        if len(dict_vals) == 0:
            return None

        return dict_vals

    @staticmethod
    def print_data(data):
        """
        Print the DataFrame in a readable format.
        Parameters
        ----------
        data : Any
            The data to print.
        """
        if isinstance(data, pd.DataFrame):
            if data.empty or data is None:
                print("DataFrame is empty.")
                return
            else:
                print(f"DataFrame shape: {data.shape}")
                print("-" * 50)
                pprint(data)
                print("-" * 50)
        else:
            if data is None:
                print("Data is None.")
            else:
                print("-" * 50)
                pprint(data)
                print("-" * 50)

    @staticmethod
    def collapse_df(df, unique_val_column=None, sep="|"):
        """
        Collapse a DataFrame by grouping on a unique value column and aggregating other columns.
        Parameters
        ----------
        df : DataFrame
            The DataFrame to collapse.
        unique_val_column : str
            The name of the column to collapse on.
        sep : str
            The separator to use for collapsing the values in collapsed columns.
        """
        if unique_val_column not in df.columns:
            if unique_val_column not in df.index.name:
                raise ValueError(
                    f"Column {unique_val_column} not found in the dataframe"
                )
            else:
                df[unique_val_column] = df.index
                df = df.reset_index(drop=True)

        if df[unique_val_column].empty:
            raise ValueError(f"Column {unique_val_column} is empty")

        pdf = pl.from_pandas(df)

        exploded_cols = [c for c in df.columns if c != unique_val_column]

        pdf_collapsed = pdf.group_by(unique_val_column).agg([
            pl.col(c).drop_nulls().cast(pl.String).str.join(sep)
            for c in exploded_cols
        ])

        pdf_collapsed = pdf_collapsed.to_pandas()

        print(f"Collapsed column {unique_val_column} using separator {sep}")

        return pdf_collapsed

    @staticmethod
    def convert_excel_date_to_gene(symbol):
        """
        Converts Excel-corrupted gene names (e.g. '03-Mar', 'Mar-03', '1-Sep', '01-Mar-23')
        into their correct gene symbols (e.g. 'MARCH3', 'SEPT1').

        If the input does not match a date-like corruption pattern, it
        returns None (so it can be easily filtered out).

        Args:
            symbol (str): The potentially corrupted gene symbol.

        Returns:
            str or None: Corrected gene symbol, or None if not date corruption.
        """

        month_to_gene = {
            "JAN": "JAN",
            "FEB": "FEB",
            "MAR": "MARCH",
            "APR": "APRIL",
            "MAY": "MAY",
            "JUN": "JUN",
            "JUL": "JULY",
            "AUG": "AUG",
            "SEP": "SEPT",
            "OCT": "OCT",
            "NOV": "NOV",
            "DEC": "DEC",
        }

        if not isinstance(symbol, str) or not symbol.strip():
            return None  # skip empty or non-string inputs

        s = symbol.strip().upper()

        # Pattern 1: day-month (e.g., '03-Mar', '3-Mar', '03_Mar', '3/Mar')
        match1 = re.match(r"^(\d{1,2})[-_/ ]([A-Z]{3})$", s, re.IGNORECASE)

        # Pattern 2: month-day (e.g., 'Mar-03', 'Sep-1')
        match2 = re.match(r"^([A-Z]{3})[-_/ ](\d{1,2})$", s, re.IGNORECASE)

        # Case 1: numeric day before month abbreviation
        if match1:
            num, month = match1.groups()
            if month in month_to_gene:
                return f"{month_to_gene[month]}{int(num)}"

        # Case 2: month before numeric day
        if match2:
            month, num = match2.groups()
            if month in month_to_gene:
                return f"{month_to_gene[month]}{int(num)}"

        # Case 3: full Excel-style date (e.g. '01-Mar-23')
        try:
            dt = datetime.strptime(s, "%d-%b-%y")
            month = dt.strftime("%b").upper()
            day = dt.day
            if month in month_to_gene:
                return f"{month_to_gene[month]}{day}"
        except Exception:
            pass

        # Default: not a corrupted date pattern
        return None

    @staticmethod
    def get_schema_dtype_map(schema_cls):
        """
        Create a mapping from pandera schema field names to pandas dtypes.
        Args:
            schema_cls: The schema class to extract annotations from.
        """
        dtype_map = {}
        for attr, annotation in schema_cls.__annotations__.items():
            # Get type info: Series[String], Series[Int64], Int64, etc.
            t = annotation
            # Handle Series[...] (pandera.typing.Series)
            if hasattr(t, "__origin__") and t.__origin__ == Series:
                dtype = t.__args__[0]
            else:
                dtype = t
            ## Map to pandas-accepted strings
            if dtype in [pa.String, str]:
                dtype_map[attr] = "string"
            elif dtype in [pa.Int64, Int64, int]:
                dtype_map[attr] = "Int64"
            elif dtype in [pa.Float64, float]:
                dtype_map[attr] = "float"
            elif dtype in [pa.Bool, bool]:
                dtype_map[attr] = "boolean"
            else:
                dtype_map[attr] = "object"
        return dtype_map


def parquet_to_bq_type(parquet_dtype):
    """Map parquet datatypes to BigQuery SQL types."""
    # This mapping can be extended based on your schema
    parquet_str = str(parquet_dtype).lower()
    if "int" in parquet_str:
        return "INT64"
    if "float" in parquet_str or "double" in parquet_str or "decimal" in parquet_str:
        return "FLOAT64"
    if "string" in parquet_str or "text" in parquet_str:
        return "STRING"
    if "boolean" in parquet_str or "bool" in parquet_str:
        return "BOOL"
    if "timestamp" in parquet_str:
        return "TIMESTAMP"
    if "date" in parquet_str:
        return "DATE"
    if "time" in parquet_str:
        return "TIME"
    # fallback
    return "STRING"


def generate_create_table_sql(
    table_name: str,
    schema: ibis.expr.schema.Schema,
    dataset_name: str = None,
    partition_column: str = None,
    partition_range_start: int = None,
    partition_range_end: int = None,
    partition_range_interval: int = None,
    cluster_columns: list = None,
):
    # Compose the full table name with dataset if provided
    full_table_name = f"{dataset_name}.{table_name}" if dataset_name else table_name
    # Generate column definitions
    columns_sql = []
    for col_name, col_type in schema.items():
        bq_type = parquet_to_bq_type(col_type)
        # BigQuery reserved keywords or spaces require backticks
        safe_col_name = f"`{col_name}`" if re.match(r"\W", col_name) else col_name
        columns_sql.append(f"{safe_col_name} {bq_type}")
    columns_def = ",\n  ".join(columns_sql)

    # Prepare partition clause
    partition_clause = ""
    if partition_column:
        # Validate partition range parameters
        if (
            partition_range_start is None
            or partition_range_end is None
            or partition_range_interval is None
        ):
            raise ValueError(
                "For integer range partitioning, you must specify start, end, and interval."
            )
        # BigQuery requires partition column identifier to be backticked if needed
        partition_col_safe = (
            f"`{partition_column}`"
            if re.match(r"\W", partition_column)
            else partition_column
        )
        partition_clause = (
            f"\nPARTITION BY RANGE_BUCKET({partition_col_safe}, GENERATE_ARRAY("
            f"{partition_range_start}, {partition_range_end}, {partition_range_interval}))"
        )

    # Prepare clustering clause
    cluster_clause = ""
    if cluster_columns:
        cluster_cols_safe = []
        for ccol in cluster_columns:
            cluster_cols_safe.append(f"`{ccol}`" if re.match(r"\W", ccol) else ccol)
        cluster_clause = f"\nCLUSTER BY {', '.join(cluster_cols_safe)}"

    create_table_sql = (
        f"CREATE TABLE IF NOT EXISTS {full_table_name} (\n"
        f"  {columns_def}\n"
        f"){partition_clause}{cluster_clause};"
    )
    return create_table_sql


def create_bq_table(
    project_id=None,
    dataset_name=None,
    table_name=None,
    schema=None,
    partition_column="perturbed_target_chromosome_encoding",
    partition_range_start=0,
    partition_range_end=25,
    partition_range_interval=1,
    cluster_columns=["dataset_id", "sample_id", "perturbed_target_symbol"],
):
    """Create a BigQuery table using the provided DDL SQL."""
    client = ibis.bigquery.connect(
        project_id=project_id, dataset_id=dataset_name, location="europe-west2"
    )
    ddl_sql = generate_create_table_sql(
        dataset_name=dataset_name,
        table_name=table_name,
        schema=schema,
        partition_column=partition_column,
        partition_range_start=partition_range_start,
        partition_range_end=partition_range_end,
        partition_range_interval=partition_range_interval,
        cluster_columns=cluster_columns,
    )
    try:
        client.raw_sql(ddl_sql)
        print(
            f"Table {dataset_name} created successfully in {project_id}.{dataset_name}."
        )
    except Exception as e:
        print(f"Error creating table: {e}")


def add_bq_upload_timestamp(bq_dest_table):
    """Add a timestamp column to the BigQuery table to track when data was ingested."""
    client = bigquery.Client()
    queries = [
        f"ALTER TABLE `{bq_dest_table}` ADD COLUMN IF NOT EXISTS ingested_at TIMESTAMP",
        f"ALTER TABLE `{bq_dest_table}` ALTER COLUMN ingested_at SET DEFAULT CURRENT_TIMESTAMP()",
        f"UPDATE `{bq_dest_table}` SET ingested_at = CURRENT_TIMESTAMP() WHERE TRUE",
    ]
    for sql in queries:
        client.query(sql).result()


def _upload_parquet_to_bq(
    parquet_path,
    project_id,
    bq_dataset_id: Literal["perturb_seq", "crispr", "mavedb"],
    bq_table_name: Literal["data", "metadata"],
    key_columns,
    verbose=True,
):
    """Upload a parquet file to BigQuery, defaulting project_id to BQ_PROJECT."""
    if project_id is None:
        if "BQ_PROJECT" not in os.environ:
            raise ValueError(
                "project_id must be provided or set in the BQ_PROJECT environment variable."
            )
        else:
            project_id = os.environ["BQ_PROJECT"]
            print(
                f"Using project_id from environment variable BQ_PROJECT: {project_id}"
            )
    else:
        print(f"Using provided project_id: {project_id}")

    if bq_dataset_id is None:
        raise ValueError("bq_dataset_id must be provided.")

    if bq_table_name is None:
        raise ValueError("bq_table_name must be provided.")
    if not key_columns:
        raise ValueError("key_columns must contain at least one column.")

    client = bigquery.Client()
    target_table_base = f"{project_id}.{bq_dataset_id}.{bq_table_name}"
    staging_table_id = f"{target_table_base}_staging"
    
    # get the target table schema
    target_table = client.get_table(target_table_base)
    # define the staging table schema (all STRING except ingested_at - it's added later)
    target_schema = [
        bigquery.SchemaField(col.name, "STRING")
        for col in target_table.schema
        if col.name != "ingested_at"
    ]


    if verbose:
        print(
            f"Staging table: loading `.parquet` file {parquet_path} to {staging_table_id}..."
        )
    
    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.PARQUET,
        write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
    )

    # create the staging table
    with open(parquet_path, "rb") as parquet_file:
        load_job = client.load_table_from_file(
            file_obj=parquet_file,
            destination=staging_table_id,
            job_config=job_config,
            rewind=True,
        )
    load_job.result()

    dest_table = client.get_table(staging_table_id)

    if verbose:
        print(f"Staging table: loaded {dest_table.num_rows} rows to {staging_table_id}")

    # add a timestamp column to the staging table
    add_bq_upload_timestamp(staging_table_id)
    if verbose:
        print(
            f"Staging table: added ingested_at timestamp column to {staging_table_id}"
        )

    # proceed to merge
    dest_table = client.get_table(target_table_base)

    # merge staging to target
    key_columns = [key.lower() for key in key_columns]
    update_columns = [col for col in dest_table.schema if col.name not in key_columns]
    update_columns = [col.name for col in update_columns if col.name != "row_id"]
    merge_staging_to_target(
        client, staging_table_id, target_table_base, key_columns, update_columns
    )

    # delete the staging table
    client.delete_table(staging_table_id, not_found_ok=True)
    if verbose:
        print(f"Staging table: deleted {staging_table_id}")


def merge_staging_to_target(
    client, staging_table_id, target_table_id, key_columns, update_columns
):
    """
    Merge staging table (all STRING columns) into target table (typed columns).

    Staging table is assumed to contain only STRING-typed columns.
    Columns are CAST into the correct types defined in the target table schema
    during INSERT and UPDATE.
    """

    # Fetch target table schema from BigQuery
    target_schema = {field.name.lower(): field for field in client.get_table(target_table_id).schema}

    def cast_expression(col):
        """Return CAST(S.col AS <typename>) based on target schema."""
        target_field = target_schema[col.lower()]
        bq_type = target_field.field_type  # e.g., STRING, INT64, FLOAT64, BOOL, DATE
        return f"CAST(S.{col} AS {bq_type})"

    # Join condition always cast key columns to their target types
    join_condition = " AND ".join(
        [f"T.{col} = {cast_expression(col)}" for col in key_columns]
    )

    # Update set — cast each column into its target type
    update_set = ", ".join(
        [f"T.{col} = {cast_expression(col)}" for col in update_columns]
    )

    # Insert columns and properly casted values
    all_columns = key_columns + update_columns
    insert_columns = ", ".join(all_columns)
    insert_values = ", ".join([cast_expression(col) for col in all_columns])

    merge_sql = f"""
    MERGE `{target_table_id}` T
    USING `{staging_table_id}` S
    ON {join_condition}
    WHEN MATCHED THEN
      UPDATE SET {update_set}
    WHEN NOT MATCHED THEN
      INSERT ({insert_columns})
      VALUES ({insert_values})
    """

    query_job = client.query(merge_sql)
    query_job.result()
    print(f"Merge completed: staging → {target_table_id} with type-safe casting.")




def download_file(
    url: str = None, dest_path: str = None, overwrite=False, unarchive: bool = False
) -> None:
    """
    Download a file from a URL to a local destination.

    Parameters:
        url: str
            URL of the file to download.
        dest_path: str
            Destination path for the downloaded file.
        overwrite: bool
            Whether to overwrite the existing file.
        unarchive: bool
            Whether to unarchive the file if it's an archive (zip/tar.gz/tgz)
    """
    # check if the file already exists
    if os.path.exists(dest_path):
        if not overwrite:
            print(f"File {dest_path} already exists. Skipping download.")
            return
        else:
            print(f"File {dest_path} already exists. Overwriting...")
    response = requests.get(url, stream=True)
    response.raise_for_status()  # Raise an error for bad responses
    # if the destination directory does not exist, create it
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    # write the content to the destination file
    with open(dest_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    if unarchive:
        if dest_path.endswith(".zip"):
            subprocess.run(["unzip", "-o", dest_path, "-d", os.path.dirname(dest_path)])
        elif dest_path.endswith((".tar.gz", ".tgz", ".tar")):
            subprocess.run(["tar", "-xzf", dest_path, "-C", os.path.dirname(dest_path)])
        else:
            print(f"Unsupported archive format for {dest_path}. Skipping unarchive.")

    print(f"Downloaded {url} to {dest_path}")

def concatenate_parquet_files(
    parquet_dir: str,
    output_path: str,
    pattern: str = "*_curated_metadata.parquet",
    verbose: bool = True
) -> None:
    """
    Stream-concatenate multiple Parquet files in `parquet_dir` matching `pattern` into a single file at `output_path`
    without loading entire tables into memory.

    Parameters
    ----------
    parquet_dir : str
        Directory containing the Parquet files to concatenate.
    output_path : str
        Path to save the concatenated Parquet file.
    pattern : str
        Glob pattern to match Parquet files (default "*_curated_metadata.parquet").
    verbose : bool
        Whether to print progress.
    """
    parquet_files = sorted(glob.glob(f"{parquet_dir}/{pattern}"))
    if not parquet_files:
        raise ValueError(f"No parquet files found with pattern {pattern} in {parquet_dir}")

    if verbose:
        print(f"Found {len(parquet_files)} files. Initializing writer...")

    # Base schema from first file
    first_pf = pq.ParquetFile(parquet_files[0])
    base_schema = first_pf.schema_arrow
    writer = pq.ParquetWriter(output_path, base_schema)

    total_rows = 0
    try:
        for idx, fpath in enumerate(parquet_files, start=1):
            pf = pq.ParquetFile(fpath)
            if pf.schema_arrow != base_schema:
                raise ValueError(f"Schema mismatch in file {fpath}. Aborting to avoid misaligned output.")
            for batch in pf.iter_batches():
                writer.write_batch(batch)
                total_rows += batch.num_rows
            if verbose:
                print(f"[{idx}/{len(parquet_files)}] Wrote {pf.metadata.num_rows} rows from {os.path.basename(fpath)} (cumulative {total_rows})")
    finally:
        writer.close()
        if verbose:
            print(f"Completed write. Total rows: {total_rows}. Output: {output_path}")

def fetch_latest_ensg_id(ensg_list: list = None):
    """
    Fetch the latest ENSG id from Ensembl REST API for a list of ENSG ids.

    Parameters
    ----------
    ensg_list : list of str
        List of ensembl IDs to query.

    Returns
    -------
    dict
        Decoded JSON response from Ensembl containing gene information.
    """
    import requests
    import sys

    if ensg_list is None or len(ensg_list) == 0:
        raise ValueError("ensg_list cannot be empty.")

    server = "https://rest.ensembl.org"
    ext = "/archive/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # if more than 500 ids, split into chunks of 500
    chunk_size = 500
    if len(ensg_list) > chunk_size:
        df_list = []
        print(f"{len(ensg_list)} unmapped ENSG IDs identified; splitting into chunks of {chunk_size}.")
        for i in range(0, len(ensg_list), chunk_size):
            chunk = ensg_list[i : i + chunk_size]
            print(f"Processing IDs {i+1} to {min(i + chunk_size, len(ensg_list))}...")
            
            data = {"id": chunk}
            
            r = requests.post(server + ext, headers=headers, json=data)
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            df_list.append(pd.DataFrame.from_dict(r.json()))
        df = pd.concat(df_list, ignore_index=True)
    else:
        data = {"id": ensg_list}
        r = requests.post(server + ext, headers=headers, json=data)

        if not r.ok:
            r.raise_for_status()
            sys.exit()
        df = pd.DataFrame.from_dict(r.json())
        
    df = df.explode("possible_replacement")
    df = pd.concat([df, df["possible_replacement"].apply(pd.Series)], axis=1).drop(
        columns=["possible_replacement"]
    )
    df = df[df["is_current"] == ""]
    df = df[["id", "stable_id"]]

    mapping_dict = {k: v for k, v in zip(df["id"], df["stable_id"])}

    return mapping_dict
