import os
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import polars as pl
import pyarrow.parquet as pq
import json
from pprint import pprint

from pydantic import ValidationError
from typing import Literal
import pandera as pa
from pandera.typing import Series, Int64, String

from libchebipy import search
import scanpy as sc
from scanpy import AnnData
import anndata as ad
from gprofiler import GProfiler

import ibis
import re
from google.cloud import bigquery

from .perturbseq_anndata_schema import ObsSchema, VarSchema
from .unified_metadata_schema.unified_metadata_schema import Experiment

# function to add a new synonym to the ontology
def add_synonym(ontology_type=Literal["genes", "cell_types", "cell_lines", "tissues", "diseases"], ref_column=str, syn_column=str, syn_map=dict, save=True):
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
    
    if ontology_type not in ["genes", "cell_types", "cell_lines", "tissues", "diseases"]:
        raise ValueError("ontology_type must be one of 'genes', 'cell_types', 'cell_lines', 'tissues', 'diseases'")
    
    # Get the path to the ontologies directory relative to this file
    ONTOLOGIES_DIR = Path(__file__).parent / "ontologies"
    
    ont = pd.read_parquet(ONTOLOGIES_DIR / f"{ontology_type}.parquet").drop_duplicates()
    
    if ref_column not in ont.columns:
        raise ValueError(f"Column `{ref_column}` not found in `{ontology_type}` ontology")
    if syn_column not in ont.columns:
        raise ValueError(f"Column `{syn_column}` not found in `{ontology_type}` ontology")
    
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
    display(ont.loc[ont[ref_column].isin(syn_map.keys()),:])
                
    # Save the updated ontology
    if save:
        ont.to_parquet(ONTOLOGIES_DIR / f"{ontology_type}.parquet", index=False)


class CuratedDataset:

    # Get the path to the ontologies directory relative to this file
    ONTOLOGIES_DIR = Path(__file__).parent / "ontologies"

    gene_ont = pd.read_parquet(ONTOLOGIES_DIR / "genes.parquet").drop_duplicates()
    ctype_ont = pd.read_parquet(ONTOLOGIES_DIR / "cell_types.parquet").drop_duplicates()
    cline_ont = pd.read_parquet(ONTOLOGIES_DIR / "cell_lines.parquet").drop_duplicates()
    tis_ont = pd.read_parquet(ONTOLOGIES_DIR / "tissues.parquet").drop_duplicates()
    dis_ont = pd.read_parquet(ONTOLOGIES_DIR / "diseases.parquet").drop_duplicates()

    def __init__(
        self,
        obs_schema,
        var_schema,
        exp_metadata_schema,
        data_source_link=None,
        noncurated_path=None,
        curated_path=None
    ):
        """Initialize the CuratedDataset class.
        Parameters
        ----------
        obs_schema : ObsSchema
            The schema for the obs data.
        var_schema : VarSchema
            The schema for the var data.
        exp_metadata_schema : Experiment
            The schema for the experimental metadata.
        data_source_link : str, optional
            The link to the data source. The default is None.
        noncurated_path : str, optional
            The path to the non-curated h5ad data. The default is None.
        curated_path : str, optional
            The path to the curated h5ad data. The default is None.
        """

        self.obs_schema = obs_schema
        self.var_schema = var_schema

        # Initialise the experiment metadata schema and its sub-schemas
        self.exp_metadata_schema = exp_metadata_schema
        self.study_schema = exp_metadata_schema.model_fields["study"].annotation
        self.experiment_schema = exp_metadata_schema.model_fields[
            "experiment"
        ].annotation
        self.perturbation_schema = exp_metadata_schema.model_fields[
            "perturbation"
        ].annotation
        self.assay_schema = exp_metadata_schema.model_fields["assay"].annotation
        self.model_system_schema = exp_metadata_schema.model_fields[
            "model_system"
        ].annotation
        self.associated_datasets_schema = exp_metadata_schema.model_fields[
            "associated_datasets"
        ].annotation
        self.associated_diseases_schema = exp_metadata_schema.model_fields[
            "associated_diseases"
        ].annotation

        self.data_source_link = data_source_link
        self.noncurated_path = noncurated_path
        self.curated_path = curated_path

        if self.noncurated_path and not self.curated_path:
            self.curated_path = noncurated_path.replace("non_curated", "curated").replace(
                ".h5ad", "_curated.h5ad"
            )
        elif self.curated_path:
            self.noncurated_path = curated_path.replace("curated", "non_curated").replace(
                "_curated.h5ad", ".h5ad"
            )
            
        self.curated_parquet_data_path = self.curated_path.replace(".h5ad", "_data.parquet").replace('h5ad', 'parquet')
        self.curated_parquet_metadata_path = self.curated_path.replace(".h5ad", "_metadata.parquet").replace('h5ad', 'parquet')

        # Initialise adata object
        self.adata = None
        # Initialise a template for experiment metadata
        self.exp_metadata = {
            k: {} for k in self.exp_metadata_schema.model_fields.keys()
        }
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
        Load the adata from curated or non-curated path.
        """
        if curated:
            if os.path.exists(self.curated_path):
                print(f"Loading data from {self.curated_path}")
                self.adata = sc.read_h5ad(self.curated_path)
            else:
                raise ValueError(f"File {self.curated_path} does not exist. Check the path.")
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

    def add_exp_metadata_as_uns(self):
        """
        Add the experiment metadata to the adata.uns dictionary.
        """
        if self.adata is None:
            raise ValueError("adata is not loaded. Please load the data first.")    

        # Add the experiment metadata to the adata.uns dictionary

        self.adata.uns['experiment_metadata'] = json.dumps(self.exp_metadata)

        print("Experiment metadata added to adata.uns")

    def save_curated_data_h5ad(self):
        """Save the curated data to a h5ad file."""

        ad.settings.allow_write_nullable_strings = True  # Allow nullable strings in adata
        
        adata = self.adata

        if adata is None:
            raise ValueError("adata is not loaded. Please load the data first.")

        # check if the base directory exists, if not create it
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

    def save_curated_data_parquet(self, split_metadata=False, save_metadata_only=False):
        """Save the curated data to a parquet file ready for BigQuery ingestion.

        Parameters
        ----------
        split_metadata : bool, optional
            Whether to split the data and metadata into two separate files (default is False).
        save_metadata_only : bool, optional
            Whether to save only the metadata and skip saving the data (default is False).
        """

        adata = self.adata

        if adata is None:
            raise ValueError("adata is not loaded. Please load the data first.")

        # check if the base directory exists, if not create it
        if not os.path.exists(os.path.dirname(self.curated_path)):
            os.makedirs(os.path.dirname(self.curated_path))

        polars_schema = self.polars_schema_from_pandera_model()

        # Concatenate the adata.obs and uns_df DataFrames
        full_metadata_df = adata.obs
        
        metadata_columns = full_metadata_df.columns.to_list()
        id_columns = metadata_columns[0:2]
        
        # Process features (e.g. genes or scores) in chunks
        feature_colnames = adata.var_names.tolist()

        if not split_metadata:
            parquet_path = self.curated_path.replace(".h5ad", "_unified.parquet").replace('h5ad', 'parquet')
            if os.path.exists(parquet_path):
                raise FileExistsError(f"File {parquet_path} already exists. Skipping write.")
            if not os.path.exists(os.path.dirname(parquet_path)):
                os.makedirs(os.path.dirname(parquet_path))
            else:
                X_df = adata.to_df()

                full_data_df = pd.concat([full_metadata_df.reset_index(drop=True), X_df.reset_index(drop=True)], axis=1, ignore_index=False)
                full_data_df = pl.from_pandas(full_data_df, schema_overrides=polars_schema)

                # convert to pyarrow for efficient streaming to parquet
                full_data_df = full_data_df.to_arrow()
                #TODO: if this works, try also using polars native parquet writer pl.write_parquet
                writer = pq.ParquetWriter(parquet_path, full_data_df.schema)
                print(f"Created ParquetWriter and started writing full data to {parquet_path}")
                writer.write_table(full_data_df)
                writer.close()
                print(f"✅ Unified data saved to {parquet_path}")

        # if split_metadata is True, save metadata and data separately
        else:
            if not os.path.exists(os.path.dirname(self.curated_parquet_data_path)):
                os.makedirs(os.path.dirname(self.curated_parquet_data_path))
            if not os.path.exists(os.path.dirname(self.curated_parquet_metadata_path)):
                os.makedirs(os.path.dirname(self.curated_parquet_metadata_path))

            if os.path.exists(self.curated_parquet_data_path) or os.path.exists(self.curated_parquet_metadata_path):
                print(f"Files {self.curated_parquet_data_path} or {self.curated_parquet_metadata_path} already exist. Skipping write.")
                return
            # if save_metadata_only is True, save only the metadata and skip saving the data
            if save_metadata_only:
                # Convert metadata_only_df to Polars DataFrame
                metadata_only_df = pl.from_pandas(full_metadata_df, schema_overrides=polars_schema)
                # Write metadata_only_df to parquet
                metadata_only_df.write_parquet(self.curated_parquet_metadata_path)
                print(f"✅ Metadata saved to {self.curated_parquet_metadata_path}")
                return
                
            else:
                # Convert metadata_only_df to Polars DataFrame
                metadata_only_df = pl.from_pandas(full_metadata_df, schema_overrides=polars_schema)
                # Write metadata_only_df to parquet
                metadata_only_df.write_parquet(self.curated_parquet_metadata_path)
                print(f"✅ Metadata saved to {self.curated_parquet_metadata_path}")
                
                print("Processing data...")
                X_df = adata.to_df()

                full_data_df = pd.concat([full_metadata_df.reset_index(drop=True), X_df.reset_index(drop=True)], axis=1, ignore_index=False)
                full_data_df = pl.from_pandas(full_data_df, schema_overrides=polars_schema)

                # Select only the ID and feature columns for data file
                data_subset_df = full_data_df.select(id_columns + feature_colnames)

                data_subset_df = data_subset_df.unpivot(
                    on=feature_colnames,
                    index=id_columns,
                    variable_name="score_name",
                    value_name="score_value"
                )

                # save data_subset_df to parquet
                print(f"Saving data to {self.curated_parquet_data_path}...")
                data_subset_df.write_parquet(self.curated_parquet_data_path)
                print(f"✅ Data saved to {self.curated_parquet_data_path}")

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
            "1": 1, "2": 2, "3": 3, "4": 4, "5": 5,
            "6": 6, "7": 7, "8": 8, "9": 9, "10": 10,
            "11": 11, "12": 12, "13": 13, "14": 14, "15": 15,
            "16": 16, "17": 17, "18": 18, "19": 19, "20": 20,
            "21": 21, "22": 22, "X": 23, "Y": 24
        }
        
        # Apply the mapping to the chromosome column
        self.adata.obs["perturbed_target_chromosome_encoding"] = [chromosome_encoding_dict[x] if x in chromosome_encoding_dict else 0 for x in self.adata.obs[chromosome_col]]
        
        print(f"Chromosome encoding applied to {chromosome_col} in adata.obs and stored as 'perturbed_target_chromosome_encoding'.")
    
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
            Whether to overwrite existing columns. If False, raises an error if any column already exists.
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

    def replace_entries(
        self,
        slot=Literal["var", "obs"],
        column=None,
        map_dict=None
    ):
        """
        Replace entries in a column of the named slot of adata object. Note that values are replaced in the defined order.
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
                df[column] = df[column].str.upper().str.replace(old_val.upper(), new_val, regex=True)
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

        # Ensure target_col is string type
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
        Remove entries in a column of the named slot of adata object.
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
        Remove NA entries in a column of the named slot of adata object.
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

    def map_symbol_to_ensg(self, symbol_column, ensg_column):
        """
        Map gene symbols to ENSEMBL IDs using the gene ontology.
        Parameters
        ----------
        symbol_column : str
            The name of the column containing gene symbols to be mapped.
        ensg_column : str
            The name of the column to store the mapped ENSEMBL IDs.
        """
        if symbol_column in self.adata.obs.columns:
            self.adata.obs[ensg_column] = self.adata.obs[symbol_column].map(
                self.gene_ont.drop_duplicates(subset=["symbol"]).set_index("symbol")[
                    "ensembl_gene_id"
                ]
            )
            print(f"Mapped {symbol_column} to ENSEMBL IDs")
            if self.adata.obs[ensg_column].isnull().any():
                print(
                    f"Warning: Some values in {symbol_column} could not be mapped to ENSEMBL IDs."
                )
                print(
                    f"Unmapped values: {self.adata.obs[symbol_column][self.adata.obs[ensg_column].isnull()]}"
                )
        else:
            print(f"Column {symbol_column} not found in adata.obs")

    def remove_version_from_genes(self, slot, column, sep="."):
        """
        Remove version numbers from gene symbols or ENSG IDs in a column of adata.var or adata.obs.

        Args:
            slot: Which AnnData attribute to use: "var" or "obs".
            column: Name of column containing gene symbols/ENSG IDs.
            sep: Separator used between the gene symbols/ENSG IDs and the version (default is ".").
        """
        if slot not in ["var", "obs"]:
            raise ValueError('slot must be either "var" or "obs"')
        df = getattr(self.adata, slot)
        if column not in df.columns:
            raise ValueError(f"Column {column} not found in adata.{slot}")
        if df[column].empty:
            raise ValueError(f"Column {column} is empty in adata.{slot}")

        df[column] = df[column].str.split(sep).str[0]

        setattr(self.adata, slot, df)

        print(f"Removed version numbers from {column} in adata.{slot}")

    def count_entries(
        self,
        slot=Literal["var", "obs"],
        input_column=None,
        count_column_name=None,
        sep="|",
    ):
        """
        Count the number of entries (e.g. number of perturbations in a cell) in a column of the named slot of adata object.
        Parameters
        ----------
        slot : str
            The slot to count entries in. Can be either "obs" or "var".
        input_column : str
            The name of the column to count entries in.
        count_column_name : str
            The name of the column to store the count of entries.
        sep : str
            The separator used to split the entries in the column. Default is '|'.
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
        df.loc[df[input_column].str.contains("untreated", na=False), count_column_name] = 0

        setattr(self.adata, slot, df)
        print(
            f"Counted entries in column {input_column} of adata.{slot} and stored in {count_column_name}"
        )

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

        # Iterate over each compound name
        for compound_name in compound_names:

            # Search for the compound in ChEBI
            chebi_results = search(compound_name)

            # Iterate over the search results
            for chebi_entity in chebi_results:
                for chebi_name_entry in chebi_entity.get_names():
                    standardized_name = chebi_name_entry.get_name()
                    # Check if the standardized name matches the input compound name (case-insensitive)
                    if standardized_name.lower() == compound_name.lower():
                        # Append the result to the list
                        compound_result_df = pd.DataFrame(
                            {
                                "original_name": compound_name,  # Original input compound name
                                "treatment_label": chebi_entity.get_name(),  # Standardized compound name
                                "treatment_id": chebi_entity.get_id(),  # ChEBI ID
                            },
                            index=[0],
                        )
                        search_results_df = pd.concat(
                            [search_results_df, compound_result_df], ignore_index=True
                        )

                        # If a match is found, print the standardized name and break the loop
                        if not compound_result_df.empty:
                            print(
                                f"Found standardized name for compound '{compound_name}': {chebi_entity.get_name()} (ChEBI ID: {chebi_entity.get_id()})"
                            )
                            break
                        else:
                            print(
                                f"No standardized name found for compound '{compound_name}'"
                            )

        # if any results were found, merge them back to the original DataFrame
        if not search_results_df.empty:
            if overwrite:
                if "treatment_label" in df.columns or "treatment_id" in df.columns:
                    print(
                        "Warning: Overwriting existing 'treatment_label' and 'treatment_id' columns."
                    )
                    temp_col = f"temp_{column}"
                    df[temp_col] = df[column]
                    for col in ["treatment_label", "treatment_id"]:
                        if col in df.columns:
                            df = df.drop(columns=[col])
                    column = temp_col
            else:
                if "treatment_label" in df.columns or "treatment_id" in df.columns:
                    raise ValueError(
                        "'treatment_label' and/or 'treatment_id' columns already exist. Set overwrite=True to replace them."
                    )
            # Merge the search results with the original DataFrame
            df = df.merge(
                search_results_df,
                how="left",
                left_on=column,
                right_on="original_name",
            ).drop(columns=["original_name"])

            # Update the adata.obs with the standardized DataFrame
            setattr(self.adata, "obs", df)

            print(
                f"Standardized compound names in column '{column}' and added 'treatment_label' and 'treatment_id' columns."
            )
        else:
            raise ValueError(
                f"No standardized names found for compounds in column '{column}'."
            )

    def standardize_genes(
        self,
        slot=Literal["var", "obs"],
        input_column=None,
        input_column_type=Literal["gene_symbol", "ensembl_gene_id"],
        remove_version=False,
        version_sep=".",
        multiple_entries=False,
        multiple_entries_sep=None,
    ):
        """
        Standardize gene symbols or ENSG in a DataFrame column using gprofiler.
        Args:
            slot: Which AnnData attribute to use: "var" or "obs".
            input_column: Column name containing gene symbols/ENSG IDs
            input_column_type: Type of the input column, either 'gene_symbol' or 'ensembl_gene_id'
            remove_version: Boolean indicating whether to remove version numbers from gene symbols/ENSG IDs (default is False)
            version_sep: Separator used between the gene symbols/ENSG IDs and the version (default is ".")
        Returns:
            DataFrame with standardized gene symbols and ENSG IDs
        """

        df = getattr(self.adata, slot)

        # Check if the column is gene symbol or ENSG
        if input_column_type not in ["gene_symbol", "ensembl_gene_id"]:
            raise ValueError(
                "Input column type must be either 'gene_symbol' or 'ensembl_gene_id'"
            )
        # Check if the column exists in the DataFrame
        if input_column not in df.columns:
            raise ValueError(f"Column {input_column} not found in DataFrame")
        # Check if the column is empty
        if df[input_column].empty:
            raise ValueError(f"Column {input_column} is empty")

        # inititalise the converted DataFrame
        conv_df = df[[input_column]].copy()
        conv_df.index = conv_df.index.rename("index")

        # Explode the column if it contains multiple entries
        if multiple_entries:
            if multiple_entries_sep is None:
                raise ValueError(
                    "multiple_entries_sep must be provided if multiple_entries is True"
                )
            conv_df[input_column] = conv_df[input_column].str.split(
                multiple_entries_sep
            )
            conv_df = conv_df.explode(input_column)

        # Remove version numbers from gene symbols/ENSG IDs
        if remove_version:
            self.remove_version_from_genes(
                slot=slot, column=input_column, sep=version_sep
            )

        # map the synonyms to the gene symbols
        if input_column_type == "gene_symbol":
            conv_df = self.map_synonyms_to_symbols(
                df=conv_df, symbol_column=input_column, gene_ont=self.gene_ont
            )

        # Convert the gene symbols/ENSG IDs to the target namespace (ENSG)
        gp = GProfiler(return_dataframe=True)
        gp_result = gp.convert(
            organism="hsapiens",
            query=list(set(conv_df[input_column])),
            target_namespace="ENSG",
        )

        # Filter the converted DataFrame to keep only the rows with a single conversion
        gp_result = gp_result[gp_result["n_converted"] == 1].drop(
            columns=["n_incoming", "n_converted"]
        )

        # Replace "None" with control entries
        gp_result.loc[gp_result["incoming"].str.contains("control"), "converted"] = gp_result.loc[gp_result["incoming"].str.contains("control"), "incoming"]
        gp_result.loc[gp_result["incoming"].str.contains("control"), "name"] = gp_result.loc[gp_result["incoming"].str.contains("control"), "incoming"]

        print(
            f"Converted {len(gp_result[gp_result['converted'] != 'None'])}/{len(gp_result)} gene symbols/ENSG IDs to standardized gene symbols/ENSG IDs\n{'-' * 50}"
        )

        # merge the converted DataFrame with the original DataFrame
        conv_df_idx = conv_df.index
        conv_df = (
            conv_df.merge(
                gp_result[["incoming", "converted", "name"]].drop_duplicates(),
                how="left",
                left_on=input_column,
                right_on="incoming",
                indicator=True,
            )
            .drop(columns=[input_column])
            .set_index(conv_df_idx)
        )

        # add the biotype and gene_coord columns
        conv_df = (
            conv_df.merge(
                self.gene_ont[["ensembl_gene_id", "biotype", "gene_coord", "chromosome_name"]]
                .drop_duplicates()
                .dropna(subset=["ensembl_gene_id"]),
                how="left",
                left_on="converted",
                right_on="ensembl_gene_id",
            )
            .drop(columns=["ensembl_gene_id"])
            .set_index(conv_df_idx)
        )

        if multiple_entries:
            # collapse the DataFrame to get the original column back
            conv_df = self.collapse_df(conv_df, unique_val_column="index")

        # ensure the length of the converted DataFrame is the same as the original DataFrame
        if len(conv_df) != len(df):
            raise ValueError(
                f"Length of converted DataFrame ({len(conv_df)}) does not match length of original DataFrame ({len(df)})"
            )

        # if some incoming ENSGs were not converted, keep them as is
        if (
            input_column_type == "ensembl_gene_id"
            and (conv_df["converted"] == "None").any()
        ):
            conv_df.loc[conv_df["converted"] == "None", "converted"] = conv_df.loc[
                conv_df["converted"] == "None", "incoming"
            ]

        # rename the columns depending on the slot
        if slot == "obs":
            new_colnames_map = {
                "converted": "perturbed_target_ensg",
                "name": "perturbed_target_symbol",
                "biotype": "perturbed_target_biotype",
                "gene_coord": "perturbed_target_coord",
                "chromosome_name": "perturbed_target_chromosome",
            }
        elif slot == "var":
            new_colnames_map = {"converted": "ensembl_gene_id", "name": "gene_symbol"}

        conv_df = conv_df.rename(columns=new_colnames_map)

        # drop irrelevant columns
        conv_df = conv_df[new_colnames_map.values()]

        # drop overlapping columns in original df to avoid conflicts when merging
        df = df[list(set(df.columns) - set(conv_df.columns))]
        df = df.merge(conv_df, "left", left_index=True, right_index=True)

        # replace "None" strings returned by gprofiler with None
        df = df.replace("None", None)

        # rename index to index
        df.index = df.index.rename("index")

        # replace slot with the converted DataFrame
        setattr(self.adata, slot, df)

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

        # get the original column values for mapping
        conv_df = df[[input_column]].drop_duplicates().reset_index(drop=True).copy()

        # rename the input column to avoid naming conflicts
        conv_df = conv_df.rename(columns={input_column: "input_column"})

        # convert the input column to lowercase for case-insensitive matching
        conv_df["input_column_lower"] = conv_df["input_column"].str.lower()

        if column_type == "term_id":
            # create lower `ontology_id` column for case-insensitive matching
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

            # create pluralised version of the term names
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

    def populate_exp_metadata(self):
        """
        Populate the available experiment metadata fields with values from the adata obs DataFrame.
        """
        if self.adata is None:
            raise ValueError("adata is not loaded. Please load the data first.")

        exp_metadata = self.exp_metadata

        # these can be automatically populated from adata.obs
        auto_populated_fields = {
            "experiment": {
                "treatments": self._get_dict_vals("treatment_id", "treatment_label"),
                "timepoints": self._get_vals("timepoint"),
                "perturbation_type": self._get_dict_vals(
                    "perturbation_type_id", "perturbation_type_label"
                ),
                "perturbed_target_biotype": self._get_vals("perturbed_target_biotype"),
                "number_of_perturbed_targets": len(
                    self._get_vals("perturbed_target_coord")
                ),
                "perturbed_targets": self._get_vals("perturbed_target_ensg"),
                "number_of_perturbed_entities": self.adata.obs.shape[0],
            },
            "model_system": {
                "model_system": self._get_dict_vals(
                    "model_system_id", "model_system_label"
                ),
                "tissue": self._get_dict_vals("tissue_id", "tissue_label"),
                "cell_type": self._get_dict_vals("cell_type_id", "cell_type_label"),
                "cell_line": self._get_dict_vals("cell_line_id", "cell_line_label"),
                "sex": self._get_dict_vals("sex_id", "sex_label"),
                "developmental_stage": self._get_dict_vals(
                    "developmental_stage_id", "developmental_stage_label"
                ),
            },
            "associated_diseases": self._get_dict_vals("disease_id", "disease_label"),
        }

        exp_metadata.update(auto_populated_fields)

        print("Experiment metadata populated with available fields from adata.obs:")
        self.print_data(auto_populated_fields)

    def add_exp_metadata(
        self,
        metadata_slot=Literal[
            "study",
            "experiment",
            "perturbation",
            "assay",
            "model_system",
            "associated_datasets",
        ],
        metadata: dict | list = None,
    ):
        """
        Add metadata to the experiment metadata schema.
        Parameters
        ----------
        metadata_slot : str
            The slot to add the metadata to. Can be either "study", "experiment", "perturbation", "assay", "model_system", or "associated_datasets".
        metadata : dict | list
            The metadata to add. It can be a dictionary or a list of dictionaries corresponding to the fields in the metadata schema.
        """

        if metadata_slot not in self.exp_metadata_schema.model_fields.keys():
            raise ValueError(
                f"metadata_slot must be one of {list(self.exp_metadata.keys())}"
            )
        if metadata is None:
            raise ValueError("Please provide metadata to add.")

        if not isinstance(metadata, dict) and not isinstance(metadata, list):
            raise ValueError("metadata must be a dictionary or a list of dictionaries.")

        # retrieve the relevant schema
        schema = getattr(self, f"{metadata_slot}_schema", None)

        # treat 'associated_datasets' and 'associated_diseases' differently, because they are lists of dictionaries
        if metadata_slot not in ["associated_datasets", "associated_diseases"]:

            # if metadata is partially automatically populated, combine it with the submitted metadata
            if metadata_slot in ["experiment", "model_system"]:
                metadata.update(self.exp_metadata[metadata_slot])

            # validate the metadata chunk against the sub-schema
            try:
                validated = schema.model_validate(metadata)

                # update the exp_metadata with the chunk
                self.exp_metadata[metadata_slot].update(validated.model_dump())

                print(f"Metadata for '{metadata_slot}' successfully validated:")
                self.print_data(validated.model_dump())

            except ValidationError as e:
                print(e)

        else:

            if metadata_slot == "associated_datasets":

                subschema = schema.__args__[0]

            elif metadata_slot == "associated_diseases":
                subschema = schema.__args__[0].__args__[0]
                metadata = self.exp_metadata[metadata_slot]

            if not isinstance(subschema, type):
                raise ValueError(
                    f"There is something wrong with the schema for {metadata_slot}."
                )

            # since 'associated_datasets' and 'associated_diseases' are lists of dictionaries, we need to validate each dictionary in the list
            validated_list = []
            for entry in metadata:
                if not isinstance(entry, dict):
                    raise ValueError(
                        f"Each element in {metadata_slot} must be a dictionary."
                    )

                try:
                    # validate the entry against the sub-schema
                    validated = subschema.model_validate(entry)
                    validated_list.append(validated.model_dump())

                except ValidationError as e:
                    print(f"Validation error for entry {entry}: {e}")

            # update the exp_metadata with the validated list
            self.exp_metadata[metadata_slot] = validated_list
            print(f"Metadata for '{metadata_slot}' successfully validated:")
            self.print_data(validated_list)

    def validate_exp_metadata(self):
        """
        Validate the experiment metadata against the schema.
        """
        if self.exp_metadata is None:
            raise ValueError(
                "Experiment metadata is not populated. Please populate it first by calling populate_exp_metadata() and add_exp_metadata() methods."
            )

        try:
            validated = self.exp_metadata_schema.model_validate(self.exp_metadata)
            print("Experiment metadata successfully validated:")
            self.print_data(validated.model_dump())
        except ValidationError as e:
            print(f"Validation error: {e}")
        
    def validate_data(
        self,
        slot=Literal["var", "obs"],
    ):
        """
        Validate the data in the specified slot of the adata object against the schema.
        Parameters
        ----------
        slot : str
            The slot to validate. Can be either "obs" or "var".
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')

        df = getattr(self.adata, slot)
        if df.empty:
            raise ValueError(f"adata.{slot} is empty")        
        
        if slot == "obs":
            schema = self.obs_schema
            dtype_map = self.get_schema_dtype_map(schema_cls = schema)
            # Only cast columns present in the dataframe
            dtype_map = {col: dtype for col, dtype in dtype_map.items() if col in df.columns}
            df = df.astype(dtype_map)
        else:
            schema = self.var_schema            

        try:
            validated_obs = schema.validate(df, lazy=True)

            setattr(self.adata, slot, validated_obs)

            print(f"adata.{slot} is valid according to the {slot}_schema.")
            print("Validated data:")
            display(validated_obs)

        except pa.errors.SchemaErrors as e:
            print(json.dumps(e.message, indent=2))
        
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

        # get the values from adata object
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
    def map_synonyms_to_symbols(df, symbol_column=None, gene_ont=None):
        """
        Map synonyms to gene symbols using the gene ontology.
        Parameters
        ----------
        symbol_column : str
            The name of the column containing gene symbols to be mapped.
        """
        if symbol_column not in df.columns:
            raise ValueError(f"Column {symbol_column} not found in the dataframe")
        if df[symbol_column].empty:
            raise ValueError(f"Column {symbol_column} is empty in the dataframe")
        if gene_ont is None:
            raise ValueError("gene_ont must be provided")

        # remove entries with no synonyms
        map_df = gene_ont[~gene_ont["synonyms"].isna()]
        # explode the synonyms column
        map_df["synonyms"] = map_df["synonyms"].str.split("|")
        map_df = map_df.explode("synonyms")
        # remove unnecessary columns
        map_df = (
            map_df[["synonyms", "symbol", "ensembl_gene_id"]]
            .drop_duplicates(subset=["synonyms"])
            .dropna(subset=["synonyms"])
        )
        # map the synonyms to the gene symbols
        df[symbol_column] = (
            df[symbol_column]
            .map(map_df.set_index("synonyms")["symbol"])
            .fillna(df[symbol_column])
        )

        print(
            f"Mapped potential synonyms in {symbol_column} of the provided dataframe to gene symbols"
        )

        return df

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

        exploded_cols = ["incoming", "converted", "name", "biotype", "gene_coord", "chromosome_name"]

        df = df.groupby([unique_val_column]).agg(
            {
                col: lambda x: sep.join(x) if x.notna().all() else None
                for col in exploded_cols
            }
        )

        print(f"Collapsed column {unique_val_column} using separator {sep}")

        return df

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

def map_metadata_score_columns_biogrid(screen_df, metadata_df, screen_id):

    # score columns in metadata are named SCORE.1_TYPE, SCORE.2_TYPE, etc.
    score_cols_meta = [f"SCORE.{e}_TYPE" for e in range(1, 6)]

    # create a DataFrame with screen id and score columns
    meta_subset_df = metadata_df[["#SCREEN_ID"] + score_cols_meta]

    # melt the DataFrame to have a long format with screen id, score type, and score value
    melt_meta_subset_df = meta_subset_df.melt(
        id_vars="#SCREEN_ID",
        value_vars=score_cols_meta,
        var_name="source_col",
        value_name="score_type",
    )

    # remove the '_TYPE' suffix from the source_col because SCORE columns in screen data do not have this suffix
    melt_meta_subset_df["source_col"] = melt_meta_subset_df["source_col"].str.rstrip(
        "_TYPE"
    )

    # filter the melted DataFrame for the given screen id and create a mapping dictionary
    # where the keys are the source_col (SCORE.1, SCORE.2, etc.) and the values are the score_type (e.g. 'MaGeCK Score')
    map_dict = (
        melt_meta_subset_df.loc[melt_meta_subset_df["#SCREEN_ID"] == screen_id]
        .set_index("source_col")["score_type"]
        .to_dict()
    )
    # new columns to be returned by the function
    new_columns = map_dict.values()

    # rename the columns in the screen DataFrame using the mapping dictionary
    screen_df = screen_df.rename(columns=map_dict)

    # drop the '-' empty column
    if "-" in screen_df.columns:
        screen_df = screen_df.drop(columns=["-"])
        new_columns = [e for e in new_columns if e != "-"]

    return [screen_df, new_columns]


def compare_metadata_biogrid(bg_metadata_df=None, gc_metadata_df=None, biogrid_screen_id=None, gemini_id=None):
    """
    Compare the metadata from biogrid and gemini-curated dataframes
    to identify any discrepancies or missing information.

    Parameters
    ----------
    bg_metadata_df : pd.DataFrame
        Biogrid metadata dataframe
    gc_metadata_df : pd.DataFrame
        Gemini-curated metadata dataframe

    Returns
    -------
    comparison_df : pd.DataFrame
        DataFrame showing the comparison between biogrid and gemini metadata
    """
    if bg_metadata_df is None or gc_metadata_df is None:
        print("Both bg_metadata_df and gc_metadata_df must be provided")
        return None
    if biogrid_screen_id is None or gemini_id is None:
        print("⚠️Provide biogrid_screen_id and gemini_id to show a subset of metadata comparison")

    gc_compar = gc_metadata_df.T
    gc_compar['col'] = gc_compar.index

    bg_compar = bg_metadata_df.T
    bg_compar.columns = ['biogrid_' + str(e) for e in bg_compar.loc['#SCREEN_ID', :]]
    bg_compar['col'] = bg_compar.index

    comparison_df = pd.merge(
        bg_compar,
        gc_compar,
        left_on='col', right_on='col', how='outer', indicator=True
    )
    comparison_df.index = comparison_df['col']
    comparison_df = comparison_df.drop(columns=['col']).sort_values(by='_merge', ascending=False)

    if biogrid_screen_id and gemini_id:
        comparison_df = comparison_df[[biogrid_screen_id, gemini_id, '_merge']]
    elif biogrid_screen_id and not gemini_id:
        comparison_df = comparison_df[
            [biogrid_screen_id] + [e for e in comparison_df.columns if e.startswith('gemini')] + ['_merge']]
    elif gemini_id and not biogrid_screen_id:
        comparison_df = comparison_df[[e for e in comparison_df.columns if e.startswith('biogrid')] + ['_merge']]

    return comparison_df

def make_adata_biogrid(biogrid_screen_path=None, bg_metadata_df=None, gc_metadata_df=None, biogrid_screen_id=None, data_modality=None, gemini_id=None, save_h5ad=True, save_dir=None):
    """
    Process the original biogrid CRISPR screen data to construct the adata object ready for downstream curation

    Parameters
    ----------
    biogrid_screen_path : Path
        Path to the biogrid screen data file
    bg_metadata_df : pd.DataFrame
        DataFrame containing the biogrid metadata for all screens from the same publication
    gc_metadata_df : pd.DataFrame
        DataFrame containing the gemini-curated metadata for the publication
    biogrid_screen_id : str
        Screen ID to be processed
    data_modality : str
        Data modality to be assigned in the adata.obs['data_modality']
    gemini_id : str
        Gemini ID of the publication for cross-referencing metadata
    save_h5ad : bool
        Whether to save the resulting AnnData object as an h5ad file

    Returns
    -------

    """

    if not biogrid_screen_path.exists():
        print(f"⚠️File {biogrid_screen_path} does not exist. Skipping...")
        return
    if bg_metadata_df.empty:
        print(f"⚠️File {bg_metadata_df} does not exist. Skipping...")
        return
    if gc_metadata_df.empty:
        print(f"⚠️File {gc_metadata_df} does not exist. Skipping...")
        return
    if biogrid_screen_id is None:
        print("⚠️Provide biogrid_screen_id to process the screen")
        return
    if data_modality is None:
        print("⚠️Provide data_modality to be assigned in the adata.obs['data_modality']")
        return
    if gemini_id is None:
        print("⚠️Provide gemini_id to cross-reference metadata")
        return

    # read the screen data
    screen_df = pd.read_csv(biogrid_screen_path, sep="\t")
    screen_df['significant'] = screen_df['HIT'].apply(lambda x: True if x == 'YES' else False)
    screen_df = screen_df[["#SCREEN_ID", "IDENTIFIER_ID", "OFFICIAL_SYMBOL", "significant"] + [e for e in screen_df.columns if e.startswith("SCORE.")]]

    screen_id = int(biogrid_screen_id.lstrip('biogrid_'))
    # get the column mappings for the score of interest
    screen_df, score_columns = map_metadata_score_columns_biogrid(
        screen_df, bg_metadata_df, screen_id
    )

    screen_significance_criteria = bg_metadata_df.loc[bg_metadata_df["#SCREEN_ID"] == screen_id, 'significance_criteria'].values[0]

    # convert the new columns to numeric
    screen_df[score_columns] = screen_df[score_columns].apply(pd.to_numeric, errors='coerce')

    # group_by 'IDENTIFIER_ID' and average the counts
    screen_df = (
        screen_df.groupby([e for e in screen_df.columns if e not in score_columns])
        .agg({e: "mean" for e in score_columns})
        .reset_index()
    )
    # create unique index for the adata.obs
    screen_df['index'] = "biogrid_" + screen_df['#SCREEN_ID'].astype(str) + '_' + screen_df['OFFICIAL_SYMBOL']
    # remove duplicates if any
    screen_df = screen_df[~screen_df['index'].duplicated()]
    # get X data for adata object based on the score columns
    X_df = screen_df[score_columns]
    # set the index to the unique index created above
    X_df.index = screen_df['index']
    # create adata.obs DataFrame
    OBS_df = pd.DataFrame(index=X_df.index)
    OBS_df['data_modality'] = data_modality
    OBS_df["perturbation_name"] = OBS_df.index
    OBS_df["perturbed_target_symbol"] = screen_df["OFFICIAL_SYMBOL"].to_list()
    OBS_df["significant"] = screen_df["significant"].to_list()
    OBS_df['significance_criteria'] = screen_significance_criteria
    OBS_df = OBS_df.drop_duplicates()

    # add all columns from gemini-curated gc_metadata_df as default values for the obs columns
    for col in gc_metadata_df.columns:
        if col not in OBS_df.columns:
            # OBS_df[col] = gc_metadata_df[col].values[0]
            # OBS_df[col] = gc_compar.loc[col, gemini_id]
            OBS_df[col] = gc_metadata_df.loc[gemini_id, col]

    # create adata.var DataFrame
    VAR_df = pd.DataFrame(index=score_columns)
    VAR_df["score_name"] = score_columns

    # replace None with np.nan to avoid issues with AnnData writing
    OBS_df = OBS_df.replace({None: np.nan})
    VAR_df = VAR_df.replace({None: np.nan})

    adata = sc.AnnData(X=X_df, obs=OBS_df, var=VAR_df)

    if save_h5ad:
        adata_h5ad_save_path = save_dir / f"{biogrid_screen_id}.h5ad"
        adata.write(adata_h5ad_save_path)
        print(f"✅Saved adata to {adata_h5ad_save_path}")
        return adata, adata_h5ad_save_path
    else:
        return adata


def parquet_to_bq_type(parquet_dtype):
    """Map parquet datatypes to BigQuery SQL types."""
    # This mapping can be extended based on your schema
    parquet_str = str(parquet_dtype).lower()
    if 'int' in parquet_str:
        return 'INT64'
    if 'float' in parquet_str or 'double' in parquet_str or 'decimal' in parquet_str:
        return 'FLOAT64'
    if 'string' in parquet_str or 'text' in parquet_str:
        return 'STRING'
    if 'boolean' in parquet_str or 'bool' in parquet_str:
        return 'BOOL'
    if 'timestamp' in parquet_str:
        return 'TIMESTAMP'
    if 'date' in parquet_str:
        return 'DATE'
    if 'time' in parquet_str:
        return 'TIME'
    # fallback
    return 'STRING'

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
    # Compose full table name with dataset if provided
    full_table_name = f"{dataset_name}.{table_name}" if dataset_name else table_name
    # Generate column definitions
    columns_sql = []
    for col_name, col_type in schema.items():
        bq_type = parquet_to_bq_type(col_type)
        # BigQuery reserved keywords or spaces require backticks
        safe_col_name = f"`{col_name}`" if re.match(r'\W', col_name) else col_name
        columns_sql.append(f"{safe_col_name} {bq_type}")
    columns_def = ",\n  ".join(columns_sql)

    # Prepare partition clause
    partition_clause = ""
    if partition_column:
        # Validate partition range parameters
        if (partition_range_start is None or partition_range_end is None or partition_range_interval is None):
            raise ValueError("For integer range partitioning, you must specify start, end, and interval.")
        # BigQuery requires partition column identifier to be backticked if needed
        partition_col_safe = f"`{partition_column}`" if re.match(r'\W', partition_column) else partition_column
        partition_clause = (
            f"\nPARTITION BY RANGE_BUCKET({partition_col_safe}, GENERATE_ARRAY("
            f"{partition_range_start}, {partition_range_end}, {partition_range_interval}))"
        )

    # Prepare clustering clause
    cluster_clause = ""
    if cluster_columns:
        cluster_cols_safe = []
        for ccol in cluster_columns:
            cluster_cols_safe.append(f"`{ccol}`" if re.match(r'\W', ccol) else ccol)
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
    cluster_columns=['dataset_id', 'sample_id', 'perturbed_target_symbol']
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
    queries = [f"ALTER TABLE `{bq_dest_table}` ADD COLUMN IF NOT EXISTS ingested_at TIMESTAMP",
                f"ALTER TABLE `{bq_dest_table}` ALTER COLUMN ingested_at SET DEFAULT CURRENT_TIMESTAMP()",
                f"UPDATE `{bq_dest_table}` SET ingested_at = CURRENT_TIMESTAMP() WHERE TRUE"]
    for sql in queries:
        client.query(sql).result()

def upload_parquet_to_bq(parquet_path, bq_dataset_id, bq_table_name, key_columns, verbose=True):
    client = bigquery.Client()
    target_table_base = f"{bq_dataset_id}.{bq_table_name}"
    staging_table_id = f"{target_table_base}_staging"

    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.PARQUET,
        write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
        )
    
    if verbose:
        print(f"Staging table: loading `.parquet` file {parquet_path} to {staging_table_id}...")

    # create the staging table
    with open(parquet_path, "rb") as parquet_file:
        load_job = client.load_table_from_file(
            file_obj=parquet_file,
            destination=staging_table_id,
            job_config=job_config,
            rewind=True
        )
    load_job.result()
    
    dest_table = client.get_table(staging_table_id)
    
    if verbose:
        print(f"Staging table: loaded {dest_table.num_rows} rows to {staging_table_id}")
        
    # add a timestamp column to the staging table
    add_bq_upload_timestamp(staging_table_id)
    if verbose:
        print(f"Staging table: added ingested_at timestamp column to {staging_table_id}")

    # proceed to merge
    dest_table = client.get_table(target_table_base)

    # merge staging to target
    key_columns = [key.lower() for key in key_columns]
    update_columns = [col for col in dest_table.schema if col.name not in key_columns]
    update_columns = [col.name for col in update_columns if col.name != "row_id"]
    merge_staging_to_target(client, staging_table_id, target_table_base, key_columns, update_columns)

    # delete the staging table
    client.delete_table(staging_table_id, not_found_ok=True)
    if verbose:
        print(f"Staging table: deleted {staging_table_id}")
    


def merge_staging_to_target(client, staging_table_id, target_table_id, key_columns, update_columns):
    """
    Merge staging table into target table using BigQuery MERGE statement.

    Args:
        client: BigQuery client
        staging_table_id (str): fully qualified staging table id (project.dataset.table)
        target_table_id (str): fully qualified target table id
        key_columns (list[str]): list of column names used as unique keys
        update_columns (list[str]): list of columns to update on match
    """
    join_condition = " AND ".join([f"T.{col} = S.{col}" for col in key_columns])
    update_set = ", ".join([f"T.{col} = S.{col}" for col in update_columns])
    insert_columns = ", ".join(key_columns + update_columns)
    insert_values = ", ".join([f"S.{col}" for col in key_columns + update_columns])

    merge_sql = f"""
    MERGE `{target_table_id}` T
    USING `{staging_table_id}` S
    ON {join_condition}
    WHEN MATCHED THEN
      UPDATE SET {update_set}
    WHEN NOT MATCHED THEN
      INSERT ({insert_columns}) VALUES ({insert_values})
    """

    query_job = client.query(merge_sql)
    query_job.result()  # Wait for completion
    print(f"Merge completed: staging data merged into {target_table_id}")


def process_biogrid_screen(
    biogrid_dataset_id: str = None,
    biogrid_screen_path: Path = None,
    biogrid_metadata_df: pd.DataFrame = None,
    curated_metadata_dict: dict = None,
    non_curated_h5ad_dir: str = None,
    upload_to_bq: bool = False,
    bq_dataset_id: str = None,
    bq_metadata_table_name: str = "metadata",
    bq_data_table_name: str = "data"
):
    """
    Process a BioGRID screen dataset: create AnnData, curate it, save as h5ad and Parquet, and optionally upload to BigQuery.

    Parameters:
    - biogrid_dataset_id: Identifier for the BioGRID dataset.
    - biogrid_screen_path: Path to the BioGRID screen data file.
    - biogrid_metadata_df: DataFrame containing metadata for the BioGRID dataset.
    - curated_metadata_dict: Dictionary containing curated metadata for the dataset.
    - non_curated_h5ad_dir: Directory to save non-curated h5ad files.
    - upload_to_bq: Boolean flag indicating whether to upload to BigQuery.
    - bq_dataset_id: BigQuery dataset ID for uploading data.
    - bq_metadata_table_name: Table name for metadata in BigQuery.
    - bq_data_table_name: Table name for data in BigQuery.
    """

    if upload_to_bq and bq_dataset_id is None:
        raise ValueError("bq_dataset_id must be provided if upload_to_bq is True.")

    # Step 1: Create AnnData from BioGRID screen data
    _adata, h5ad_path = make_adata_biogrid(
        biogrid_dataset_id=biogrid_dataset_id,
        biogrid_screen_path=biogrid_screen_path,
        biogrid_metadata_df=biogrid_metadata_df,
        curated_metadata_dict=curated_metadata_dict,
        save_h5ad_dir=non_curated_h5ad_dir,
    )

    # Step 2: Curate the AnnData object
    cur_data = curate_biogrid_screen(
        adata_h5ad_path=h5ad_path,
        save_curated_h5ad=True,
        save_curated_parquet=True,
        split_parquet=True,
    )

    # Step 3: Upload to BigQuery if specified
    if upload_to_bq:
        print(f"Uploading metadata Parquet to BigQuery: {cur_data.curated_parquet_metadata_path}")
        upload_parquet_to_bq(
            parquet_path=cur_data.curated_parquet_metadata_path,
            bq_dataset_id=bq_dataset_id,
            bq_table_name=bq_metadata_table_name,
            key_columns=["dataset_id", "sample_id"]
        )

        print(f"Uploading data Parquet to BigQuery: {cur_data.curated_parquet_data_path}")
        upload_parquet_to_bq(
            parquet_path=cur_data.curated_parquet_data_path,
            bq_dataset_id=bq_dataset_id,
            bq_table_name=bq_data_table_name,
            key_columns=["dataset_id", "sample_id"]
        )


    return cur_data


def curate_biogrid_screen(
        adata_h5ad_path: Path = None,
        save_curated_h5ad: bool = True,
        save_curated_parquet: bool = True,
        split_parquet=True,
):
    """
    Curate a BioGRID screen AnnData object.

    Parameters:
    - biogrid_dataset_id: Identifier of the BioGRID dataset (e.g., "biogrid_5").
    - adata_h5ad_path: Path to the non-curated AnnData h5ad file.
    - save_curated_h5ad: Whether to save the curated AnnData object.
    - save_curated_parquet: Whether to save the curated data as Parquet files.
    - split_parquet: Whether to save separate Parquet files for data and metadata.
    """

    # Create a CuratedDataset object from the non-curated AnnData file
    cur_data = CuratedDataset(
        obs_schema=ObsSchema,
        var_schema=VarSchema,
        exp_metadata_schema=Experiment,
        noncurated_path=adata_h5ad_path.as_posix(),
    )

    # Load the non-curated data
    cur_data.load_data(curated=False)

    # remove rows with NA in the 'perturbation_name' column
    cur_data.remove_na(slot="obs", column="perturbation_name")

    # Standardize gene symbols in the 'perturbed_target_symbol' column
    cur_data.standardize_genes(
        slot="obs",
        input_column="perturbed_target_symbol",
        input_column_type="gene_symbol",
        multiple_entries=False,
    )

    # count number of perturbations in each sample
    cur_data.count_entries(
        slot="obs",
        input_column="perturbed_target_symbol",
        count_column_name="perturbed_target_number",
        sep="|",
    )

    # add chromosome encoding
    cur_data.chromosome_encoding()

    # match the order of columns in obs to the schema
    cur_data.match_schema_columns(slot="obs")

    # validate the data against the schema
    cur_data.validate_data(slot="obs")

    # save the curated AnnData object
    if save_curated_h5ad:
        cur_data.save_curated_data_h5ad()

    # save the curated data as Parquet files
    if save_curated_parquet:
        cur_data.save_curated_data_parquet(split_metadata=split_parquet)

    return cur_data


def make_adata_biogrid(
        biogrid_dataset_id: str = None,
        biogrid_screen_path: Path = None,
        biogrid_metadata_df: pd.DataFrame = None,
        curated_metadata_dict: dict = None,
        save_h5ad_dir: str = None,
):
    """
    Create an AnnData object from a BioGRID screen file and curated metadata.

    Parameters:
    - biogrid_dataset_id: Identifier for the BioGRID dataset.
    - biogrid_screen_path: Path to the BioGRID screen file.
    - biogrid_metadata: DataFrame containing BioGRID metadata.
    - curated_metadata_dict: Dictionary containing curated metadata.
    - save_h5ad_dir: Directory to save the AnnData object as an h5ad file.
    """

    # Extract the numeric part of the dataset ID
    biogrid_dataset_id_num = int(biogrid_dataset_id.lstrip('biogrid_'))

    # Load the BioGRID screen data
    biogrid_screen_df = pd.read_csv(biogrid_screen_path, sep="\t")

    # create a boolean 'significant' column based on the 'HIT' column
    biogrid_screen_df['significant'] = biogrid_screen_df['HIT'].apply(lambda x: True if x == 'YES' else False)

    # get the column mappings for the score of interest
    biogrid_screen_df, score_columns = map_metadata_score_columns_biogrid(
        biogrid_screen_df, biogrid_metadata_df, biogrid_dataset_id_num
    )

    # convert the new columns to numeric
    biogrid_screen_df[score_columns] = biogrid_screen_df[score_columns].apply(pd.to_numeric, errors='coerce')

    # group_by 'IDENTIFIER_ID' and average the counts
    biogrid_screen_df = (
        biogrid_screen_df.groupby([e for e in biogrid_screen_df.columns if e not in score_columns])
        .agg({e: "mean" for e in score_columns})
        .reset_index()
    )

    # create unique index for the adata.obs
    biogrid_screen_df['index'] = "biogrid_" + biogrid_screen_df['#SCREEN_ID'].astype(str) + '_' + biogrid_screen_df[
        'OFFICIAL_SYMBOL']

    # remove duplicates if any
    biogrid_screen_df = biogrid_screen_df[~biogrid_screen_df['index'].duplicated()]

    # get X data for adata object based on the score columns
    X_df = biogrid_screen_df[score_columns]

    # set the index to the unique index created above
    X_df.index = biogrid_screen_df['index']

    # create adata.obs DataFrame
    OBS_df = pd.DataFrame(index=X_df.index)
    OBS_df["sample_id"] = OBS_df.index
    OBS_df["perturbation_name"] = OBS_df.index
    OBS_df["perturbed_target_symbol"] = biogrid_screen_df["OFFICIAL_SYMBOL"].to_list()
    OBS_df["significant"] = biogrid_screen_df["significant"].to_list()
    OBS_df['significance_criteria'] = biogrid_metadata_df.loc[
        biogrid_metadata_df["#SCREEN_ID"] == biogrid_screen_df['#SCREEN_ID'].values[0], 'significance_criteria'].values[
        0]
    OBS_df['guide_sequence'] = None
    OBS_df = OBS_df.drop_duplicates()

    # add all columns from gemini-curated metadata
    for key, value in curated_metadata_dict.items():
        OBS_df[key] = value

    # create adata.var DataFrame
    VAR_df = pd.DataFrame(index=score_columns)
    VAR_df["score_name"] = score_columns

    # replace None with np.nan to avoid issues with AnnData writing
    OBS_df = OBS_df.replace({None: np.nan})
    VAR_df = VAR_df.replace({None: np.nan})

    # create AnnData object
    adata = AnnData(X=X_df, obs=OBS_df, var=VAR_df)

    # save the AnnData object as an h5ad file
    if save_h5ad_dir:
        h5ad_path = Path(save_h5ad_dir) / f"{biogrid_dataset_id}.h5ad"
        adata.write_h5ad(h5ad_path)
        print(f"✅ Saved AnnData object to {h5ad_path}")

    return adata, h5ad_path