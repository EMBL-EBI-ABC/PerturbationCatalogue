import os
from pathlib import Path
import numpy as np
import pandas as pd
import json
from pprint import pprint

from pydantic import ValidationError
from typing import Literal
import pandera as pa

from libchebipy import search
import scanpy as sc
from gprofiler import GProfiler


class CuratedDataset:

    # Get the path to the ontologies directory relative to this file
    ONTOLOGIES_DIR = Path(__file__).parent.parent.parent / "ontologies"

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
            The path to the non-curated data. The default is None.
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
        self.curated_path = noncurated_path.replace("non_curated", "curated").replace(
            ".h5ad", "_curated.h5ad"
        )

        # Initialise adata object
        self.adata = None
        # Initialise a template for experiment metadata
        self.exp_metadata = {
            k: {} for k in self.exp_metadata_schema.model_fields.keys()
        }

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

    def load_data(self):
        """
        Load the adata from the specified source.
        """
        if os.path.exists(self.noncurated_path):
            print(f"Loading data from {self.noncurated_path}")
            self.adata = sc.read_h5ad(self.noncurated_path)
            for col in self.adata.obs.columns:
                if self.adata.obs[col].dtype == "object":
                    self.adata.obs[col] = self.adata.obs[col].astype(
                        "str"
                    )  # .astype("category")
        else:
            print(
                f"File {self.noncurated_path} does not exist. Run download_data() first."
            )

    def save_curated_data(self):
        """Save the curated data to a h5ad file."""

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
        print(f"Curated data saved to {self.curated_path}")

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
        to_replace=None,
        replace_value=None,
    ):
        """
        Replace entries in a column of the named slot of adata object.
        Parameters
        ----------
        slot : str
            The slot to replace entries in. Can be either "obs" or "var".
        column : str
            The name of the column to replace entries in.
        to_replace : str
            The value to replace. Must be a regex-like string.
        replace_value : str
            The value to replace with.
        """
        if slot not in ["obs", "var"]:
            raise ValueError('slot must be either "obs" or "var"')
        df = getattr(self.adata, slot)
        if column not in df.columns:
            raise ValueError(f"Column {column} not found in adata.{slot}")
        if df[column].empty:
            raise ValueError(f"Column {column} is empty in adata.{slot}")

        df[column] = df[column].str.replace(to_replace, replace_value, regex=True)

        setattr(self.adata, slot, df)

        print(
            f"Replaced entries {to_replace} -> {replace_value} in column {column} of adata.{slot}"
        )

    def map_values_from_column(self, ref_col, target_col, ref_value, target_value):
        """
        Replace values in target_col based on corresponding values in ref_col using ref_value and target_value.
        """
        
        df = self.adata.obs
        
        if ref_col not in df.columns:
            raise ValueError(f"Column {ref_col} not found in adata.obs")
        if target_col not in df.columns:
            raise ValueError(f"Column {target_col} not found in adata.obs")

        if df[ref_col].empty:
            raise ValueError(f"Column {ref_col} is empty in adata.obs")
        if df[target_col].empty:
            raise ValueError(f"Column {target_col} is empty in adata.obs")
        
        # Check if the ref_value exists in the ref_col
        if not (df[ref_col] == ref_value).any():
            raise ValueError(
                f"Value {ref_value} not found in column {ref_col} of adata.obs"
            )
        
        # convert the target_col to string if it is not already
        df[target_col] = df[target_col].astype(str)
            
        # Replace the ref_value with target_value in target_col where ref_col matches ref_value
        df.loc[df[ref_col] == ref_value, target_col] = target_value
        print(
            f"Replaced '{ref_value}' -> '{target_value}' in column '{target_col}' based on values in '{ref_col}'"
        )
        
    
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
        else:
            raise ValueError(f"Column {column} has no NA entries in adata.{slot}")

        print(
            f"Removed {sum(na_entries)} NA entries from column {column} of adata.{slot}"
        )

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

    def standardize_compounds(self, column=None):
        """
        Standardize compound names in a DataFrame column using ChEBI.

        Parameters:
            column (str): The name of the column containing compound names to be standardized.
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
            query=conv_df[input_column].to_list(),
            target_namespace="ENSG",
        )

        # Filter the converted DataFrame to keep only the rows with a single conversion
        gp_result = gp_result[gp_result["n_converted"] == 1].drop(
            columns=["n_incoming", "n_converted"]
        )

        gp_result.loc[gp_result["incoming"] == "control", "converted"] = "control"
        gp_result.loc[gp_result["incoming"] == "control", "name"] = "control"

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

        # add the biotype column
        conv_df = (
            conv_df.merge(
                self.gene_ont[["ensembl_gene_id", "biotype"]]
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
            raise ValueError(
                f"No {ontology_type} ontology terms could be mapped from `{input_column}` column to ontology terms. Map the terms manually or check the input column for errors."
            )

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
                    df = df.drop(columns=col)

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
                    self._get_vals("perturbed_target_ensg")
                ),
                "perturbed_targets": self._get_vals("perturbed_target_ensg"),
                "number_of_perturbed_cells": self.adata.obs.shape[0],
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

        schema = self.obs_schema if slot == "obs" else self.var_schema

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

        exploded_cols = ["incoming", "converted", "name", "biotype"]

        df = df.groupby([unique_val_column]).agg(
            {
                col: lambda x: sep.join(x) if x.notna().all() else None
                for col in exploded_cols
            }
        )

        print(f"Collapsed column {unique_val_column} using separator {sep}")

        return df
