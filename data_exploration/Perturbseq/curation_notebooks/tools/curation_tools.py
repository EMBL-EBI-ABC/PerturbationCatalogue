import pandas as pd
import os
import requests
from typing import Optional, Literal
from libchebipy import search
from pathlib import Path
import scanpy as sc
from Bio import SeqIO
from io import StringIO
from gprofiler import GProfiler

# Get the path to the ontologies directory relative to this file
ONTOLOGIES_DIR = Path(__file__).parent.parent.parent / "ontologies"


# function to save a new term to ctype ontology
def add_new_term(
    ontology_df_name: Literal["ctype_ont", "cline_ont", "tis_ont", "dis_ont"],
    term_name: str,
    ontology_id: str,
    synonyms: str,
    description: str,
    save_ontology: bool = False,
):
    """
    Function to save a new term to ontology.

    Args:

    ontology_df_name: str
        Name of the ontology dataframe to update.
    term_name: str
        Name of the new term.
    ontology_id: str
        Ontology ID of the new term.
    synonyms: str
        Synonyms of the new term, separated by '|'.
    description: str
        Description of the new term.
    save_ontology: bool
        If True, save the updated ontology to disk.
    """

    global ctype_ont, cline_ont, tis_ont, dis_ont
    # load the ontology dataframe based on the name
    if ontology_df_name == "ctype_ont":
        ont = ctype_ont
    elif ontology_df_name == "cline_ont":
        ont = cline_ont
    elif ontology_df_name == "tis_ont":
        ont = tis_ont
    elif ontology_df_name == "dis_ont":
        ont = dis_ont
    else:
        raise ValueError("Invalid ontology dataframe name.")
    # check if the term name is valid

    new_term = pd.DataFrame(
        {
            "name": [term_name],
            "ontology_id": [ontology_id],
            "synonyms": [synonyms],
            "description": [description],
            "created_at": [datetime.now()],
        }
    )
    # schema for validation
    ont_schema = pa.DataFrameSchema(
        {
            "name": pa.Column(
                pa.String, nullable=False, checks=pa.Check.str_length(min_value=1)
            ),
            "ontology_id": pa.Column(
                pa.String,
                nullable=False,
                checks=[pa.Check.str_contains(":"), pa.Check.str_length(min_value=1)],
            ),
            "synonyms": pa.Column(
                pa.String, nullable=True, checks=pa.Check.str_length(min_value=1)
            ),
            "description": pa.Column(
                pa.String, nullable=True, checks=pa.Check.str_length(min_value=1)
            ),
            "created_at": pa.Column(pa.DateTime, nullable=False),
        }
    )
    # validate the new term
    try:
        validated_term = ont_schema.validate(new_term)
        print("New term is valid.")
    except pa.errors.SchemaError as e:
        print(f"Validation error: {e}")
        print(json.dumps(e.message, indent=2))

    # check if the term already exists
    if ontology_id in ont["ontology_id"].values:
        print(f"Term with ontology ID {ontology_id} already exists.")
        return
    if term_name in ont["name"].values:
        print(f"Term with name {term_name} already exists.")
        return

    # append the new term to the ontology
    ont = pd.concat([ont, validated_term], ignore_index=True)

    if save_ontology:
        # save the updated ontology

        if ontology_df_name == "ctype_ont":
            ont.to_parquet("Ontologies/cell_types.parquet", index=False)
        elif ontology_df_name == "cline_ont":
            ont.to_parquet("Ontologies/cell_lines.parquet", index=False)
        elif ontology_df_name == "tis_ont":
            ont.to_parquet("Ontologies/tissues.parquet", index=False)
        elif ontology_df_name == "dis_ont":
            ont.to_parquet("Ontologies/diseases.parquet", index=False)
        print(f"Ontology {ontology_df_name} saved with new term.")

    # update the global variable
    if ontology_df_name == "ctype_ont":
        ctype_ont = ont
    elif ontology_df_name == "cline_ont":
        cline_ont = ont
    elif ontology_df_name == "tis_ont":
        tis_ont = ont
    elif ontology_df_name == "dis_ont":
        dis_ont = ont

    print(f"New term {term_name} is added to {ontology_df_name} ontology.")


# function to remove a term from ontology
def remove_term(
    ontology_df_name: Literal["ctype_ont", "cline_ont", "tis_ont", "dis_ont"],
    term_name: str,
    ontology_id: str,
    save_ontology: bool = False,
):
    """
    Function to remove a term from ontology.

    Args:
    ontology_df_name: str
        Name of the ontology dataframe to update.
    term_name: str
        Name of the term to remove.
    ontology_id: str
        Ontology ID of the term to remove.
    save_ontology: bool
        If True, save the updated ontology to disk.
    """
    global ctype_ont, cline_ont, tis_ont, dis_ont

    # load the ontology dataframe based on the name
    if ontology_df_name == "ctype_ont":
        ont = ctype_ont
    elif ontology_df_name == "cline_ont":
        ont = cline_ont
    elif ontology_df_name == "tis_ont":
        ont = tis_ont
    elif ontology_df_name == "dis_ont":
        ont = dis_ont
    else:
        raise ValueError("Invalid ontology dataframe name.")

    # check if the term exists
    if ontology_id not in ont["ontology_id"].values:
        print(f"Term with ontology ID {ontology_id} does not exist.")
        return
    elif term_name not in ont["name"].values:
        print(f"Term with name {term_name} does not exist.")
        return

    # remove the term from the ontology
    ont = ont[ont["ontology_id"] != ontology_id]
    ont = ont[ont["name"] != term_name]
    # reset the index
    ont.reset_index(drop=True, inplace=True)

    if save_ontology:
        # save the updated ontology
        if ontology_df_name == "ctype_ont":
            ont.to_parquet("Ontologies/cell_types.parquet", index=False)
        elif ontology_df_name == "cline_ont":
            ont.to_parquet("Ontologies/cell_lines.parquet", index=False)
        elif ontology_df_name == "tis_ont":
            ont.to_parquet("Ontologies/tissues.parquet", index=False)
        elif ontology_df_name == "dis_ont":
            ont.to_parquet("Ontologies/diseases.parquet", index=False)
        print(f"Ontology {ontology_df_name} saved with removed term.")

    # update the global variable
    if ontology_df_name == "ctype_ont":
        ctype_ont = ont
    elif ontology_df_name == "cline_ont":
        cline_ont = ont
    elif ontology_df_name == "tis_ont":
        tis_ont = ont
    elif ontology_df_name == "dis_ont":
        dis_ont = ont

    print(f"Term with ID {ontology_id} is removed from {ontology_df_name} ontology.")


# function to add synonyms to existing terms in ontologies
def add_synonym(
    ontology_df_name: Literal[
        "gene_ont", "ctype_ont", "cline_ont", "tis_ont", "dis_ont"
    ],
    ontology_id: str,
    synonym: str,
    save_ontology: bool = False,
):
    """
    Function to add synonyms to existing terms in ontologies.

    Args:
    ontology_df_name: str
        Name of the ontology dataframe to update.
    ontology_id: str
        Ontology ID of the term to update.
    synonym: str
        Synonyms to add, separated by '|'.
    save_ontology: bool
        If True, save the updated ontology to disk.
    """
    global ctype_ont, gene_ont, cline_ont, tis_ont, dis_ont

    # load the ontology dataframe based on the name
    if ontology_df_name == "gene_ont":
        ontology_df = gene_ont
    elif ontology_df_name == "ctype_ont":
        ontology_df = ctype_ont
    elif ontology_df_name == "cline_ont":
        ontology_df = cline_ont
    elif ontology_df_name == "tis_ont":
        ontology_df = tis_ont
    elif ontology_df_name == "dis_ont":
        ontology_df = dis_ont
    else:
        raise ValueError("Invalid ontology dataframe name.")
    # check if the term already exists
    if ontology_id not in ontology_df["ontology_id"].values:
        print(f"Term with ontology ID {ontology_id} does not exist.")
        return

    # check if supplied synonyms already exist
    existing_synonyms = (
        ontology_df[ontology_df["ontology_id"] == ontology_id]["synonyms"]
        .str.split("|")
        .values[0]
    )
    new_synonyms = synonym.split("|")
    for new_syn in new_synonyms:
        if len(new_syn) == 0 or new_syn == "None":
            print("Empty synonym provided. Skipping.")
            continue
        if new_syn in existing_synonyms:
            print(f"Synonym {new_syn} already exists for term {ontology_id}.")
        else:
            # add the new synonym
            existing_synonyms.append(new_syn)
            print(f"Synonym {new_syn} added for term {ontology_id}.")
    # update the synonyms column in the dataframe
    ontology_df.loc[ontology_df["ontology_id"] == ontology_id, "synonyms"] = "|".join(
        existing_synonyms
    )

    if ontology_df_name == "gene_ont":
        gene_ont = ontology_df
    elif ontology_df_name == "ctype_ont":
        ctype_ont = ontology_df
    elif ontology_df_name == "cline_ont":
        cline_ont = ontology_df
    elif ontology_df_name == "tis_ont":
        tis_ont = ontology_df
    elif ontology_df_name == "dis_ont":
        dis_ont = ontology_df

    if save_ontology:
        # save the updated ontology
        if ontology_df_name == "gene_ont":
            gene_ont.to_parquet("Ontologies/genes.parquet", index=False)
        elif ontology_df_name == "ctype_ont":
            ctype_ont.to_parquet("Ontologies/cell_types.parquet", index=False)
        elif ontology_df_name == "cline_ont":
            cline_ont.to_parquet("Ontologies/cell_lines.parquet", index=False)
        elif ontology_df_name == "tis_ont":
            tis_ont.to_parquet("Ontologies/tissues.parquet", index=False)
        elif ontology_df_name == "dis_ont":
            dis_ont.to_parquet("Ontologies/diseases.parquet", index=False)
        print(f"Ontology {ontology_df_name} saved with new synonym.")
    else:
        print(f"Ontology {ontology_df_name} not saved. Set save_ontology=True to save.")


# function to remove synonyms from existing terms in ontologies
def remove_synonym(
    ontology_df_name: Literal[
        "gene_ont", "ctype_ont", "cline_ont", "tis_ont", "dis_ont"
    ],
    ontology_id: str,
    synonym: str,
    save_ontology: bool = False,
):
    """
    Function to remove synonyms from existing terms in ontologies.

    Args:
    ontology_df_name: str
        Name of the ontology dataframe to update.
    ontology_id: str
        Ontology ID of the term to update.
    synonym: str
        Synonyms to remove, separated by '|'.
    save_ontology: bool
        If True, save the updated ontology to disk.
    """
    global ctype_ont, gene_ont, cline_ont, tis_ont, dis_ont

    # load the ontology dataframe based on the name
    if ontology_df_name == "gene_ont":
        ontology_df = gene_ont
    elif ontology_df_name == "ctype_ont":
        ontology_df = ctype_ont
    elif ontology_df_name == "cline_ont":
        ontology_df = cline_ont
    elif ontology_df_name == "tis_ont":
        ontology_df = tis_ont
    elif ontology_df_name == "dis_ont":
        ontology_df = dis_ont
    else:
        raise ValueError("Invalid ontology dataframe name.")

    # check if the term already exists
    if ontology_id not in ontology_df["ontology_id"].values:
        print(f"Term with ontology ID {ontology_id} does not exist.")
        return

    # check if supplied synonyms already exist
    existing_synonyms = (
        ontology_df[ontology_df["ontology_id"] == ontology_id]["synonyms"]
        .str.split("|")
        .values[0]
    )

    new_synonyms = synonym.split("|")

    for new_syn in new_synonyms:
        if new_syn not in existing_synonyms:
            print(f"Synonym {new_syn} does not exist for term {ontology_id}.")
        else:
            # remove the synonym
            existing_synonyms.remove(new_syn)
            print(f"Synonym {new_syn} removed for term {ontology_id}.")

    # update the synonyms column in the dataframe
    if len(existing_synonyms) > 0:
        ontology_df.loc[ontology_df["ontology_id"] == ontology_id, "synonyms"] = (
            "|".join(existing_synonyms)
        )
        print(f"New synonyms for term {ontology_id}: {existing_synonyms}")
    else:
        ontology_df.loc[ontology_df["ontology_id"] == ontology_id, "synonyms"] = None
        print(f"All synonyms removed for term {ontology_id}.")

    # update the ontology dataframe
    if ontology_df_name == "gene_ont":
        gene_ont = ontology_df
    elif ontology_df_name == "ctype_ont":
        ctype_ont = ontology_df
    elif ontology_df_name == "cline_ont":
        cline_ont = ontology_df
    elif ontology_df_name == "tis_ont":
        tis_ont = ontology_df
    elif ontology_df_name == "dis_ont":
        dis_ont = ontology_df

    if save_ontology:
        # save the updated ontology
        if ontology_df_name == "gene_ont":
            gene_ont.to_parquet("Ontologies/genes.parquet", index=False)
        elif ontology_df_name == "ctype_ont":
            ctype_ont.to_parquet("Ontologies/cell_types.parquet", index=False)
        elif ontology_df_name == "cline_ont":
            cline_ont.to_parquet("Ontologies/cell_lines.parquet", index=False)
        elif ontology_df_name == "tis_ont":
            tis_ont.to_parquet("Ontologies/tissues.parquet", index=False)
        elif ontology_df_name == "dis_ont":
            dis_ont.to_parquet("Ontologies/diseases.parquet", index=False)
        print(f"Ontology {ontology_df_name} saved with removed synonym.")


# Function to standardize the data
def standardize_data(
    obs_df: pd.DataFrame,
    obs_column: str,
    ref_df: pd.DataFrame,
    ref_column: Literal["name", "symbol"],
    return_fuzzy: bool = False,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Standardize the data in the specified column of the obs DataFrame.
    Args:
        obs_df (pd.DataFrame): The obs DataFrame containing the data to be standardized.
        obs_column (str): The name of the column in the obs DataFrame to be standardized.
        ref_df (pd.DataFrame): The ontology DataFrame containing the reference data.
        ref_column (str): The name of the column in the ontology DataFrame to be used as a reference.
    Returns:
        pd.DataFrame: The standardized data.
    """
    # from rapidfuzz import fuzz, process

    # Check if the obs_column exists in the obs DataFrame
    if obs_column not in obs_df.columns:
        raise ValueError(f"Column '{obs_column}' not found in the obs DataFrame.")
    # Check if the reference_column exists in the reference DataFrame
    if ref_column not in ref_df.columns:
        raise ValueError(f"Column '{ref_column}' not found in the reference DataFrame.")

    # explode the synonyms column in the reference DataFrame
    # to create a new row for each synonym
    exploded_ref_df = ref_df.copy()
    exploded_ref_df["ref_synonyms"] = exploded_ref_df["synonyms"].str.split("|")
    exploded_ref_df = exploded_ref_df.explode("ref_synonyms")
    exploded_ref_df["ref_synonyms_lower"] = exploded_ref_df["ref_synonyms"].str.lower()

    standardized_values = {}
    possible_match = {}
    unstandardized_values = []

    # filter out the values that are perfectly matched
    unique_values_obs = (
        obs_df[~obs_df[obs_column].isin(ref_df[ref_column].values)][obs_column]
        .dropna()
        .unique()
    )

    # create a new DataFrame with the standardized values
    standardized_df = obs_df.copy()

    # "name" is for cell types, cell lines, tissues, diseases
    if ref_column == "name":

        exploded_ref_df["name_lower"] = exploded_ref_df[ref_column].str.lower()
        exploded_ref_df["ontology_id_lower"] = exploded_ref_df[
            "ontology_id"
        ].str.lower()

        for value in unique_values_obs:
            value_lower = value.lower()

            # if value_lower is a match in the reference DataFrame, get the corresponding name
            if value_lower in exploded_ref_df["name_lower"].values:
                matched_values = exploded_ref_df[
                    exploded_ref_df["name_lower"] == value_lower
                ]
                # display(matched_values)
                replacement_value = matched_values[ref_column].values[0]
                standardized_values[value] = replacement_value
                if verbose:
                    print(
                        f"Value '{value}' is a match for '{replacement_value}'. Standardizing to '{replacement_value}'."
                    )
                continue

            # if value is a match in the reference DataFrame synonyms, get the corresponding name
            elif value_lower in exploded_ref_df["ref_synonyms_lower"].values:
                matched_values = exploded_ref_df[
                    exploded_ref_df["ref_synonyms_lower"] == value_lower
                ]
                # display(matched_values)
                replacement_value = matched_values[ref_column].values[0]
                standardized_values[value] = replacement_value
                if verbose:
                    print(
                        f"Value '{value}' is a synonym of '{replacement_value}'. Standardizing to '{replacement_value}'."
                    )
                continue

            else:
                # use thefuzz.process to get the best match
                best_match_name = process.extract(
                    value_lower, exploded_ref_df["name_lower"], scorer=fuzz.ratio
                )
                best_match_name_df = pd.DataFrame(
                    best_match_name, columns=["match", "score", "index"]
                )

                best_match_name_df["source"] = "name"
                best_match_synonym = process.extract(
                    value_lower,
                    exploded_ref_df["ref_synonyms_lower"],
                    scorer=fuzz.ratio,
                )
                best_match_synonym_df = pd.DataFrame(
                    best_match_synonym, columns=["match", "score", "index"]
                )
                best_match_synonym_df["source"] = "synonym"
                best_match = pd.concat(
                    [best_match_name_df, best_match_synonym_df], axis=0
                )
                best_match = best_match.sort_values(
                    by="score", ascending=False
                ).drop_duplicates(subset=["match"])
                possible_match_values = list(best_match["match"].values)
                possible_match[value] = possible_match_values
                if verbose:
                    print(
                        f"Value '{value}' is a possible match for {possible_match_values}:"
                    )
                    display(best_match)

                unstandardized_values.append(value)

                continue

        # map the standardized values to the obs DataFrame
        standardized_df[obs_column] = standardized_df[obs_column].map(
            standardized_values
        )

    # "symbol" is for genes
    elif ref_column == "symbol":

        exploded_ref_df = exploded_ref_df.rename(
            columns={
                "symbol": "ref_symbol",
            }
        )

        exploded_ref_df["ref_symbol_lower"] = exploded_ref_df["ref_symbol"].str.lower()

        # vectorised lowercase matching
        standardized_df[f"{obs_column}_lower"] = standardized_df[obs_column].str.lower()

        # check for perfect matches
        standardized_df = standardized_df.merge(
            exploded_ref_df[["ref_symbol_lower", "ref_symbol"]].drop_duplicates(),
            how="left",
            left_on=f"{obs_column}_lower",
            right_on="ref_symbol_lower",
        ).set_axis(standardized_df.index)

        print(
            f"Number of perfect matches: {standardized_df['ref_symbol_lower'].notna().sum()}"
        )

        standardized_df.loc[
            (standardized_df[obs_column].notna())
            & (standardized_df["ref_symbol"].notna()),
            obs_column,
        ] = standardized_df.loc[
            (standardized_df[obs_column].notna())
            & (standardized_df["ref_symbol"].notna()),
            "ref_symbol",
        ]

        standardized_df = standardized_df.drop("ref_symbol", axis=1)

        # check for synonym matches
        standardized_df = standardized_df.merge(
            exploded_ref_df[["ref_synonyms_lower", "ref_symbol"]].drop_duplicates(
                subset=["ref_synonyms_lower"]
            ),
            how="left",
            left_on=f"{obs_column}_lower",
            right_on="ref_synonyms_lower",
        ).set_axis(standardized_df.index)

        print(
            f"Number of synonym matches: {standardized_df['ref_symbol'].notna().sum()}"
        )

        standardized_df.loc[
            (standardized_df[obs_column].notna())
            & (standardized_df["ref_symbol"].notna()),
            obs_column,
        ] = standardized_df.loc[
            (standardized_df[obs_column].notna())
            & (standardized_df["ref_symbol"].notna()),
            "ref_symbol",
        ]

        standardized_df = standardized_df.drop("ref_symbol", axis=1)

        if return_fuzzy:
            unique_ref_symbols = (
                exploded_ref_df["ref_symbol"].dropna().drop_duplicates()
            )

            standardized_df.loc[
                (
                    ~standardized_df[obs_column].isin(
                        exploded_ref_df["ref_symbol"].dropna().drop_duplicates()
                    )
                )
                & (
                    ~standardized_df["ensembl_gene_id"].isin(
                        exploded_ref_df["ensembl_gene_id"].dropna().drop_duplicates()
                    )
                ),
                "fuzz_name",
            ] = standardized_df.loc[
                (
                    ~standardized_df[obs_column].isin(
                        exploded_ref_df["ref_symbol"].dropna().drop_duplicates()
                    )
                )
                & (
                    ~standardized_df["ensembl_gene_id"].isin(
                        exploded_ref_df["ensembl_gene_id"].dropna().drop_duplicates()
                    )
                ),
                obs_column,
            ].apply(
                lambda x: (
                    process.extract(x, choices=unique_ref_symbols, limit=1)
                    if pd.notna(x)
                    else None
                )[0]
            )

            standardized_df["fuzz_gene"] = standardized_df["fuzz_name"].apply(
                lambda x: x[0] if isinstance(x, tuple) else None
            )
            standardized_df["fuzz_score"] = standardized_df["fuzz_name"].apply(
                lambda x: x[1] if isinstance(x, tuple) else None
            )

            standardized_df.drop(columns=["fuzz_name"], inplace=True)

            print(
                f"Number of fuzzy matches: {standardized_df['fuzz_gene'].notna().sum()}"
            )

        standardized_df.drop(
            columns=[f"{obs_column}_lower", "ref_symbol_lower", "ref_synonyms_lower"],
            inplace=True,
        )
        # standardized_df = standardized_df.set_index('ensembl_gene_id', drop=False)

    if standardized_values:
        print("Standardized values:")
        print(standardized_values)
    if possible_match:
        print("Possible matches (fix these manually):")
        print(possible_match)
    if unstandardized_values:
        print("Unstandardized values:")
        print(unstandardized_values)

    return standardized_df


def get_vals(series, dtype):
    """
    Function to get the values of a series and convert them to a specific data type.
    Args:
        series: pd.Series
            The series to get the values from.
        dtype: str
            The data type to convert the values to.
    """
    series = series.dropna().unique()

    if len(series) == 0:
        return None

    if dtype == "list":
        series = [str(x) for x in series]
        # if len(series) == 0:
        #     series = None
    elif dtype == "str":
        series = str(series[0])
    elif dtype == "int":
        series = int(series[0])
    elif dtype == "float":
        series = float(series[0])

    return series


def get_dict_vals(term_id, term_label, adata):
    """
    Function to get the values from adata object and convert them to a list with dictionaries.
    Args:
        term_id: str
            The name of the term ID column.
        term_label: str
            The name of the term label column.
        adata: AnnData object
            The AnnData object to get the values from.
    """
    # get the values from adata object
    df = adata.obs[[term_id, term_label]].dropna().drop_duplicates()
    # rename the columns
    df = df.rename(columns={term_id: "term_id", term_label: "term_label"})
    # convert the values to a list of dictionaries
    dict_vals = []
    for index, row in df.iterrows():
        dict_vals.append({"term_id": row["term_id"], "term_label": row["term_label"]})
    
    if len(dict_vals) == 0:
        return None

    return dict_vals


def search_compounds_in_chebi(compound_names):
    """
    Search for compound names in ChEBI and return a DataFrame with the results.

    Parameters:
        compound_names (str or list of str): A single compound name or a list of compound names.

    Returns:
        pd.DataFrame: A DataFrame with columns ['original_name', 'standardized_name', 'chebi_id'].
    """
    # Ensure compound_names is a list
    if isinstance(compound_names, str):
        compound_names = [compound_names]

    # Initialize a list to store the search results
    search_results = []

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
                    search_results.append({
                        "original_name": compound_name,  # Original input compound name
                        "standardized_name": chebi_entity.get_name(),  # Standardized compound name
                        "chebi_id": chebi_entity.get_id()  # ChEBI ID
                    })
                    break  # Stop after finding the first match for the compound

    # Convert the results list into a pandas DataFrame
    return pd.DataFrame(search_results, columns=["original_name", "standardized_name", "chebi_id"])


def standardize_gene_symbols(obs_df: pd.DataFrame, column: str) -> pd.DataFrame:
    """
    Standardize gene symbols using genes.parquet ontology.

    Args:
        obs_df: Input DataFrame from adata.obs or adata.var slots containing gene symbols
        column: Name of column containing gene symbols

    Returns:
        DataFrame with standardized values and mapping info
    """
    try:
        # Load gene ontology
        gene_ont = pd.read_parquet(ONTOLOGIES_DIR / "genes.parquet")
        print(f"Loaded gene ontology with {len(gene_ont)} entries\n{'-' * 50}")

        # add a column with lowercase gene symbols for matching
        # and rename the original symbol column to standardized_symbol
        gene_ont = gene_ont.assign(
            lower_symbol=lambda x: x["symbol"].str.lower()
        ).rename(columns={"symbol": "standardized_symbol"})

        # Create a DataFrame with unique perturbed target symbols
        # and their lowercase versions
        result_df = (
            obs_df[[column]]
            .reset_index(drop=True)
            .drop_duplicates()
            .assign(lower_symbol=lambda x: x[column].str.lower())
        )

        # Left-merge with gene ontology to find direct matches with standardized symbols
        result_df = result_df.merge(
            gene_ont[["lower_symbol", "standardized_symbol"]],
            how="left",
            left_on="lower_symbol",
            right_on="lower_symbol",
            indicator=True,
        ).drop_duplicates(subset="lower_symbol")

        # Separate matched and unmatched gene symbols
        result_df_matched = result_df[result_df["_merge"] == "both"].drop(
            columns=["_merge"]
        )
        result_df_unmatched = result_df[result_df["_merge"] != "both"].drop(
            columns=["_merge", "standardized_symbol"]
        )

        print(
            f"{len(result_df_matched)} out of {len(result_df)} gene symbols mapped to standardized symbols\n{'-' * 50}"
        )
        print(
            f"{len(result_df_unmatched)} gene symbols could not be mapped to standardized symbols\n{'-' * 50}"
        )

        # If there are unmatched gene symbols, try to match them against the known synonyms
        if len(result_df_unmatched) > 0:

            print(f"Trying to match the unmatched gene symbols against known synonyms\n{'-' * 50}")

            # Explode the synonyms column in the gene ontology
            # and create a DataFrame with lowercase synonyms
            # and their corresponding standardized symbols
            exploded_gene_ont = (
                gene_ont.assign(synonyms=lambda x: x["synonyms"].str.split("|"))
                .explode("synonyms")
                .assign(lower_synonyms=lambda x: x["synonyms"].str.lower())[
                    ["lower_synonyms", "standardized_symbol"]
                ]
                .drop_duplicates()
            )

            # Left-merge the unmatched gene symbols with the exploded gene ontology
            # to find matches based on synonyms
            result_df_synonyms = result_df_unmatched.merge(
                exploded_gene_ont,
                how="left",
                left_on="lower_symbol",
                right_on="lower_synonyms",
                indicator=True,
            ).drop_duplicates(subset="lower_symbol")

            # Separate matched and unmatched gene symbols
            result_df_synonyms_matched = result_df_synonyms[
                result_df_synonyms["_merge"] == "both"
            ].drop(columns=["_merge", "lower_synonyms"])

            # update the unmatched dataframe
            result_df_unmatched = result_df_synonyms[
                result_df_synonyms["_merge"] != "both"
            ].drop(columns=["_merge", "lower_synonyms"])

            print(
                f"{len(result_df_synonyms_matched)} gene symbols mapped to standardized symbols using synonyms\n{'-' * 50}"
            )

            # Add the matched synonyms to the result DataFrame
            result_df_matched = pd.concat(
                [result_df_matched, result_df_synonyms_matched], ignore_index=True
            )

        # If there are still unmatched gene symbols, try to match them against the ENA identifiers using gprofiler
        if len(result_df_unmatched) > 0:
            print(f"Trying to match the unmatched gene symbols against ENA identifiers\n{'-' * 50}")

            gp = GProfiler(return_dataframe=True)
            converted = gp.convert(organism='hsapiens',
                        query=result_df_unmatched['gene_symbol'].to_list(),
                        target_namespace='ENSG')

            converted = converted[converted['n_converted'] == 1]
            converted = converted[converted['name'] != 'None']
            converted = converted[["incoming", "name"]].rename(
                columns={"incoming": "gene_symbol", "name": "standardized_symbol"}
            )
            
            result_df_converted = result_df_unmatched.drop(columns=['standardized_symbol']).merge(
                converted,
                how="left",
                left_on="gene_symbol",
                right_on="gene_symbol",
                indicator=True,
            ).drop_duplicates(subset="lower_symbol")#.drop(columns=['lower_symbol'])
            
            # Separate matched and unmatched gene symbols
            result_df_converted_matched = result_df_converted[
                result_df_converted["_merge"] == "both"
            ].drop(columns=["_merge"])
            # update the unmatched dataframe
            result_df_unmatched = result_df_converted[
                result_df_converted["_merge"] != "both"
            ].drop(columns=["_merge"])
            print(
                f"{len(result_df_converted_matched)} gene symbols mapped to standardized symbols using ENA identifiers\n{'-' * 50}"
            )
            # Add the matched synonyms to the result DataFrame
            result_df_matched = pd.concat(
                [result_df_matched, result_df_converted_matched], ignore_index=True
            )

        # If there are still unmatched gene symbols, keep them as is
        # and assign them to the standardized_symbol column
        if len(result_df_unmatched) > 0:
            print(
                f"{len(result_df_unmatched)} gene symbols could not be mapped to standardized symbols using ENA identifiers\n{'-' * 50}"
            )
            print("Unmatched genes will be kept as is in the final DataFrame")
            print("Unmatched gene symbols:", result_df_unmatched[column].unique())

            # Assign the unmatched gene symbols to the standardized_symbol column
            result_df_unmatched = result_df_unmatched.assign(
                standardized_symbol=result_df_unmatched[column]
            )

            # Add the unmatched gene symbols to the result DataFrame
            result_df_matched = pd.concat(
                [result_df_matched, result_df_unmatched], ignore_index=True
            )
        else:
            print(
                f"All unmatched gene symbols have been mapped to standardized symbols using synonyms\n{'-' * 50}"
            )

        result_df_matched = result_df_matched.drop(columns=["lower_symbol"])

        display(result_df_matched)

        # Map the standardized symbols back to the original DataFrame
        obs_df = (
            obs_df.reset_index(drop=True).merge(
                result_df_matched,
                how="left",
                left_on=column,
                right_on=column
            )
            .set_index(column)
            # .drop(columns=[column])
            .rename(columns={"standardized_symbol": column})
        )

        print(
            f"Mapped the standardized symbols in column {column} back to the original DataFrame"
        )

        return obs_df

    except Exception as e:
        print(f"Error standardizing gene symbols: {str(e)}")
        raise


def standardize_ontology(
    obs_df: pd.DataFrame,
    column: str,
    ont_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Standardize values in a DataFrame column using an ontology DataFrame.

    Args:
        obs_df: Input DataFrame containing values to be standardized
        column: Name of column to be standardized
        ont_df: Ontology DataFrame containing standardized values
    Returns:
        DataFrame with standardized values
    """
    
    try:
        # Rename the original name column to standardized_name
        ont_df = ont_df.rename(columns={"name": "standardized_name"}).assign(
            lower_standardized_name=lambda x: x["standardized_name"].str.lower()
            )
        
        # Make a version of the ontology DataFrame with pluralized names
        ont_df_plural = ont_df.assign(
            lower_standardized_name=lambda x: x["standardized_name"].str.lower()+"s"
        )
        
        # Concatenate the original and pluralized DataFrames
        ont_df = pd.concat([ont_df, ont_df_plural], ignore_index=True)

        # Create a DataFrame with unique ontology labels
        # and their lowercase versions
        result_df = (
            obs_df[[column]]
            .reset_index(drop=True)
            .drop_duplicates()
            .assign(lower_name=lambda x: x[column].str.lower())
        )

        # Left-merge with gene ontology to find direct matches with standardized names
        result_df = result_df.merge(
            ont_df[["lower_standardized_name", "standardized_name"]],
            how="left",
            left_on="lower_name",
            right_on="lower_standardized_name",
            indicator=True,
        ).drop_duplicates(subset="lower_name")

        # Separate matched and unmatched ontology labels
        result_df_matched = result_df[result_df["_merge"] == "both"].drop(
            columns=["_merge", "lower_standardized_name"]
        )
        result_df_unmatched = result_df[result_df["_merge"] != "both"].drop(
            columns=["_merge", "standardized_name", "lower_standardized_name"]
        )

        print(
            f"{len(result_df_matched)} out of {len(result_df)} ontology labels mapped to standardized names\n{'-' * 50}"
        )
        print(
            f"{len(result_df_unmatched)} ontology label could not be mapped to standardized names\n{'-' * 50}"
        )

        # If there are unmatched ontology labels, try to match them against the known synonyms
        if len(result_df_unmatched) > 0:

            print(f"Trying to match the unmatched ontology labels against known synonyms\n{'-' * 50}")

            # Explode the synonyms column in the gene ontology
            # and create a DataFrame with lowercase synonyms
            # and their corresponding standardized names
            exploded_ont = (
                ont_df.dropna(subset=["synonyms"])
                .assign(synonyms=lambda x: x["synonyms"].str.split("|"))
                .explode("synonyms")
                .assign(lower_synonyms=lambda x: x["synonyms"].str.lower())[
                    ["lower_synonyms", "standardized_name"]
                ]
                .drop_duplicates()
            )
            
            # Make a version of the exploded DataFrame with pluralized names
            exploded_ont_plural = (
                exploded_ont.dropna(subset=["lower_synonyms"])
                .assign(lower_synonyms=lambda x: x["lower_synonyms"].str.lower()+"s")
            )
            # Concatenate the original and pluralized DataFrames
            exploded_ont = pd.concat([exploded_ont, exploded_ont_plural], ignore_index=True)

            # Left-merge the unmatched ontology labels with the exploded gene ontology
            # to find matches based on synonyms
            result_df_synonyms = result_df_unmatched.merge(
                exploded_ont,
                how="left",
                left_on="lower_name",
                right_on="lower_synonyms",
                indicator=True,
            ).drop_duplicates(subset="lower_name")

            # Separate matched and unmatched ontology labels
            result_df_synonyms_matched = result_df_synonyms[
                result_df_synonyms["_merge"] == "both"
            ].drop(columns=["_merge", "lower_synonyms"])

            # update the unmatched dataframe
            result_df_unmatched = result_df_synonyms[
                result_df_synonyms["_merge"] != "both"
            ].drop(columns=["_merge", "lower_synonyms"])

            print(
                f"{len(result_df_synonyms_matched)} ontology label mapped to standardized names using synonyms\n{'-' * 50}"
            )

            # Add the matched synonyms to the result DataFrame
            result_df_matched = pd.concat(
                [result_df_matched, result_df_synonyms_matched], ignore_index=True
            )

            # If there are still unmatched ontology labels, keep them as is
            # and assign them to the standardized_name column
            if len(result_df_unmatched) > 0:
                print(
                    f"{len(result_df_unmatched)} ontology labels could not be mapped to standardized names using synonyms\n{'-' * 50}"
                )
                print("These ontology labels will be kept as is in the final DataFrame")
                print(f"Unmatched ontology labels: {list(result_df_unmatched[column].unique())}")

                # Assign the unmatched ontology label to the standardized_name column
                result_df_unmatched = result_df_unmatched.assign(
                    standardized_name=result_df_unmatched[column]
                ).drop(columns=["lower_name"])

                # Add the unmatched gene symbols to the result DataFrame
                result_df_matched = pd.concat(
                    [result_df_matched, result_df_unmatched], ignore_index=True
                )
            else:
                print(
                    f"All unmatched ontology labels have been mapped to standardized names using synonyms\n{'-' * 50}"
                )

        result_df_matched = result_df_matched.drop(columns=["lower_name"])

        display(result_df_matched)

        # Map the standardized symbols back to the original DataFrame
        obs_df = (
            obs_df.merge(
                result_df_matched,
                how="left",
                left_on=column,
                right_on=column,
            )
            .set_index(obs_df.index)
            .drop(columns=[column])
            .rename(columns={"standardized_name": column})
        )

        print(
            f"Mapped the standardized ontology labels in column {column} back to the original DataFrame"
        )

        return obs_df

    except Exception as e:
        print(f"Error standardizing ontology terms: {str(e)}")
        raise

def remove_version_from_genes(df, column, sep = "."):
    """
    Remove version numbers from gene symbols or ENSG IDs in a DataFrame column.
    Args:
        df: Input DataFrame containing gene symbols/ENSG IDs
        column: Name of column containing gene symbols/ENSG IDs
        sep: Separator used between the gene symbols/ENSG IDs and the version (default is ".")
    Returns:
        DataFrame with version numbers removed from gene symbols/ENSG IDs
    """
    # Check if the column exists in the DataFrame
    if column not in df.columns:
        raise ValueError(f"Column {column} not found in DataFrame")
    # Check if the column is empty
    if df[column].empty:
        raise ValueError(f"Column {column} is empty")

    # Remove version numbers from gene symbols/ENSG IDs
    df[column] = df[column].str.split(sep).str[0]

    return df


def standardize_var_genes(df, column = Literal['gene_symbol', 'ensembl_gene_id'], remove_version = False, sep = '.'):
    """
    Standardize gene symbols or ENSG in a DataFrame column using gprofiler.
    Args:
        obs_df: Input DataFrame containing gene symbol/ENSG IDs
        column: Name of column containing gene symbols/ENSG IDs
        remove_version: Boolean indicating whether to remove version numbers from gene symbols/ENSG IDs (default is False)
        sep: Separator used between the gene symbols/ENSG IDs and the version (default is ".")
    Returns:
        DataFrame with standardized gene symbols and ENSG IDs
        
    Example:
        df = standardize_genes(df, 'gene_symbol', remove_version=True, sep = '.')
        df = standardize_genes(df, 'ensembl_gene_id')
    """

    # Check if the column is gene symbol or ENSG
    if column not in ["gene_symbol", "ensembl_gene_id"]:
        raise ValueError("Column must be either 'gene_symbol' or 'ensembl_gene_id'")
    # Check if the column exists in the DataFrame
    if column not in df.columns:
        raise ValueError(f"Column {column} not found in DataFrame")
    # Check if the column is empty
    if df[column].empty:
        raise ValueError(f"Column {column} is empty")

    # keep original query column (will be added back at the end of the processing)
    orig_column = df[column]
    orig_column_name = f"original_{column}"
    # rename df index to avoid naming conflicts
    df.index = df.index.rename('index')

    # Remove version numbers from gene symbols/ENSG IDs
    if remove_version:
        df = remove_version_from_genes(df, column, sep)

    # Convert the gene symbols/ENSG IDs to the target namespace (ENSG)
    gp = GProfiler(return_dataframe=True)
    converted = gp.convert(organism='hsapiens',
                            query=df[column].to_list(),
                            target_namespace='ENSG')

    # Filter the converted DataFrame to keep only the rows with a single conversion
    converted = (
        converted[converted["n_converted"] == 1]
        .drop(columns=["n_incoming", "n_converted"])
    )

    matched = converted[converted['converted'] != 'None']

    print(f"Converted {len(matched)}/{len(converted)} gene symbols/ENSG IDs to standardized gene symbols/ENSG IDs\n{'-' * 50}")
    
    df[f"original_{column}"] = orig_column

    # merge the converted DataFrame with the original DataFrame
    converted = (
        df.merge(
            converted[["incoming", "converted", "name"]],
            how="left",
            left_on=column,
            right_on="incoming",
            indicator=True,
        ).drop(columns=[column])
        .drop_duplicates(subset=[orig_column_name])
    )

    # ensure the length of the converted DataFrame is the same as the original DataFrame
    if len(converted) != len(df):
        raise ValueError(f"Length of converted DataFrame ({len(converted)}) does not match length of original DataFrame ({len(df)})")

    # add the original index
    converted.index = df.index

    # rename the columns
    new_colnames_map = {"converted": "ensembl_gene_id", "name": "gene_symbol"}
    converted = converted.rename(columns=new_colnames_map)

    # keep only the relevant columns
    converted = converted[list(new_colnames_map.values()) + [orig_column_name]]
    
    # replace "None" with None
    converted = converted.replace("None", None)

    return converted
