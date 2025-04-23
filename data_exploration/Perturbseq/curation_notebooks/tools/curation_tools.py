import pandas as pd
import os
from neofuzz import char_ngram_process, Process
from typing import Optional, Literal



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
