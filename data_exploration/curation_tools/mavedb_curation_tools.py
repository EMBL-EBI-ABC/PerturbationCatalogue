import json
import anndata
import numpy as np
import pandas as pd
from pathlib import Path
import re

from curation_tools.curation_tools import (
    CuratedDataset,
    ObsSchema,
    VarSchema,
    Experiment,
)


def make_gene_mapping_dict(entries_list: list) -> dict:
    """
    Create a dictionary mapping the provided entries to standardised gene symbols, e.g. "alpha-synuclein" -> "SNCA"

    Parameters:
    ----------
    entries_list: list, pd.Series
        List of MaveDB gene name entries to map

    Returns:
    -------
    dict:
        Dictionary mapping original gene names to standardized gene symbols
    """
    out_gene_mapping = {}
    # custom_mappings dict was compiled manually by examining all human entries in MAVEDB
    custom_mappings = {
        "alpha-synuclein": "SNCA",
        "Ras": "HRAS",
        "ErbB2": "ERBB2",
        "p53": "TP53",
        "COMT_ROI1_2": "COMT",
        "Aβ42": "APP",
        "hYAP65 WW domain": "YAP1",
        "BRCA1 translation start through RING domain": "BRCA1",
        "S505N MPL": "MPL",
        "human L-Selectin": "SELL",
        "PSD95 PDZ3": "DLG4",
        "Glycophorin A": "GYPA",
        "Minigene exon - Wilms’ tumor gene": "WT1",
        "Minigene exon - SRSF1 (ASF/SF2) binding site": "WT1",
        "Minigene exon - SRSF7 (9G8) binding site": "WT1",
        "Minigene exon - hnRNPA1 binding site": "WT1",
        "Minigene exon - hnRNP D binding site": "WT1",
        "Minigene exon - hnRNP I (PTB) binding site": "WT1",
        "Minigene exon - hnRNP L binding site": "WT1",
        "Minigene exon - CG-containing enhancer": "WT1",
        "Minigene exon - AC-rich enhancer": "WT1",
        "Minigene exon - pyrimidine sequence": "WT1",
        "NA Transcription factor IIS, N-terminal domain": "TCEA1",
        "NA Ubiquitin-like domain": "RPS27AP5",
    }

    for e in entries_list:
        # if in custom mappings use that
        if e in custom_mappings.keys():
            out_gene_mapping[e] = custom_mappings[e]
        # if the entry is a single word, uppercase it
        elif len(e.split(" ")) == 1:
            out_gene_mapping[e] = e.upper()
        # if the entry consists of multiple words, take the first word and uppercase it
        else:
            out_gene_mapping[e] = re.split(" ", e)[0].upper().strip()

    return out_gene_mapping


def edit_mavedb_metadata_df_columns(metadata_df: pd.DataFrame = None) -> pd.DataFrame:
    """
    Edit columns of the mavedb metadata dataframe to match the unified schema.

    Parameters:
    ----------
    metadata_df: pd.DataFrame
        flattened mavedb metadata dataframe generated from zenodo mavedb json

    Returns:
    -------
    pd.DataFrame
        metadata dataframe matching columns of the unified schema
    """
    # columns that exist and need to be mapped
    columns_mapping_dict = {
        "dataset_id": "index_urn3",
        "perturbation_name": "target_name",
        "perturbed_target_symbol": "gene_name",
        "perturbed_target_biotype": "gene_category",
        "species": "gene_targetSequence.taxonomy.organismName",
        "study_title": "publication_title",
        "study_uri": "publication_doi",
        "study_year": "publication_publicationYear",
        "first_author": "publication_first_author",
        "last_author": "publication_last_author",
        "experiment_title": "experiment_title",
        "experiment_summary": "experiment_shortDescription",
        "number_of_perturbed_targets": "scoreset_numVariants",
    }

    # get all columns in the unified schema
    columns_adding_list = list(ObsSchema.to_schema().columns.keys())
    # define a list of columns to drop - comes from data
    columns_to_drop = [
        "sample_id",
        "significant",
        "significance_criteria",
        "perturbed_target_coord",
        "perturbed_target_chromosome",
        "perturbed_target_chromosome_encoding",
        "perturbed_target_number",
        "perturbed_target_ensg",
        "guide_sequence",
    ]

    # create empty dataframe with all columns and drop unnecessary ones
    output_df = pd.DataFrame(columns=columns_adding_list).drop(columns=columns_to_drop)
    # map existing columns
    for out_col, in_col in columns_mapping_dict.items():
        output_df[out_col] = metadata_df[in_col]

    output_df["experiment_title"] = (
        output_df["experiment_title"] + "; " + metadata_df["scoreset_title"]
    )

    return output_df


def process_mavedb_metadata(
    mavedb_metadata_json_path: str = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process MaveDB metadata JSON into a flat DataFrame.

    Parameters:
    ----------
    mavedb_metadata_json_path : str
        Path to the MaveDB metadata JSON file.

    Returns:
    -------
    pd.DataFrame
        Flattened DataFrame containing MaveDB metadata.

    """
    with open(mavedb_metadata_json_path) as f:
        mdb = json.load(f)["experimentSets"]

    processed_df = pd.DataFrame()
    for d in mdb:
        urn1 = d["urn"]
        experiment_set = d["experiments"]
        experiment_set_df = pd.json_normalize(experiment_set)
        urn2 = experiment_set_df["urn"]

        if (
            experiment_set_df["primaryPublicationIdentifiers"][0] == []
            and experiment_set_df["secondaryPublicationIdentifiers"][0] == []
        ):
            publication_df = pd.DataFrame(
                columns=[
                    "publication_title",
                    "publication_abstract",
                    "publication_doi",
                    "publication_publicationYear",
                    "publication_url",
                    "publication_first_author",
                    "publication_last_author",
                    "index_urn1",
                ]
            )
        else:
            if experiment_set_df["primaryPublicationIdentifiers"][0] != []:
                publication_df = pd.json_normalize(
                    experiment_set_df["primaryPublicationIdentifiers"][0][0]
                )
            else:
                publication_df = pd.json_normalize(
                    experiment_set_df["secondaryPublicationIdentifiers"][0][0]
                )
            publication_df["first_author"] = " ".join(
                publication_df["authors"][0][0]["name"].split(", ")[::-1]
            )
            publication_df["last_author"] = " ".join(
                publication_df["authors"][0][-1]["name"].split(", ")[::-1]
            )
            publication_df = publication_df[
                [
                    "title",
                    "abstract",
                    "doi",
                    "publicationYear",
                    "url",
                    "first_author",
                    "last_author",
                ]
            ]
            publication_df.columns = [
                "publication_" + col for col in publication_df.columns
            ]
            publication_df["index_urn1"] = urn1

        experiment_set_df = experiment_set_df[
            ["title", "shortDescription", "methodText"]
        ]
        experiment_set_df.columns = [
            "experiment_" + col for col in experiment_set_df.columns
        ]
        experiment_set_df["index_urn1"] = urn1
        experiment_set_df["index_urn2"] = urn2
        experiment_set_df = pd.merge(
            experiment_set_df, publication_df, how="left", on="index_urn1"
        )

        score_set_df_full = pd.DataFrame()
        for i, exp in enumerate(experiment_set):
            score_set = exp["scoreSets"]
            score_set_df = pd.json_normalize(score_set)
            urn3 = score_set_df["urn"]
            score_set_df = score_set_df[["title", "methodText", "numVariants"]]
            score_set_df.columns = ["scoreset_" + col for col in score_set_df.columns]
            # score_set_df['index_urn1'] = urn1
            score_set_df["index_urn2"] = urn2[i]
            score_set_df["index_urn3"] = urn3

            target_genes_df_full = pd.DataFrame()
            for score in score_set:
                urn3_i = score["urn"]
                target_genes = score["targetGenes"]
                target_genes_df = pd.json_normalize(target_genes)
                if "externalIdentifiers" not in target_genes_df.columns:
                    target_genes_df["externalIdentifiers"] = [[]]
                if target_genes_df["externalIdentifiers"][0] != []:
                    identifiers_df = pd.json_normalize(
                        target_genes_df["externalIdentifiers"][0]
                    )
                    identifiers_df = identifiers_df[
                        ["identifier.dbName", "identifier.identifier"]
                    ]
                    identifiers_df = identifiers_df.set_index(
                        "identifier.dbName"
                    ).T.reset_index(drop=True)
                    identifiers_df.columns = [
                        "identifier_" + col for col in identifiers_df.columns
                    ]
                    identifiers_df["index_urn3"] = urn3_i
                else:
                    identifiers_df = pd.DataFrame(
                        columns=[
                            "identifier_Ensembl",
                            "identifier_RefSeq",
                            "identifier_UniProt",
                            "index_urn3",
                        ]
                    )

                for col in ["name", "category", "targetSequence.taxonomy.organismName"]:
                    if col not in target_genes_df.columns:
                        target_genes_df[col] = None
                target_genes_df = target_genes_df[
                    ["name", "category", "targetSequence.taxonomy.organismName"]
                ]
                target_genes_df.columns = [
                    "gene_" + col for col in target_genes_df.columns
                ]
                target_genes_df["index_urn3"] = urn3_i

                target_genes_df = pd.merge(
                    target_genes_df, identifiers_df, how="left", on="index_urn3"
                )
                target_genes_df_full = pd.concat(
                    [target_genes_df_full, target_genes_df], axis=0
                )

            score_set_df = pd.merge(
                score_set_df, target_genes_df_full, how="left", on="index_urn3"
            )
            score_set_df_full = pd.concat([score_set_df_full, score_set_df], axis=0)
        experiment_set_df = pd.merge(
            experiment_set_df, score_set_df_full, how="left", on="index_urn2"
        )
        processed_df = pd.concat([processed_df, experiment_set_df], axis=0)

    # filter for human only
    processed_df = processed_df[
        processed_df["gene_targetSequence.taxonomy.organismName"] == "Homo sapiens"
    ]
    # rename gene_name to target_name
    processed_df = processed_df.rename(columns={"gene_name": "target_name"})
    # map gene names
    mapping_dict = make_gene_mapping_dict(processed_df["target_name"])
    processed_df["gene_name"] = processed_df["target_name"].map(mapping_dict)
    # reorder columns
    processed_df = processed_df[
        ["index_urn1", "index_urn2", "index_urn3"]
        + [
            col
            for col in processed_df.columns
            if col not in ["index_urn1", "index_urn2", "index_urn3"]
        ]
    ]
    # order by index_urn3
    processed_df = processed_df.sort_values(by=["index_urn3"]).reset_index(drop=True)

    # edit the columns of the metadata based on the ObsSchema unified schema
    final_df = edit_mavedb_metadata_df_columns(processed_df)

    # add constant columns
    final_df["data_modality"] = "MAVE"

    return final_df, processed_df


def make_adata_mavedb(
    mavedb_id: str = None,
    mavedb_csv_dir: str = "../Dump/mavedb-dump.20250612164404/csv",
    curated_metadata_df: pd.DataFrame = None,
    save_h5ad_dir=None,
):
    """
    Create an AnnData object for a specific MaveDB experiment id.

    Parameters:
    ----------
    mavedb_id : str
        MaveDB experiment id (e.g. "urn:mavedb:00000001-b-1")
    mavedb_csv_dir : str
        Path to the directory containing MaveDB CSV files.
    curated_metadata_df : pd.DataFrame
        Curated metadata dataframe for MaveDB experiments.
    save_h5ad_dir : Path
        Directory to save the AnnData h5ad file. If None, the file is not saved.

    Returns:
    -------
    anndata.AnnData
        AnnData object containing the MaveDB data.
    """
    mavedb_data_path = f"{mavedb_csv_dir}/{mavedb_id.replace(':', '-')}.scores.csv"
    if not Path(mavedb_data_path).exists():
        raise FileNotFoundError(f"MaveDB data file not found at {mavedb_data_path}")
    mavedb_data = pd.read_csv(mavedb_data_path)

    X_df = mavedb_data.copy()
    X_df.index = X_df["accession"].str.replace("#", "-").str.replace(":", "-")
    X_df = X_df.drop(columns=["accession"])
    X_df = X_df.iloc[
        :, 3::
    ]  # first three columns are always hgvs ids, take all columns from 4th onwards as these are the scores

    metadata_subset_dict = curated_metadata_df[
        curated_metadata_df["dataset_id"] == mavedb_id
    ].to_dict(orient="records")[0]

    OBS_df = pd.DataFrame(index=X_df.index, data=metadata_subset_dict)

    OBS_df["sample_id"] = OBS_df.index

    # create perturbation_name column based on hgvs_nt and hgvs_pro
    hgvs_str_list = []
    for nt, pro in zip(mavedb_data["hgvs_nt"], mavedb_data["hgvs_pro"]):
        if isinstance(nt, float) and isinstance(pro, float):
            hgvs_str = None
        elif isinstance(nt, float) and isinstance(pro, str):
            hgvs_str = pro
        elif isinstance(nt, str) and isinstance(pro, float):
            hgvs_str = nt
        else:
            hgvs_str = f"{nt}|{pro}"
        hgvs_str_list.append(hgvs_str)
    OBS_df["perturbation_name"] = hgvs_str_list

    OBS_df["significant"] = None
    OBS_df["significance_criteria"] = None
    OBS_df["guide_sequence"] = None

    VAR_df = pd.DataFrame(index=X_df.columns, data={"score_name": X_df.columns})

    # replace None with np.nan to avoid issues with AnnData writing
    OBS_df = OBS_df.replace({None: np.nan})
    VAR_df = VAR_df.replace({None: np.nan})

    # Create AnnData object
    adata = anndata.AnnData(X=X_df, obs=OBS_df, var=VAR_df)

    # save the anndata object as an h5ad file
    h5ad_path = None
    if save_h5ad_dir:
        save_h5ad_dir.mkdir(parents=True, exist_ok=True)
        h5ad_path = save_h5ad_dir / f"{mavedb_id.replace(':', '-')}.h5ad"
        adata.write_h5ad(h5ad_path)
        print(f"✅ Saved AnnData object to {h5ad_path}")

    return adata, h5ad_path


def curate_mavedb(
    adata_h5ad_path: Path = None,
    save_curated_h5ad: bool = True,
    save_curated_parquet: bool = True,
    split_parquet: bool = True,
):
    """Curate DepMap AnnData object using curation tools.

    Parameters:
    ----------
        adata_h5ad_path: Path
            Path to the non-curated AnnData h5ad file.
        save_curated_h5ad: bool
            Whether to save the curated AnnData object as an h5ad file. Defaults to True.
        save_curated_parquet: bool
            Whether to save the curated data as a parquet file. Defaults to True.
        split_parquet: bool
            Whether to save separate Parquet files for data and metadata. Defaults to True.

    Returns:
    -------
        CuratedDataset
            CuratedDataset object containing the curated data.
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
    cur_data.validate_data(slot="obs", verbose=False)

    if save_curated_h5ad:
        cur_data.save_curated_data_h5ad()

    # save the curated data as Parquet files
    if save_curated_parquet:
        cur_data.save_curated_data_parquet(split_metadata=split_parquet)

    return cur_data


def process_mavedb(
    mavedb_dataset_id: str = None,
    mavedb_csv_dir: str = "../Dump/mavedb-dump.20250612164404/csv",
    curated_metadata_df: pd.DataFrame = None,
    non_curated_h5ad_dir: Path = Path("../non_curated/h5ad"),
    overwrite: bool = False,
):
    """
    Process and curate MaveDB dataset.
    
    Parameters:
    ----------
    mavedb_dataset_id : str
        MaveDB experiment id (e.g. "urn:mavedb:00000001-b-1")
    mavedb_csv_dir : str
        Path to the directory containing MaveDB CSV files.
    curated_metadata_df : pd.DataFrame
        Curated metadata dataframe for MaveDB experiments.
    non_curated_h5ad_dir : Path
        Directory to save the non-curated AnnData h5ad file.
    overwrite : bool
        Whether to overwrite existing curated data. Defaults to False.
        
    Returns:
    -------
    CuratedDataset
        CuratedDataset object containing the curated data.
    """

    # check if the data has been processed already, if yes, skip processing
    curated_h5ad_path = (
        Path(non_curated_h5ad_dir.as_posix().replace("non_curated", "curated"))
        / f"{mavedb_dataset_id}_curated.h5ad"
    )
    if curated_h5ad_path.exists():
        if overwrite:
            print(
                f"♻️ Curated DepMap data for {mavedb_dataset_id} already exists at {curated_h5ad_path}. Overwriting..."
            )
        else:
            print(
                f"✅ Curated DepMap data for {mavedb_dataset_id} already exists at {curated_h5ad_path}. Skipping processing."
            )
            return

    # Step 1: Create AnnData object for a specific MaveDB experiment id
    _adata, h5ad_path = make_adata_mavedb(
        mavedb_id=mavedb_dataset_id,
        mavedb_csv_dir=mavedb_csv_dir,
        curated_metadata_df=curated_metadata_df,
        save_h5ad_dir=non_curated_h5ad_dir,
    )

    # Step 2: Curate the AnnData object
    cur_data = curate_mavedb(
        adata_h5ad_path=h5ad_path,
        save_curated_h5ad=True,
        save_curated_parquet=True,
        split_parquet=True,
    )

    return cur_data
