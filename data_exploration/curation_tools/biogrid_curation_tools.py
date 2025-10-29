from pathlib import Path

from anndata import AnnData
import numpy as np
import pandas as pd

from data_exploration.curation_tools.curation_tools import (
    upload_parquet_to_bq,
    CuratedDataset,
    ObsSchema,
    VarSchema,
    Experiment
)


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
    
    Returns:
    - AnnData: The created AnnData object.
    - str: Path to the saved h5ad file.
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
    h5ad_path = None
    if save_h5ad_dir:
        h5ad_path = Path(save_h5ad_dir) / f"{biogrid_dataset_id}.h5ad"
        adata.write_h5ad(h5ad_path)
        print(f"✅ Saved AnnData object to {h5ad_path}")

    return adata, h5ad_path


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
    biogrid_screen_id : str
        Specific biogrid screen ID to compare metadata for
    gemini_id : str
        Specific gemini ID to compare metadata for

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
    save_dir : str
        Directory to save the h5ad file if save_h5ad is True

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

    adata = AnnData(X=X_df, obs=OBS_df, var=VAR_df)

    if save_h5ad:
        adata_h5ad_save_path = save_dir / f"{biogrid_screen_id}.h5ad"
        adata.write(adata_h5ad_save_path)
        print(f"✅Saved adata to {adata_h5ad_save_path}")
        return adata, adata_h5ad_save_path
    else:
        return adata

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