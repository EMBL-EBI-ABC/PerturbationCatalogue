import numpy as np
import pandas as pd
import json
import re

from pathlib import Path
from anndata import AnnData

from CRISPR.supplementary.depmap import depmap_mapping as dm
from curation_tools.curation_tools import (
    CuratedDataset,
    ObsSchema,
    VarSchema,
    Experiment,
)


def process_depmap_metadata_df(
    df: pd.DataFrame,
    age_mapping: dict = None,
    sex_mapping: dict = None,
    tissue_mapping: dict = None,
    disease_mapping: dict = None,
) -> pd.DataFrame:
    """Process DepMap metadata dataframe to match Perturbation Catalogue unified schema.

    Parameters:
    -----------
        df (pd.DataFrame): DepMap metadata dataframe.
        age_mapping (dict): Mapping dictionary for age categories. Defaults to None.
        sex_mapping (dict): Mapping dictionary for sex categories. Defaults to None.
        tissue_mapping (dict): Mapping dictionary for tissue categories. Defaults to None.
        disease_mapping (dict): Mapping dictionary for disease categories. Defaults to None.

    Returns:
    -------
        pd.DataFrame: Processed metadata dataframe.
    """
    # Define a mapping dictionary
    column_mapping = {
        "ModelID": "ModelID",
        "Library": "Library",
        "CellLineName": "cell_line_label",
        "OncotreeSubtype": "disease_label",
        "AgeCategory": "developmental_stage_label",
        "Sex": "sex_label",
        "ModelType": "model_system_label",
        "Days": "timepoint",
    }

    model_system_mapping = {"Cell Line": "cell_line"}

    mapping_keys = list(column_mapping.keys())
    mapping_values = list(column_mapping.values())

    # Subset the DataFrame to include only the mapped columns
    df_mapped = df[mapping_keys].copy()
    df_mapped.columns = mapping_values

    # Replace values in specific columns to match unified schema requirements
    df_mapped["developmental_stage_label"] = df_mapped["developmental_stage_label"].map(
        age_mapping
    )
    df_mapped["sex_label"] = df_mapped["sex_label"].map(sex_mapping)
    df_mapped["tissue_label"] = df_mapped["disease_label"].map(
        tissue_mapping
    )  # mapping from disease is correct!
    df_mapped[["disease_label", "disease_id"]] = (
        df_mapped["disease_label"].map(disease_mapping).to_list()
    )
    df_mapped["model_system_label"] = df_mapped["model_system_label"].map(
        model_system_mapping
    )

    df_mapped["timepoint"] = df_mapped["timepoint"].map(
        {0: "P0D", 14: "P14D", 15: "P15D", 21: "P21D", 22: "P22D"}
    )

    df_days_library_for_merge = pd.DataFrame(
        data={
            "Library": ["Avana", "KY", "Humagne-CD"],
            "library_name": [
                "Avana",
                "Human Improved Genome-wide Knockout CRISPR Library v1 (KY)",
                "Humagne Set C and Set D Human CRISPR Knockout Libraries",
            ],
            "library_uri": [
                "https://www.addgene.org/pooled-library/yusa-crispr-knockout-human-v1/",
                "https://www.addgene.org/pooled-library/yusa-crispr-knockout-human-v1/",
                "https://www.addgene.org/pooled-library/broadgpp-human-knockout-humagne/",
            ],
            "library_generation_method_label": ["SpCas9", "SpCas9", "AsCas12a"],
            "library_generation_method_id": [
                "EFO:0022876",
                "EFO:0022876",
                "EFO:0022878",
            ],
            "library_manufacturer": [
                "Broad Institute",
                "Kosuke Yusa lab",
                "Doench and Root labs",
            ],
            "library_grnas_per_target": ["6", "5", "2"],
            "library_total_grnas": ["110257", "90709", "40710"],
            "library_lentiviral_generation": ["3", "3", "2"],
        }
    )

    df_mapped = df_mapped.merge(
        df_days_library_for_merge, on="Library", how="left"
    ).drop(columns=["Library"])

    return df_mapped


def get_cellosaurus_ontology(
    output_json_path: str,
    overwrite: bool = False,
    cellosaurus_ontology_csv_path: str = "../supplementary/depmap/cellosaurus_derived_ontology.csv",
) -> pd.DataFrame:
    """
    Download and convert Cellosaurus OBO ontology to Obograph JSON format.

    Parameters:
    ----------
        output_json_path (str): Path to save the converted Cellosaurus ontology in Obograph JSON format.
        overwrite (bool): Whether to overwrite the existing JSON file if it exists. Defaults to False.
        cellosaurus_ontology_csv_path (str): Path to save the processed Cellosaurus ontology dataframe as a CSV file. Defaults to "../supplementary/depmap/cellosaurus_derived_ontology.csv".
    Returns:
    -------
        pd.DataFrame: Processed Cellosaurus ontology dataframe.
    """
    from bioontologies import convert_to_obograph
    import json

    # if the cell line ontology mapping already exists, load it and return
    if Path(cellosaurus_ontology_csv_path).exists():
        print(
            f"✅ Cellosaurus ontology mapping already exists at {cellosaurus_ontology_csv_path}. Loading..."
        )
        codf = pd.read_csv(cellosaurus_ontology_csv_path)
        return codf

    # download and convert cellosaurus ontology to obograph json format if it doesn't exist
    if not Path(output_json_path).exists():
        print(
            f"⬇️ Downloading and converting Cellosaurus ontology to {output_json_path}..."
        )
        convert_to_obograph(
            input_path="https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo",
            input_flag="-I",
            json_path=output_json_path,
        )
    else:
        if overwrite:
            print(
                f"♻️ Cellosaurus ontology already exists at {output_json_path}. Overwriting..."
            )
            print(
                f"⬇️ Downloading and converting Cellosaurus ontology to {output_json_path}..."
            )
            convert_to_obograph(
                input_path="https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo",
                input_flag="-I",
                json_path=output_json_path,
            )
        else:
            print(
                f"✅ Cellosaurus ontology already exists at {output_json_path}. Skipping download."
            )

    # read json with cellosaurus ontology
    with open(output_json_path) as f:
        co = json.load(f)["graphs"][0]["nodes"]
        print(
            f"✅ Loaded Cellosaurus ontology from {output_json_path} with {len(co)} entries."
        )

    # convert to dataframe
    codf = pd.DataFrame(co)
    # filter out irrelevant entries
    codf = codf[codf["id"].str.contains("Cellosaurus")]
    # drop rows with missing labels
    codf = codf.dropna(subset=["lbl"])
    # drop irrelevant columns
    codf = codf.drop(columns=["type", "propertyType"])

    # extract relevant information from existing columns
    print("ℹ️ Extracting relevant information from Cellosaurus ontology...")
    codf["comments"] = codf["meta"].apply(
        lambda x: x.get("comments", []) if isinstance(x, dict) else []
    )
    codf["subsets"] = codf["meta"].apply(
        lambda x: x.get("subsets", []) if isinstance(x, dict) else []
    )
    codf["xrefs"] = codf["meta"].apply(
        lambda x: x.get("xrefs", []) if isinstance(x, dict) else []
    )
    codf["synonyms"] = codf["meta"].apply(
        lambda x: x.get("synonyms", []) if isinstance(x, dict) else []
    )
    codf["subsets"] = codf["subsets"].apply(
        lambda x: [s.split("Cellosaurus#")[1] for s in x] if isinstance(x, list) else []
    )

    # extract synonyms
    syn_dict = {}
    for e, lbl in zip(codf["synonyms"], codf["lbl"]):
        syn_list_temp = []
        for i in e:
            for k, v in i.items():
                if k == "pred" and v == "hasRelatedSynonym":
                    syn_list_temp.append(i.get("val", []))

        # add the label itself to the list of synonyms
        syn_dict[lbl] = [lbl] + syn_list_temp

    # extract cell type (CL) information
    celltype_dict = {}
    for e, lbl in zip(codf["comments"], codf["lbl"]):
        if e != []:
            e = e[0]
            if "CL=CL" in e:
                celltype_dict[lbl] = re.search(r"(?<=CL=)CL_\d+", e).group(0)

    # extract cell line (CLO) information
    celline_dict = {}
    for lbl, xref_list in zip(codf["lbl"], codf["xrefs"]):
        if xref_list != []:
            for xref in xref_list:
                if "CLO:CLO" in xref["val"]:
                    clo_entry = xref["val"].split(":")[1]
                    if lbl in celline_dict.keys():
                        celline_dict[lbl] = celline_dict[lbl] + [clo_entry]
                    else:
                        celline_dict[lbl] = [clo_entry]

    # map extracted information to the dataframe
    print("ℹ️ Mapping extracted information to the dataframe...")
    codf["synonyms"] = codf["lbl"].map(syn_dict)
    codf["cell_type_id"] = codf["lbl"].map(celltype_dict).str.replace("_", ":")
    codf["cell_line_id"] = codf["lbl"].map(celline_dict)
    codf["cell_line_id"] = [
        e[0] if isinstance(e, list) else None for e in codf["cell_line_id"]
    ]
    codf["cell_line_id"] = codf["cell_line_id"].str.replace("_", ":")

    # explode on synonyms to have one synonym per row for simpler mapping
    codf = codf.explode("synonyms")

    # rename lbl to cell_line_label
    codf = codf.rename(columns={"lbl": "cell_line_label"})

    # keep only relevant columns
    codf = codf[["cell_line_label", "synonyms", "cell_type_id", "cell_line_id"]]

    # save the cellosaurus ontology dataframe as a csv file
    codf.to_csv(cellosaurus_ontology_csv_path, index=False)

    return codf


def make_adata_depmap(
    depmap_model_id: str,
    depmap_data: pd.DataFrame,
    depmap_metadata: pd.DataFrame,
    save_h5ad_dir: Path = None,
    cellosaurus_ont_json_path: str = "../supplementary/depmap/cellosaurus_ontology.json",
) -> tuple[AnnData, str | None]:
    """Create an AnnData object from DepMap data and metadata.

    Parameters:
    ----------
        depmap_model_id (str): DepMap model ID (e.g. ACH-000001).
        depmap_data (pd.DataFrame): DepMap data dataframe.
        depmap_metadata (pd.DataFrame): DepMap metadata dataframe.
        save_h5ad_dir (Path, optional): Directory to save the AnnData object as an h5ad file. Defaults to None.
        cellosaurus_ont_json_path (str, optional): Path to the Cellosaurus ontology JSON file. Defaults to "../supplementary/depmap/cellosaurus_ontology.json".
    Returns:
    -------
        AnnData: AnnData object for a given DepMap cell line identified by its model ID.
        str: Path to the saved h5ad file (if applicable).
    """
    dataset_id = f"depmap_{depmap_model_id.replace('-', '')}"

    score_name = ["Gene dependency probability estimate"]
    # Subset data for the given model ID
    X_df = depmap_data[[depmap_model_id]].copy()
    cleaned_index = [
        gene.split(" (")[0].replace(".", "").replace("-", "") for gene in X_df.index
    ]
    X_df.index = [
        f"depmap_{depmap_model_id.replace('-', '')}_{gene}" for gene in cleaned_index
    ]
    X_df.columns = score_name

    # Extraxt metadata from depmap metadata and manual mappings
    metadata_subset = process_depmap_metadata_df(
        df=depmap_metadata,
        age_mapping=dm.age_mapping,
        sex_mapping=dm.sex_mapping,
        tissue_mapping=dm.tissue_mapping,
        disease_mapping=dm.disease_mapping,
    )
    metadata_subset_dict = (
        metadata_subset[metadata_subset["ModelID"] == depmap_model_id]
        .copy()
        .to_dict(orient="records")[0]
    )
    OBS_df = pd.DataFrame(index=X_df.index, data=metadata_subset_dict)
    OBS_df = OBS_df.drop(columns=["ModelID"])

    OBS_df["dataset_id"] = dataset_id
    OBS_df["sample_id"] = OBS_df.index
    OBS_df["perturbation_name"] = X_df.index
    OBS_df["perturbed_target_symbol"] = cleaned_index
    # significance criteria taken from Dependent Cell Lines box of gene overview page, e.g. https://depmap.org/portal/gene/BRCA1?tab=overview
    OBS_df["significant"] = X_df[score_name] > 0.5

    metadata_dict = {
        "significance_criteria": "Gene dependency probability estimate > 0.5",
        "data_modality": "CRISPR screen",
        "guide_sequence": None,
        "perturbation_type_label": "CRISPRn",
        "perturbation_type_id": None,
        "treatment_label": None,
        "treatment_id": None,
        "species": "Homo sapiens",
        "cell_type_label": None,
        "cell_type_id": None,
        "sex_id": None,
        "developmental_stage_id": None,
        "disease_id": None,
        "model_system_id": None,
        "study_title": "DepMap, Broad (2025). DepMap Public 25Q3. Dataset.",
        "study_year": 2025,
        "study_uri": "https://depmap.org",
        "first_author": "Rand Arafeh",
        "last_author": "Francisca Vazquez",
        "experiment_title": f"DepMap whole-genome CRISPR-knockout screen in {metadata_subset_dict['cell_line_label']} cancer cell line.",
        "experiment_summary": f"""
            Whole-genome CRISPR knockout screen in {metadata_subset_dict['cell_line_label']} cancer cell line was carried out
            at Wellcome Sanger Institute (Project Score) and Broad Institute (Project Achilles) to identify essential survival genes.
            The two centres used different CRISPR libraries, with Wellcome Sanger Institute using Human Improved Genome-wide Knockout CRISPR Library v1
            and Broad Institute using GeCKOv2 and Avana libraries. Cells were cultured for 2-3 weeks post-transduction (depending on the centre).
        """,
        "number_of_perturbed_targets": X_df.shape[0],
        "number_of_perturbed_samples": "33000000-36000000",
        "library_generation_type_label": "Endogenous genetic perturbation method",
        "library_generation_type_id": "EFO:0022868",
        "enzyme_delivery_method_label": "lentivirus transduction",
        "enzyme_delivery_method_id": None,
        "library_delivery_method_label": "lentivirus transduction",
        "library_delivery_method_id": None,
        "enzyme_integration_state_label": "random locus integration",
        "enzyme_integration_state_id": None,
        "library_integration_state_label": "random locus integration",
        "library_integration_state_id": None,
        "enzyme_expression_control_label": "constitutive transgene expression",
        "enzyme_expression_control_id": None,
        "library_expression_control_label": "constitutive transgene expression",
        "library_expression_control_id": None,
        "library_format_label": "pooled",
        "library_format_id": None,
        "library_scope_label": "genome-wide",
        "library_scope_id": None,
        "library_perturbation_type_label": "knockout",
        "library_perturbation_type_id": None,
        "library_total_variants": None,
        "readout_dimensionality_label": "single-dimensional assay",
        "readout_dimensionality_id": None,
        "readout_type_label": "phenotypic",
        "readout_type_id": None,
        "readout_technology_label": "population growth assay",
        "readout_technology_id": None,
        "method_name_label": "proliferation CRISPR screen",
        "method_name_id": None,
        "method_uri": None,
        "sequencing_library_kit_label": None,
        "sequencing_library_kit_id": None,
        "sequencing_platform_label": "Illumina HiSeq 2000",
        "sequencing_platform_id": None,
        "sequencing_strategy_label": "direct sequencing",
        "sequencing_strategy_id": None,
        "software_counts_label": "custom",
        "software_counts_id": None,
        "software_analysis_label": "Achilles",
        "software_analysis_id": None,
        "reference_genome_label": "GRCh38",
        "reference_genome_id": None,
        "associated_datasets": json.dumps(
            [
                {
                    "dataset_accession": None,
                    "dataset_uri": "https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q3&filename=CRISPRGeneDependency.csv",
                    "dataset_description": "Post-Chronos gene dependency probability estimates for all models in the integrated gene effect.",
                    "dataset_file_name": "CRISPRGeneDependency.csv",
                }
            ]
        ),
    }

    for key, value in metadata_dict.items():
        if key not in OBS_df.columns:
            OBS_df[key] = value

    # download and process cellosaurus ontology
    cellosaurus_ontology = get_cellosaurus_ontology(
        cellosaurus_ont_json_path, overwrite=False
    )

    # map cellosaurus ontology information to the metadata
    print("ℹ️ Mapping Cellosaurus ontology information to the metadata...")
    co_fields = ["cell_line_label", "cell_line_id", "cell_type_id"]
    for field in co_fields:
        mapping_dict = dict(
            zip(cellosaurus_ontology["synonyms"], cellosaurus_ontology[field])
        )
        if field not in OBS_df.columns:
            OBS_df[field] = np.nan  # create the column if it doesn't exist
        OBS_df[field] = (
            OBS_df["cell_line_label"].map(mapping_dict).fillna(OBS_df[field])
        )  # fillna - if there are no matches, keep the existing value

    # VAR
    VAR_df = pd.DataFrame(index=score_name, data={"score_name": score_name})

    # replace None with np.nan to avoid issues with AnnData writing
    OBS_df = OBS_df.replace({None: np.nan})
    VAR_df = VAR_df.replace({None: np.nan})

    # Create AnnData object
    adata = AnnData(X=X_df, obs=OBS_df, var=VAR_df)

    # save the anndata object as an h5ad file
    h5ad_path = None
    if save_h5ad_dir:
        save_h5ad_dir.mkdir(parents=True, exist_ok=True)
        h5ad_path = save_h5ad_dir / f"depmap_{depmap_model_id}.h5ad"
        adata.write_h5ad(h5ad_path)
        print(f"✅ Saved AnnData object to {h5ad_path}")

    return adata, h5ad_path


def curate_depmap(
    adata_h5ad_path: Path = None,
    save_curated_h5ad: bool = True,
    save_curated_parquet: bool = True,
    split_parquet: bool = True,
):
    """Curate DepMap AnnData object using curation tools.

    Parameters:
    ----------
        adata_h5ad_path: Path to the non-curated AnnData h5ad file.
        save_curated_h5ad: Whether to save the curated AnnData object as an h5ad file. Defaults to True.
        save_curated_parquet: Whether to save the curated data as a parquet file. Defaults to True.
        split_parquet: Whether to save separate Parquet files for data and metadata. Defaults to True.

    Returns:
    -------
        CuratedDataset: CuratedDataset object containing the curated data.
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

    # add tissue information
    cur_data.standardize_ontology(
        input_column="tissue_label",
        column_type="term_name",
        ontology_type="tissue",
        overwrite=True,
    )

    cur_data.standardize_ontology(
        input_column="cell_type_id",
        column_type="term_id",
        ontology_type="cell_type",
        overwrite=True,
    )

    # if tissue_label is 'soft tissue', add an empty tissue_id
    if "soft tissue" in cur_data.adata.obs["tissue_label"].values:
        cur_data.adata.obs["tissue_id"] = None

    # if all tissue labels are empty, set tissue_id to None
    if cur_data.adata.obs["tissue_label"].isna().all():
        cur_data.adata.obs["tissue_id"] = None

    # match the order of columns in obs to the schema
    cur_data.match_schema_columns(slot="obs")

    # validate the data against the schema
    cur_data.validate_data(slot="obs", verbose=False)

    # save the curated AnnData object
    if save_curated_h5ad:
        cur_data.save_curated_data_h5ad()

    # save the curated data as Parquet files
    if save_curated_parquet:
        cur_data.save_curated_data_parquet(split_metadata=split_parquet)

    return cur_data


def process_depmap(
    depmap_dataset_id: str = None,
    depmap_data: pd.DataFrame = None,
    depmap_metadata: pd.DataFrame = None,
    non_curated_h5ad_dir: Path = Path("../non_curated/h5ad/depmap"),
    overwrite: bool = False,
):
    """Process and curate DepMap data for a specific cell line/dataset id.
    Parameters:
    ----------
        depmap_dataset_id (str): DepMap dataset ID (e.g. ACH-000001).
        depmap_data (pd.DataFrame): DepMap data dataframe.
        depmap_metadata (pd.DataFrame): DepMap metadata dataframe.
        non_curated_h5ad_dir (Path): Directory to save non-curated h5ad files.
        overwrite (bool): Whether to overwrite existing curated data. Defaults to False.
    Returns:
    -------
        CuratedDataset: CuratedDataset object containing the curated data.
    """
    # check if the data has been processed already, if yes, skip processing
    curated_h5ad_path = (
        Path(non_curated_h5ad_dir.as_posix().replace("non_curated", "curated"))
        / f"depmap_{depmap_dataset_id}_curated.h5ad"
    )
    if curated_h5ad_path.exists():
        if overwrite:
            print(
                f"♻️ Curated DepMap data for {depmap_dataset_id} already exists at {curated_h5ad_path}. Overwriting..."
            )
        else:
            print(
                f"✅ Curated DepMap data for {depmap_dataset_id} already exists at {curated_h5ad_path}. Skipping processing."
            )
            return

    # Step 1: Create AnnData object for a specific DepMap cell line
    _adata, h5ad_path = make_adata_depmap(
        depmap_model_id=depmap_dataset_id,
        depmap_data=depmap_data,
        depmap_metadata=depmap_metadata,
        save_h5ad_dir=non_curated_h5ad_dir,
    )

    # Step 2: Curate the AnnData object
    cur_data = curate_depmap(
        adata_h5ad_path=h5ad_path,
        save_curated_h5ad=True,
        save_curated_parquet=True,
        split_parquet=True,
    )

    return cur_data
