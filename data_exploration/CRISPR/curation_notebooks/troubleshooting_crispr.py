# %%

import pandas as pd
import scanpy as sc

import os
import sys
from pathlib import Path

import requests
import xml.etree.ElementTree as ET

os.chdir("/Users/zakirov/Documents/GitHub/PerturbationCatalogue/data_exploration/")

# sys.path.append("../../")
# sys.path.append("/Users/zakirov/Documents/GitHub/PerturbationCatalogue/data_exploration/")

from curation_tools.curation_tools import CuratedDataset
from curation_tools.perturbseq_anndata_schema import ObsSchema, VarSchema

from curation_tools.unified_metadata_schema.unified_metadata_schema import Experiment


# %%

data_source_link = "downloads.thebiogrid.org/Download/BioGRID-ORCS/Latest-Release/BIOGRID-ORCS-ALL-homo_sapiens-LATEST.screens.tar.gz"

data_dir = Path("CRISPR/supplementary/temp/biogrid/data")
# Download and unpack the data
# !mkdir -p {data_dir}
# !wget {data_source_link} -O {data_dir}/biogrid_data.tar.gz
# !tar -xzvf {data_dir}/biogrid_data.tar.gz -C {data_dir}
# !rm {data_dir}/biogrid_data.tar.gz

biogrid_files = os.listdir(data_dir)
metadata_file = [e for e in biogrid_files if e.startswith("BIOGRID-ORCS-SCREEN_INDEX")][0]
data_files = [e for e in biogrid_files if e != metadata_file]
metadata_df = pd.read_csv(data_dir / metadata_file, sep="\t")

# %%


def extract_screen_id_from_file(file_name):
    """
    Extracts the screen id from the file name.
    The file name format is expected to be like 'BIOGRID-ORCS-SCREEN_12345-...'.
    """
    parts = file_name.split("_")
    if len(parts) > 1:
        return int(parts[1].split("-")[0])
    return None


def extract_screen_name_from_metadata(metadata_df, screen_id):
    """
    Extracts the screen name from the metadata DataFrame based on the screen id.
    """
    metadata_df_subset = metadata_df.loc[metadata_df["#SCREEN_ID"] == screen_id]
    screen_name = f"{metadata_df_subset['F_AUTHOR_L_NAME'].values[0]}_{metadata_df_subset['YEAR'].values[0]}_{metadata_df_subset['SCREEN_NAME'].values[0]}"
    screen_name = screen_name.replace('-', '_')
    return screen_name


def map_metadata_score_columns(screen_df, metadata_df, screen_id):

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


def doi_to_pmid(doi):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "pubmed", "term": doi, "retmode": "xml"}
    response = requests.get(base_url, params=params)
    response.raise_for_status()
    root = ET.fromstring(response.content)
    idlist = root.find("IdList")
    if idlist is not None and len(idlist) > 0:
        pmid_elem = idlist.find("Id")
        if pmid_elem is not None and pmid_elem.text:
            return pmid_elem.text.strip()
    return None


def get_pubmed_info(ids):
    """
    Fetch author, title, year, abstract, and link for PMIDs or DOIs.
    Returns: dict keyed by input ID, each value is a dict with keys:
      - 'title', 'year', 'abstract', 'link', 'authors' (list)
    """
    if not ids:
        return {}

    pmid_map = {}
    unresolved_ids = []
    for original_id in ids:
        if original_id.isdigit():
            pmid_map[original_id] = original_id
        else:
            pmid = doi_to_pmid(original_id)
            if pmid:
                pmid_map[original_id] = pmid
            else:
                unresolved_ids.append(original_id)

    if not pmid_map:
        return {id_: {} for id_ in ids}

    pmids_unique = list(set(pmid_map.values()))
    pmids_str = ",".join(pmids_unique)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmids_str}&retmode=xml"
    response = requests.get(url)
    response.raise_for_status()
    xml_root = ET.fromstring(response.content)

    pmid_metadata = {}
    for article in xml_root.findall("PubmedArticle"):
        medline_citation = article.find("MedlineCitation")
        if medline_citation is None:
            continue
        pmid_elem = medline_citation.find("PMID")
        if pmid_elem is None or pmid_elem.text is None:
            continue
        pmid = pmid_elem.text.strip()

        # Title
        article_title = medline_citation.findtext("Article/ArticleTitle")

        # Year: try standard Year, then MedlineDate as fallback
        year = medline_citation.findtext("Article/Journal/JournalIssue/PubDate/Year")
        if year is None:
            medlinedate = medline_citation.findtext(
                "Article/Journal/JournalIssue/PubDate/MedlineDate"
            )
            if medlinedate:
                year = medlinedate[:4]  # Take first 4 chars as year

        # Abstract
        abstract_elem = medline_citation.find("Article/Abstract/AbstractText")
        abstract = None
        if abstract_elem is not None:
            if abstract_elem.text:
                abstract = abstract_elem.text.strip()
            else:
                # If there are multiple AbstractText sections (sometimes labeled), concatenate them.
                abstract = " ".join(
                    part.text.strip()
                    for part in medline_citation.findall(
                        "Article/Abstract/AbstractText"
                    )
                    if part.text
                )

        # Authors
        authors = []
        author_list = medline_citation.find(".//AuthorList")
        if author_list is not None:
            for author in author_list.findall("Author"):
                lastname = author.findtext("LastName")
                forename = author.findtext("ForeName")
                authors.append({"LastName": lastname, "ForeName": forename})
        # DOI extraction
        doi = None
        for eloc in medline_citation.findall("Article/ELocationID"):
            if eloc.attrib.get("EIdType") == "doi" and eloc.text:
                doi = eloc.text.strip()
                break

        # DOI link
        doi_link = f"https://doi.org/{doi}" if doi else None

        pmid_metadata[pmid] = {
            "title": article_title,
            "year": year,
            "abstract": abstract,
            "doi": doi_link,
            "authors": authors,
        }

        article_title, year, abstract, doi_link, authors = (
            None,
            None,
            None,
            None,
            None,
        )

    # Prepare final output keyed by original id
    result = {}
    for original_id in ids:
        if original_id in pmid_map:
            pmid = pmid_map[original_id]
            result[original_id] = pmid_metadata.get(pmid, {})
        else:
            result[original_id] = {}

    return result


def curate_dataset(screen_path, metadata_df, non_curated_dir=None):

    if non_curated_dir is None:
        raise ValueError("Provide a directory to the non-curated h5ad folder")

    # read the screen data
    screen_df = pd.read_csv(screen_path, sep="\t")

    # Extract screen ID from the file name
    screen_id = extract_screen_id_from_file(screen_path.name)
    if screen_id is None:
        raise ValueError(
            f"Could not extract screen ID from file name: {screen_path.name}"
        )
    # extract screen name from metadata and screen id
    screen_name = extract_screen_name_from_metadata(metadata_df, screen_id)

    # get screen-specific metadata
    screen_metadata = metadata_df.loc[metadata_df["#SCREEN_ID"] == screen_id]
    screen_metadata_dict = screen_metadata.iloc[0].to_dict()
    screen_metadata_dict = {
        k: (None if v == "-" else v) for k, v in screen_metadata_dict.items()
    }

    # get the column mappings for the score of interest
    screen_df, new_columns = map_metadata_score_columns(
        screen_df, metadata_df, screen_id
    )
    
    # group_by 'IDENTIFIER_ID' and average the counts
    screen_df = screen_df.groupby([e for e in screen_df.columns if e not in ['Read counts']]).agg({'Read counts': 'mean'}).reset_index()

    index = screen_df["OFFICIAL_SYMBOL"]+ '_' + str(screen_id) +'_' +screen_df["IDENTIFIER_ID"]

    X_df = screen_df[new_columns].set_index(index, drop=True)

    OBS_df = pd.DataFrame(index=index)
    OBS_df["perturbation_name"] = OBS_df.index
    OBS_df["perturbed_target_symbol"] = screen_df["OFFICIAL_SYMBOL"].to_list()

    VAR_df = pd.DataFrame(index=new_columns)
    VAR_df["score_name"] = new_columns

    adata = sc.AnnData(X=X_df, obs=OBS_df, var=VAR_df)

    noncurated_fname = f"{non_curated_dir}/biogrid_{screen_name}.h5ad"

    adata.write(noncurated_fname)

    cur_data = CuratedDataset(
        obs_schema=ObsSchema,
        var_schema=VarSchema,
        exp_metadata_schema=Experiment,
        noncurated_path=noncurated_fname,
    )

    cur_data.load_data()

    # replace any present ORFs
    orf_map_dict = cur_data.gene_ont[["symbol", "synonyms"]]
    orf_map_dict["synonyms"] = orf_map_dict["synonyms"].str.split("|")
    orf_map_dict = orf_map_dict.explode("synonyms").drop_duplicates().dropna()
    orf_map_dict = orf_map_dict[
        orf_map_dict["synonyms"]
        .str.upper()
        .str.contains(r"^C\d{1,2}ORF\d.+", regex=True)
    ]
    orf_map_dict = orf_map_dict[
        orf_map_dict["synonyms"].isin(
            cur_data.adata.obs["perturbed_target_symbol"].str.upper()
        )
    ]
    orf_map_dict = orf_map_dict.set_index("synonyms")["symbol"].to_dict()
    
    cur_data.adata.obs['perturbed_target_symbol'] = cur_data.adata.obs['perturbed_target_symbol'].replace(orf_map_dict)

    cur_data.standardize_genes(
        slot="obs",
        input_column="perturbed_target_symbol",
        input_column_type="gene_symbol",
        multiple_entries=False,
    )

    cur_data.count_entries(
        slot="obs",
        input_column="perturbed_target_symbol",
        count_column_name="perturbed_target_number",
        sep="|",
    )

    cur_data.create_columns(
        slot="obs", col_dict={"treatment_label": screen_metadata_dict["CONDITION_NAME"]}
    )

    display(screen_df)


# %%

# replace the biorxiv link with doi for author mapping
metadata_df["SOURCE_ID"] = metadata_df["SOURCE_ID"].str.replace(
    "https://www.biorxiv.org/content/10.1101/2020.08.27.270819v1.full",
    "10.1101/2020.08.27.270819",
)
metadata_df["SOURCE_ID"] = metadata_df["SOURCE_ID"].str.replace(
    "10.1101/2023.01.23.525275", "10.1016/j.celrep.2023.112987"
)

# get pubmed info
pm_info = get_pubmed_info(list(metadata_df["SOURCE_ID"].unique()))

# manually replace 10.1101/2024.05.08.593110
pm_info["10.1101/2024.05.08.593110"] = {
    "title": "TMEM106B-mediated SARS-CoV-2 infection allows for robust ACE2-independent infection in vitro but not in vivo",
    "year": 2024,
    "abstract": "Angiotensin converting enzyme 2 (ACE2) serves as the primary entry receptor for severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2). However, ACE2-independent entry has been observed in vitro for SARS-CoV-2 strains containing the E484D amino acid substitution in the spike protein. In this study, we conducted a whole genome CRISPR-Cas9 knockout screen using a SARS-CoV-2 strain containing the spike-E484D substitution (SARS-CoV-2MA1) to identify the ACE2-independent entry mechanisms. Our findings revealed that SARS-CoV-2MA1 infection in HEK293T cells relied on heparan sulfate and endocytic pathways, with TMEM106B emerging as the most significant contributor. While SARS-CoV-2MA1 productively infected human brain organoids and K18-hACE2 mouse brains, it did not infect C57BL/6J or Ifnar-/- mouse brains. This suggests that ACE2-independent entry via TMEM106B, which is a protein that is predominantly expressed in the brain, did not overtly increase the risk of SARS-CoV-2 neuroinvasiveness in wild-type mice. Importantly, SARS-CoV-2MA1 did not replicate in Ace2-/- mouse respiratory tracts. Overall, this suggests that robust ACE2-independent infection by SARS-CoV-2E484D is likely a phenomenon specific to in vitro conditions, with no apparent clinical implications.",
    "doi": "https://doi.org/10.1101/2024.05.08.593110",
    "authors": [
        {"LastName": "Tau", "ForeName": "Kexin"},
        {"LastName": "Rawle", "ForeName": "Daniel J"},
    ],
}


metadata_df_curated = metadata_df[["#SCREEN_ID", "SCREEN_NAME", "SOURCE_ID"]]

metadata_df_curated["F_AUTHOR_F_NAME"] = metadata_df_curated["SOURCE_ID"].map(
    {k: v["authors"][0]["ForeName"] for k, v in pm_info.items()}
)
metadata_df_curated["F_AUTHOR_L_NAME"] = metadata_df_curated["SOURCE_ID"].map(
    {k: v["authors"][0]["LastName"] for k, v in pm_info.items()}
)
metadata_df_curated["L_AUTHOR_F_NAME"] = metadata_df_curated["SOURCE_ID"].map(
    {k: v["authors"][-1]["ForeName"] for k, v in pm_info.items()}
)
metadata_df_curated["L_AUTHOR_L_NAME"] = metadata_df_curated["SOURCE_ID"].map(
    {k: v["authors"][-1]["LastName"] for k, v in pm_info.items()}
)

metadata_df_curated["YEAR"] = metadata_df_curated["SOURCE_ID"].map(
    {k: v["year"] for k, v in pm_info.items()}
)

metadata_df_curated["TITLE"] = metadata_df_curated["SOURCE_ID"].map(
    {k: v["title"] for k, v in pm_info.items()}
)

metadata_df_curated["ABSTRACT"] = metadata_df_curated["SOURCE_ID"].map(
    {k: v["abstract"] for k, v in pm_info.items()}
)

metadata_df_curated["DOI"] = metadata_df_curated["SOURCE_ID"].map(
    {k: v["doi"] for k, v in pm_info.items()}
)

# screen format
metadata_df_curated["library_format_label"] = metadata_df["SCREEN_FORMAT"].map(
    {"Pool": "pooled", "Array": "arrayed", "in vivo": "in vivo"}
)

# software analysis
metadata_df_curated["software_analysis_label"] = metadata_df["ANALYSIS"]

# library perturbation type
metadata_df_curated["library_perturbation_type_label"] = metadata_df[
    "LIBRARY_METHODOLOGY"
].map({"Knockout": "knockout", "Inhibition": "inhibition", "Activation": "activation"})

# perturbation type
metadata_df_curated["perturbation_type_label"] = metadata_df["LIBRARY_TYPE"].map(
    {
        "CRISPRn": "CRISPRn",
        "CRISPRi": "CRISPRi",
        "CRISPRa": "CRISPRa",
        "Cytosine Base Editing-Mediated Gene KnockOut": "CBE",
    }
)

# library name and manufacturer
metadata_df["LIBRARY"] = metadata_df["LIBRARY"].str.replace(
    "(Gordon 2020)", "(Gordon, 2020)"
)
metadata_df_curated["library_name"] = metadata_df["LIBRARY"]


# %%
import re

l = metadata_df["LIBRARY"].unique()

lib_manufact_map_dict = {}
lib_name_map_dict = {}

for e in l:
    match_res = re.search(
        r"\(([^()]+?,\s*\d{4})\)", e, re.IGNORECASE
    )  # detect year - if detected, library manufacturer is present
    if match_res:
        match_str = match_res.group()
        match_str = match_str.lstrip("(").rstrip(")")
        manufact = match_str.split(",")[0]
        name = e.split(manufact)[0].rstrip("(").strip()
        lib_manufact_map_dict[e] = manufact
        lib_name_map_dict[e] = name
        # print(match_str)
    else:  # if no author provided then the manufacturer is not provided and only the name is available
        lib_name_map_dict[e] = e
        lib_manufact_map_dict[e] = None
        # print(e)

metadata_df_curated["library_manufacturer"] = metadata_df_curated["library_name"].map(
    lib_manufact_map_dict
)
metadata_df_curated["library_name"] = metadata_df_curated["library_name"].map(
    lib_name_map_dict
)

# %%

# library  generation method

metadata_df["ENZYME"].unique()

enzyme_map_dict = {
    "Cas9": "SpCas9",
    "dCas9-KRAB": "dCas9-KRAB",
    "sunCas9": "sunCas9",
    "SAM (NLS-dCas9-VP64/MS2-p65-HSF1)": "dCas9-SAM",
    "AsCpf1": "AsCpf1",
    "Cas9-v1": "SpCas9",
    "dCas9-VP64 + PP7-P65-HSF1": "dCas9-VP64",
    "dCas9-BFP-KRAB": "dCas9-KRAB",
    "dCas9-SunTag-P2A-HygR": "dCas9-Suntag",
    "dCAS-VP64_Blast (Zhu, 2021)": "dCas9-VP64",
    "AncBE4max": "AncBE4max",
    "dCas9-VP64": "dCas9-VP64",
    "dCas9–VPR": "dCas9–VPR",
    "Cas12a": "Cas12a",
    "dCas9-VP64 & p65-HSF1 (CRISPR SAM)": "dCas9-SAM",
    "dCas9": "dCas9",
    "iCAS9": "SpCas9",
}

metadata_df_curated["library_generation_method_label"] = metadata_df["ENZYME"].map(
    enzyme_map_dict
)

# Vertex autocuration
# %%
from time import sleep
from dotenv import load_dotenv
load_dotenv()
from curation_tools.vertex_curation_tools import (
    intitialize_vertexai,
    extract_text_with_elsapy,
    extract_data_with_gemini,
    fetch_pmc_json,
    extract_pmc_full_text
)
intitialize_vertexai()

# %%
# extract full text with elsapy for the articles that can be accessed with the API
full_text_dict = {}
for pmid, v in pm_info.items():
    doi_link = v.get("doi", None)
    if doi_link:
        doi = doi_link.replace("https://doi.org/", "")
        full_text = extract_text_with_elsapy(doi)
        if full_text:
            full_text_dict[pmid] = full_text
            sleep(2)  # to avoid rate limiting


# %%
# clean the text
import re

url_pattern = re.compile(r'https?://\S+|www\.\S+')
docs_pattern = re.compile(r'\S+\.(jpg|png|gif|tif|sml|zip|docx|xlsx)', re.IGNORECASE | re.MULTILINE)

clean_text_dict = {}
for k,v in full_text_dict.items():
    clean_text = v
    
    # remove everything after 'References'
    ref_pos = v.rfind('References')
    if ref_pos != -1:
        clean_text = clean_text[:ref_pos]
    
    # remove everything after 'Acknowledgements'
    ack_pos = v.rfind('Acknowledgements')
    if ack_pos != -1:
        clean_text = clean_text[:ack_pos]
        
    # remove everything before 'Introduction'
    intro_pos = v.find('Introduction')
    if intro_pos != -1:
        clean_text = clean_text[intro_pos:]

    # remove URLs
    clean_text = url_pattern.sub('', clean_text)
    # remove document links
    clean_text = docs_pattern.sub('', clean_text)
    
    clean_text_dict[k] = clean_text

# %%
# check the fraction of text remaining after cleaning
for k in full_text_dict.keys():
    print(f"PMID:{k}: {round(len(clean_text_dict[k])/len(full_text_dict[k]),2)}")

# %%
# extract data with gemini
from curation_tools.vertex_curation_tools import cleanup_response
import json

curated_data_dict = {}
os.makedirs(f"curated_gemini", exist_ok=True)

for pmid, v in clean_text_dict.items():
    try:
        curated_data = extract_data_with_gemini(v)
        curated_data_clean = cleanup_response(curated_data)
        curated_data_dict[pmid] = curated_data
        with open(f"curated_gemini/curated_{pmid}.json", "w") as f:
            json.dump(curated_data_clean, f, indent=4)
        sleep(2)  # to avoid rate limiting
    except Exception as e:
        print(f"Error processing PMID {pmid}: {e}")
        continue
# %%

# list of PMIDs curated with crossref API (i.e. available through crossref API)

curated_gemini_base_dir = '/Users/zakirov/Documents/GitHub/PerturbationCatalogue/data_exploration/curated_gemini'

cf_list = os.listdir(curated_gemini_base_dir)
cf_list_ids = [e.lstrip('curated_').rstrip('.json') for e in cf_list]

all_pmids = metadata_df_curated['SOURCE_ID'].unique()

remaining_pmids = [e for e in all_pmids if e not in cf_list_ids]

for pmid in remaining_pmids:
    fetch_pmc_json(pmid, save_dir='pmc_json')

downloaded_pmids_paths = [e for e in os.listdir('pmc_json') if e.endswith('.json')]

for pmid in downloaded_pmids_paths:
    extract_pmc_full_text(filepath=f"pmc_json/{pmid}", save_dir='pmc_full_text')


# %%
from docling.datamodel.base_models import InputFormat
from docling.document_converter import DocumentConverter, PdfFormatOption
from docling.datamodel.pipeline_options import (
    PdfPipelineOptions,
)

# Docling Parse without EasyOCR
# -------------------------

output_dir = Path("full_text_from_pdfs")

input_doc_file_mapping = {
    "1-s2.0-S2211124723009981-main.pdf": "10.1016_j.celrep.2023.112987",
    "nihpp-2020.08.27.270819.pdf": "10.1101_2020.08.27.270819",
    "nihpp-2021.01.19.427194v2.pdf": "10.1101_2021.01.19.427194",
    "2024.05.08.593110v1.full.pdf": "10.1101_2024.05.08.593110",
    "nihpp-rs1910932v1.pdf": "10.21203_rs.3.rs-1910932_v1",
    "nm.4219.pdf": "27869803",
    "pnas.201617467.pdf": "28611215",
    "jcs206425.pdf": "28775154",
    "pnas.201712176.pdf": "28894007",
    "e00302-17.pdf": "29038160",
    "849.pdf": "29440296",
    "s41556-018-0088-1.pdf": "29662178",
    "s41564-018-0159-x.pdf": "29736038",
    "30317623.pdf": "30317623",
    "s41586-019-0955-3.pdf": "30787439",
    "30844312.pdf": "30844312",
    "s41586-019-1103-9.pdf": "30971826",
    "MCB.00037-19.pdf": "31010806",
    "pnas.201900867.pdf": "31019072",
    "kfz037.pdf": "31059574",
    "JVI.00211-19.pdf": "31142663",
    "pnas.201906275.pdf": "31213543",
    "EMBR-20-e48235.pdf": "31353801",
    "s41388-019-0924-1.pdf": "31406246",
    "510.pdf": "31551363",
    "scitranslmed.aax2863.pdf": "31666400",
    "pnas.201909720.pdf": "31818950",
    "JVI.01751-19.pdf": "31915280",
    "kfaa071.pdf": "32421776",
    "pnas.202011645.pdf": "33122441",
    "pnas.202010723.pdf": "33172989",
    "33186476.pdf": "33186476",
    "pnas.202022120.pdf": "33483422",
    "s41588-021-00805-2.pdf": "33686287",
    "s41388-021-01692-x.pdf": "33742119",
    "s10549-021-06213-8.pdf": "33864166",
    "s41587-021-00944-1.pdf": "34155407",
    "s41588-021-00889-w.pdf": "34253920",
    "jvi.02042-21.pdf": "35420441",
    "jvi.00056-23.pdf": "37167561",
}

pipeline_options = PdfPipelineOptions()
pipeline_options.do_ocr = False
pipeline_options.do_table_structure = True
pipeline_options.table_structure_options.do_cell_matching = True

doc_converter = DocumentConverter(
    format_options={
        InputFormat.PDF: PdfFormatOption(pipeline_options=pipeline_options)
    }
)

for pdf_file in os.listdir('paper_pdfs'):
    if pdf_file.endswith('.pdf'):
        print(f"Processing {pdf_file}...")
        full_path = os.path.join('paper_pdfs', pdf_file)
        if pdf_file in input_doc_file_mapping:
            # convert the document
            conv_result = doc_converter.convert(full_path)
            
            ## Export results
            output_dir.mkdir(parents=True, exist_ok=True)
            doc_filename = input_doc_file_mapping[pdf_file]
            with (output_dir / f"{doc_filename}.txt").open("w", encoding="utf-8") as f:
                f.write(conv_result.document.export_to_markdown())
        else:
            print(f"PDF file {pdf_file} not in mapping dictionary.")



# %%

# manually verifying studies

x = metadata_df_curated[['SOURCE_ID', 'SCREEN_NAME']].groupby('SOURCE_ID')['SCREEN_NAME'].count()
x = x.rename('study_count')

one_study = x[x == 1].index.tolist()

cf_list_fnames_one_study = [f"curated_{e}.json" for e in one_study if e in cf_list_ids]

with open(os.path.join(curated_gemini_base_dir, cf_list_fnames_one_study[0]), "r") as f:
    mdata = json.load(f)


# df from the dict

mdata_df = pd.json_normalize(mdata, sep='_').T

data_df = pd.read_csv(os.path.join(data_dir,  'BIOGRID-ORCS-SCREEN_1174-1.1.17.screen.tab.txt'), sep='\t')


# %%
curate_dataset(data_dir / 'BIOGRID-ORCS-SCREEN_1174-1.1.17.screen.tab.txt', metadata_df_curated, non_curated_dir='CRISPR/non_curated/h5ad')












