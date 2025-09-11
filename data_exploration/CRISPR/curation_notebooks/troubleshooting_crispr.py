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
metadata_file = [e for e in biogrid_files if e.startswith('BIOGRID-ORCS-SCREEN_INDEX')][0]
data_files = [e for e in biogrid_files if e != metadata_file]
metadata_df = pd.read_csv(data_dir / metadata_file, sep="\t")

# %%

def extract_screen_id_from_file(file_name):
    """
    Extracts the screen id from the file name.
    The file name format is expected to be like 'BIOGRID-ORCS-SCREEN_12345-...'.
    """
    parts = file_name.split('_')
    if len(parts) > 1:
        return int(parts[1].split('-')[0])
    return None

def extract_screen_name_from_metadata(metadata_df, screen_id):
    """
    Extracts the screen name from the metadata DataFrame based on the screen id.
    """
    screen_author = metadata_df.loc[metadata_df['#SCREEN_ID'] == screen_id, 'AUTHOR'].values[0]
    screen_author = screen_author.split(' ')
    screen_name = f"{screen_author[0]}_{screen_author[2].replace('(', '').replace(')','')}"
    return screen_name

def map_metadata_score_columns(screen_df, metadata_df, screen_id):

    # score columns in metadata are named SCORE.1_TYPE, SCORE.2_TYPE, etc.
    score_cols_meta = [f"SCORE.{e}_TYPE" for e in range(1,6)]
    
    # create a DataFrame with screen id and score columns
    meta_subset_df = metadata_df[['#SCREEN_ID']+score_cols_meta]
    
    # melt the DataFrame to have a long format with screen id, score type, and score value
    melt_meta_subset_df = meta_subset_df.melt(id_vars='#SCREEN_ID', value_vars=score_cols_meta, var_name='source_col', value_name='score_type')
    
    # remove the '_TYPE' suffix from the source_col because SCORE columns in screen data do not have this suffix
    melt_meta_subset_df['source_col'] = melt_meta_subset_df['source_col'].str.rstrip('_TYPE')
    
    # filter the melted DataFrame for the given screen id and create a mapping dictionary
    # where the keys are the source_col (SCORE.1, SCORE.2, etc.) and the values are the score_type (e.g. 'MaGeCK Score')
    map_dict = melt_meta_subset_df.loc[melt_meta_subset_df['#SCREEN_ID']==screen_id].set_index('source_col')['score_type'].to_dict()
    # new columns to be returned by the function
    new_columns = map_dict.values()
    
    # rename the columns in the screen DataFrame using the mapping dictionary
    screen_df = screen_df.rename(columns=map_dict)
    
    # drop the '-' empty column
    if '-' in screen_df.columns:
        screen_df = screen_df.drop(columns=['-'])
        new_columns = [e for e in new_columns if e != '-']

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
            medlinedate = medline_citation.findtext("Article/Journal/JournalIssue/PubDate/MedlineDate")
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
                    part.text.strip() for part in medline_citation.findall("Article/Abstract/AbstractText") if part.text
                )

        # Authors
        authors = []
        author_list = medline_citation.find(".//AuthorList")
        if author_list is not None:
            for author in author_list.findall("Author"):
                lastname = author.findtext("LastName")
                forename = author.findtext("ForeName")
                authors.append({
                    "LastName": lastname,
                    "ForeName": forename
                })
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
            "authors": authors
        }
        
        article_title, year, abstract, doi_link, authors = (None, None, None, None, None)

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
        raise ValueError('Provide a directory to the non-curated h5ad folder')

    # read the screen data
    screen_df = pd.read_csv(screen_path, sep="\t")

    # Extract screen ID from the file name
    screen_id = extract_screen_id_from_file(screen_path.name)
    if screen_id is None:
        raise ValueError(f"Could not extract screen ID from file name: {screen_path.name}")
    # extract screen name from metadata and screen id
    screen_name = extract_screen_name_from_metadata(metadata_df, screen_id)
    
    # get screen-specific metadata
    screen_metadata = metadata_df.loc[metadata_df['#SCREEN_ID']==screen_id]
    screen_metadata_dict = screen_metadata.to_dict(orient='index')[screen_metadata.index[0]]
    screen_metadata_dict = {k: (None if v == '-' else v) for k, v in screen_metadata_dict.items()}
    
    # get the column mappings for the score of interest
    screen_df, new_columns = map_metadata_score_columns(screen_df, metadata_df, screen_id)
    
    X_df = screen_df[['OFFICIAL_SYMBOL'] + new_columns].set_index('OFFICIAL_SYMBOL', drop=True)
    
    OBS_df = pd.DataFrame(index=screen_df['OFFICIAL_SYMBOL'])
    OBS_df['perturbation_name'] = list(screen_df['#SCREEN_ID'].astype(str) + '_' + OBS_df.index)
    OBS_df['perturbed_target_symbol'] = OBS_df.index
    
    VAR_df = pd.DataFrame(index=new_columns)
    VAR_df['score_name'] = new_columns

    adata = sc.AnnData(X=X_df, obs=OBS_df, var=VAR_df)
    
    noncurated_fname = f"{non_curated_dir}/biogrid_{screen_name}.h5ad"
    
    adata.write(noncurated_fname)
    
    cur_data = CuratedDataset(
        obs_schema=ObsSchema,
        var_schema=VarSchema,
        exp_metadata_schema=Experiment,
        noncurated_path = noncurated_fname
    )
    
    cur_data.load_data()
    
    # replace any present ORFs
    orf_map_dict = cur_data.gene_ont[['symbol', 'synonyms']]
    orf_map_dict['synonyms'] = orf_map_dict['synonyms'].str.split('|')
    orf_map_dict = orf_map_dict.explode('synonyms').drop_duplicates().dropna()
    orf_map_dict = orf_map_dict[orf_map_dict['synonyms'].str.upper().str.contains(r"^C\d{1,2}ORF\d.+", regex=True)]
    orf_map_dict = orf_map_dict[orf_map_dict['synonyms'].isin(cur_data.adata.obs['perturbed_target_symbol'].str.upper())]
    orf_map_dict = orf_map_dict.set_index('synonyms')['symbol'].to_dict()
    
    cur_data.replace_entries(
        slot='obs',
        column='perturbed_target_symbol',
        map_dict=orf_map_dict
    )
    
    cur_data.standardize_genes(
        slot='obs',
        input_column='perturbed_target_symbol',
        input_column_type='gene_symbol',
        multiple_entries=False
    )
    
    cur_data.count_entries(
        slot='obs',
        input_column='perturbed_target_symbol',
        count_column_name='perturbed_target_number',
        sep='|'
    )
    
    cur_data.create_columns(
        slot="obs",
        col_dict={
            "treatment_label": screen_metadata_dict['CONDITION_NAME']
        }
    )
    
    display(screen_df)

# curate_dataset(data_dir / mageck_data_files[0], metadata_df, non_curated_dir='CRISPR/non_curated/h5ad')

# %%

# replace the biorxiv link with doi for author mapping
metadata_df['SOURCE_ID'] = metadata_df['SOURCE_ID'].str.replace('https://www.biorxiv.org/content/10.1101/2020.08.27.270819v1.full','10.1101/2020.08.27.270819')
metadata_df['SOURCE_ID'] = metadata_df['SOURCE_ID'].str.replace('10.1101/2023.01.23.525275','10.1016/j.celrep.2023.112987')

# get pubmed info
pm_info = get_pubmed_info(list(metadata_df['SOURCE_ID'].unique()))

# manually replace 10.1101/2024.05.08.593110
pm_info['10.1101/2024.05.08.593110'] = {
    "title": 'TMEM106B-mediated SARS-CoV-2 infection allows for robust ACE2-independent infection in vitro but not in vivo',
    "year": 2024,
    "abstract": "Angiotensin converting enzyme 2 (ACE2) serves as the primary entry receptor for severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2). However, ACE2-independent entry has been observed in vitro for SARS-CoV-2 strains containing the E484D amino acid substitution in the spike protein. In this study, we conducted a whole genome CRISPR-Cas9 knockout screen using a SARS-CoV-2 strain containing the spike-E484D substitution (SARS-CoV-2MA1) to identify the ACE2-independent entry mechanisms. Our findings revealed that SARS-CoV-2MA1 infection in HEK293T cells relied on heparan sulfate and endocytic pathways, with TMEM106B emerging as the most significant contributor. While SARS-CoV-2MA1 productively infected human brain organoids and K18-hACE2 mouse brains, it did not infect C57BL/6J or Ifnar-/- mouse brains. This suggests that ACE2-independent entry via TMEM106B, which is a protein that is predominantly expressed in the brain, did not overtly increase the risk of SARS-CoV-2 neuroinvasiveness in wild-type mice. Importantly, SARS-CoV-2MA1 did not replicate in Ace2-/- mouse respiratory tracts. Overall, this suggests that robust ACE2-independent infection by SARS-CoV-2E484D is likely a phenomenon specific to in vitro conditions, with no apparent clinical implications.",
    "doi": "https://doi.org/10.1101/2024.05.08.593110",
    "authors": [{'LastName': 'Tau', 'ForeName': 'Kexin'}, {'LastName': 'Rawle', 'ForeName': 'Daniel J'}]
}


metadata_df_curated = metadata_df[['#SCREEN_ID', 'SCREEN_NAME', 'SOURCE_ID']]

metadata_df_curated['F_AUTHOR_F_NAME'] = metadata_df_curated['SOURCE_ID'].map(
    {k: v['authors'][0]['ForeName'] for k, v in pm_info.items()}
    )
metadata_df_curated['F_AUTHOR_L_NAME'] = metadata_df_curated['SOURCE_ID'].map(
    {k: v['authors'][0]['LastName'] for k, v in pm_info.items()}
    )
metadata_df_curated['L_AUTHOR_F_NAME'] = metadata_df_curated['SOURCE_ID'].map(
    {k: v['authors'][-1]['ForeName'] for k, v in pm_info.items()}
    )
metadata_df_curated['L_AUTHOR_L_NAME'] = metadata_df_curated['SOURCE_ID'].map(
    {k: v['authors'][-1]['LastName'] for k, v in pm_info.items()}
    )

metadata_df_curated['YEAR'] = metadata_df_curated['SOURCE_ID'].map(
    {k: v['year'] for k, v in pm_info.items()}
    )

metadata_df_curated['TITLE'] = metadata_df_curated['SOURCE_ID'].map(
    {k: v['title'] for k, v in pm_info.items()}
    )

metadata_df_curated['ABSTRACT'] = metadata_df_curated['SOURCE_ID'].map(
    {k: v['abstract'] for k, v in pm_info.items()}
    )

metadata_df_curated['DOI'] = metadata_df_curated['SOURCE_ID'].map(
    {k: v['doi'] for k, v in pm_info.items()}
    )

# screen format
metadata_df_curated['library_format_label'] = metadata_df['SCREEN_FORMAT'].map(
    {'Pool':'pooled', 'Array':'arrayed', 'in vivo':'in vivo'}
)

# software analysis
metadata_df_curated['software_analysis_label'] = metadata_df['ANALYSIS']

# library perturbation type
metadata_df_curated['library_perturbation_type_label'] = metadata_df['LIBRARY_METHODOLOGY'].map(
    {'Knockout':'knockout', 'Inhibition':'inhibition', 'Activation':'activation'}
)

# perturbation type
metadata_df_curated['perturbation_type_label'] = metadata_df['LIBRARY_TYPE'].map(
    {'CRISPRn': 'CRISPRn', 
     'CRISPRi': 'CRISPRi', 
     'CRISPRa': 'CRISPRa', 
     'Cytosine Base Editing-Mediated Gene KnockOut': 'CBE'}
)

# library name and manufacturer
metadata_df['LIBRARY'] = metadata_df['LIBRARY'].str.replace('(Gordon 2020)', '(Gordon, 2020)')
metadata_df_curated['library_name'] = metadata_df['LIBRARY']


# %%
import re
l = metadata_df['LIBRARY'].unique()

lib_manufact_map_dict = {}
lib_name_map_dict = {}

for e in l:
    match_res = re.search(r'\(([^()]+?,\s*\d{4})\)', e, re.IGNORECASE) # detect year - if detected, library manufacturer is present
    if match_res:
        match_str = match_res.group()
        match_str = match_str.lstrip('(').rstrip(')')
        manufact = match_str.split(',')[0]
        name = e.split(manufact)[0].rstrip('(').strip()
        lib_manufact_map_dict[e] = manufact
        lib_name_map_dict[e] = name
        # print(match_str)
    else: # if no author provided then the manufacturer is not provided and only the name is available
        lib_name_map_dict[e] = e
        lib_manufact_map_dict[e] = None
        # print(e)

metadata_df_curated['library_manufacturer'] = metadata_df_curated['library_name'].map(lib_manufact_map_dict)
metadata_df_curated['library_name'] = metadata_df_curated['library_name'].map(lib_name_map_dict)

# %%

# library  generation method

metadata_df['ENZYME'].unique()

enzyme_map_dict = {
    'Cas9':'SpCas9',
    'dCas9-KRAB': 'dCas9-KRAB',
    'sunCas9':'sunCas9',
    'SAM (NLS-dCas9-VP64/MS2-p65-HSF1)':'dCas9-SAM', 
    'AsCpf1':'AsCpf1', 
    'Cas9-v1':'SpCas9',
    'dCas9-VP64 + PP7-P65-HSF1':'dCas9-VP64', 
    'dCas9-BFP-KRAB':'dCas9-KRAB',
    'dCas9-SunTag-P2A-HygR':'dCas9-Suntag', 
    'dCAS-VP64_Blast (Zhu, 2021)':'dCas9-VP64',
    'AncBE4max':'AncBE4max', 
    'dCas9-VP64':'dCas9-VP64',
    'dCas9–VPR':'dCas9–VPR', 
    'Cas12a':'Cas12a',
    'dCas9-VP64 & p65-HSF1 (CRISPR SAM)':'dCas9-SAM', 
    'dCas9':'dCas9', 
    'iCAS9':'SpCas9'
}

metadata_df_curated['library_generation_method_label'] = metadata_df['ENZYME'].map(enzyme_map_dict)



# %%

pm_info = get_pubmed_info(list(metadata_df['SOURCE_ID'].unique()))










