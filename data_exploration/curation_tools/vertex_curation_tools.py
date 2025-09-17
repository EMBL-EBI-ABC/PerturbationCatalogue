import os
import json
import requests

from dotenv import load_dotenv

import vertexai
from vertexai.generative_models import GenerativeModel

from elsapy.elsclient import ElsClient
from elsapy.elsdoc import FullDoc

load_dotenv()

schema = {
    "perturbation_type_label": {
        "description": "Perturbation type ontology term label of the investigated sample.",
        "type": "string",
        "examples": ["CRISPRn", "CRISPRi", "CRISPRa"],
    },
    "timepoint": {
        "description": "Timepoint of the investigated sample in ISO 8601 format.",
        "type": "string",
        "examples": ["P1DT12H30M15S"],
    },
    "treatment_label": {
        "description": "Treatment/compound ontology term label used to stimulate the investigated sample. ChEMBL compound label.",
        "type": "string",
        "examples": ["dexamethasone", "tamoxifen", "cisplatin"],
    },
    "model_system_label": {
        "description": "Model system ontology term label of the investigated sample.",
        "type": "string",
        "examples": ["cell_line", "primary_cell", "organoid"],
    },
    "species": {
        "description": "Species name of the investigated sample.",
        "type": "string",
        "examples": ["Homo sapiens"],
    },
    "tissue_label": {
        "description": "Tissue ontology term label of the investigated sample. Must be part of the UBERON ontology.",
        "type": "string",
        "examples": ["brain", "liver", "lung", "heart", "blood"]
    },
    "cell_type_label": {
        "description": "Cell type ontology term label of the investigated sample. Must be part of the Cell Ontology (CL).",
        "type": "string",
        "examples": ["neurons", "astrocytes", "microglia", "endothelial cells", "fibroblasts"],
    },
    "cell_line_label": {
        "description": "Cell line ontology term label of the investigated sample. Must be part of the Cell Line Ontology (CLO).",
        "type": "string",
        "examples": ["HeLa", "HEK293", "A549", "MCF7", "RPE1", "HAP1"],
    },
    "sex_label": {
        "description": "Sex ontology term label of the investigated sample.",
        "type": "string",
        "examples": ["male", "female", "unknown"]
    },
    "developmental_stage_label": {
        "description": "Developmental stage ontology term label of the investigated sample.",
        "type": "string",
        "examples": [
            "embryonic",
            "fetal",
            "neonatal",
            "adolescent",
            "adult",
            "senior adult",
        ],
    },
    "disease_label": {
        "description": "Disease ontology term label of the investigated sample. Must be part of the MONDO ontology.",
        "type": "string",
        "examples": ["breast cancer", "lung cancer", "diabetes mellitus", "asthma", "alzheimer's disease", "covid-19"],
    },
    "study_title": {"description": "Title of the study/publication.", "type": "string"},
    "study_uri": {
        "description": "URI/DOI of the study/publication.",
        "type": "string",
        "format": "uri",
    },
    "study_year": {
        "description": "Publication year of the study/publication.",
        "type": "integer",
        "examples": [2020, 2021, 2022, 2023],
    },
    "first_author": {
        "description": "Full name of the first author of the study/publication.",
        "type": "string",
        "examples": ["John Doe", "Jane Smith"],
    },
    "last_author": {
        "description": "Full name of the last author of the study/publication.",
        "type": "string",
        "examples": ["Alice Johnson", "Bob Brown"],
    },
    "experiment_title": {"description": "Title of the experiment.", "type": "string"},
    "experiment_summary": {
        "description": "Summary of the experiment.",
        "type": "string",
    },
    "number_of_perturbed_targets": {
        "description": "Total number of perturbed targets in the experiment.",
        "type": "integer",
        "examples": [10, 20, 30, 20000],
    },
    "number_of_perturbed_samples": {
        "description": "Total number of perturbed samples/cells in the experiment.",
        "type": "integer",
        "examples": [100, 200, 300, 10000, 12000000],
    },
    "library_generation_type_label": {
        "description": "Library generation type ontology term label, defined in EFO under parent term EFO:0022867 (genetic perturbation). Endogenous perturbation methods directly modify the genome (e.g., CRISPR-Cas9), while exogenous perturbation methods introduce external genetic material (e.g., RNAi or introduction genetic variants within an exogenous construct).",
        "type": "string",
        "examples": ['endogenous perturbation method', 'exogenous perturbation method'],
        
    },
    "library_generation_method_label": {
        "description": "Library generation method ontology term label, defined in EFO under parent term EFO:0022868 (Endogenous genetic perturbation method)",
        "type": "string",
        "examples": ["SpCas9", "AsCas12a", "dCas9-KRAB", "dCas9-VPR", "AbeMax", "PE2", "PE3"],
    },
    "enzyme_delivery_method_label": {
        "description": "Enzyme delivery method ontology term label.",
        "type": "string",
        "examples": [
            "lipofection",
            "nucleofection",
            "retrovirus transduction",
            "lentivirus transduction",
        ],
    },
    "library_delivery_method_label": {
        "description": "Library delivery method ontology term label.",
        "type": "string",
        "examples": [
            "lipofection",
            "nucleofection",
            "retrovirus transduction",
            "lentivirus transduction",
        ],
    },
    "enzyme_integration_state_label": {
        "description": "Enzyme integration state ontology term label.",
        "type": "string",
        "examples": [
            "random locus integration",
            "targeted locus integration",
            "native locus replacement",
            "non-integrative transgene expression",
        ],
    },
    "library_integration_state_label": {
        "description": "Library integration state ontology term label.",
        "type": "string",
        "examples": [
            "random locus integration",
            "targeted locus integration",
            "native locus replacement",
            "non-integrative transgene expression",
        ],
    },
    "enzyme_expression_control_label": {
        "description": "Enzyme expression control ontology term label.",
        "type": "string",
        "examples": [
            "constitutive transgene expression",
            "inducible transgene expression",
            "native promoter-driven transgene expression",
            "degradation domain-based transgene control",
        ],
    },
    "library_expression_control_label": {
        "description": "Library expression control ontology term label.",
        "type": "string",
        "examples": [
            "constitutive transgene expression",
            "inducible transgene expression",
            "native promoter-driven transgene expression",
            "degradation domain-based transgene control",
        ],
    },
    "library_name": {
        "description": "Name of the perturbation library.",
        "type": "string",
        "examples": ["Bassik Human CRISPR Knockout Library", "Brunello", "GeCKO v2"],
    },
    "library_uri": {
        "description": "URI/accession of the perturbation library.",
        "type": "string",
        "format": "uri",
    },
    "library_format_label": {
        "description": "Perturbation library format ontology term label.",
        "type": "string",
        "examples": ["pooled", "arrayed", "in vivo"],
    },
    "library_scope_label": {
        "description": "Perturbation library scope ontology term label.",
        "type": "string",
        "examples": ["focused", "genome-wide"],
    },
    "library_perturbation_type_label": {
        "description": "Ontology term label for the library perturbation type.",
        "type": "string",
        "examples": [
            "knockout",
            "inhibition",
            "activation",
            "base editing",
            "prime editing",
        ],
    },
    "library_manufacturer": {
        "description": "Name of the library manufacturer/vendor/origin lab.",
        "type": "string",
        "examples": ["Bassik", "Weissman"],
    },
    "library_lentiviral_generation": {
        "description": "Generation number of the lentiviral library.",
        "type": "integer",
        "examples": [3],
    },
    "library_grnas_per_target": {
        "description": "Number of gRNAs per target.",
        "type": "string",
        "examples": ["4", "5-7"],
    },
    "library_total_grnas": {
        "description": "Total number of gRNAs in the library.",
        "type": "integer",
        "examples": [22, 100, 20000],
    },
    "library_total_variants": {
        "description": "Only for MAVE studies; Total number of variants in the library.",
        "type": "integer",
        "examples": [5000, 20000, 100000],
    },
    "readout_dimensionality_label": {
        "description": "Ontology term label associated with the dimensionality of the readout assay.",
        "type": "string",
        "examples": ["single-dimensional assay", "high-dimensional assay"],
    },
    "readout_type_label": {
        "description": "Ontology term label associated with the type of the readout assay.",
        "type": "string",
        "examples": ["transcriptomic", "proteomic", "phenotypic"],
    },
    "readout_technology_label": {
        "description": "Ontology term label associated with the technology used in the readout assay.",
        "type": "string",
        "examples": [
            "single-cell rna-seq",
            "population growth assay",
            "flow cytometry",
        ],
    },
    "method_name_label": {
        "description": "Ontology term label associated with the method name used in the readout assay.",
        "type": "string",
        "examples": ["Perturb-seq", "scRNA-seq", "proliferation CRISPR screen"],
    },
    "method_uri": {
        "description": "URI associated with the method used in the readout assay.",
        "type": "string",
        "format": "uri",
    },
    "sequencing_library_kit_label": {
        "description": "Ontology term label associated with the sequencing library kit.",
        "type": "string",
        "examples": [
            "10x Genomics Chromium GEM-X Single Cell 5-prime kit v3",
            "10x Genomics Single Cell 3-prime",
            "Nextera XT",
        ],
    },
    "sequencing_platform_label": {
        "description": "Ontology term label associated with the sequencing platform.",
        "type": "string",
        "examples": [
            "Illumina NovaSeq X Plus",
            "Illumina HiSeq 4000",
            "Illumina HiSeq 2500",
        ],
    },
    "sequencing_strategy_label": {
        "description": "Ontology term label associated with the sequencing strategy.",
        "type": "string",
        "examples": ["barcode sequencing", "shotgun sequencing"],
    },
    "software_counts_label": {
        "description": "Ontology term label for the software used for generating counts.",
        "type": "string",
        "examples": ["CellRanger", "Drop-seq Tools"],
    },
    "software_analysis_label": {
        "description": "Ontology term label for the software used for analysis.",
        "type": "string",
        "examples": ["custom", "MAGeCK"],
    },
    "reference_genome_label": {
        "description": "Ontology term label for the reference genome.",
        "type": "string",
        "examples": ["GRCh38", "GRCh37"],
    },
}

def intitialize_vertexai():
    """Initializes Vertex AI with project and location from environment variables."""
    try:
        vertexai.init(project=os.getenv("PROJECT_ID"), location=os.getenv("LOCATION"))
        print("Successfully connected to Google Cloud services.")
    except Exception as e:
        print(f"Error initializing Google Cloud services: {e}")
        print(
            "Please ensure you have run 'gcloud auth application-default login' in your terminal."
        )

def extract_text_with_elsapy(doi, save_dir=None):
    """Fetches full text of a document from Elsevier API using its DOI."""
    client = ElsClient(os.getenv("ELSEVIER_API_KEY"))
    try:
        doi_doc = FullDoc(doi=doi)
        if doi_doc.read(client):
            print("Successfully retrieved document from Elsevier API.")
            # save the full JSON data to a file
            if save_dir:
                if not os.path.exists(save_dir):
                    os.makedirs(save_dir, exist_ok=True)
                    print(f"Created directory: {save_dir}")
                doi_save_name = doi.replace('/', '_')
                with open(f"{save_dir}/{doi_save_name}.json", 'w') as f:
                    json.dump(doi_doc.data, f, indent=2)
                    print(f"Saved document data to {save_dir}/{doi_save_name}.json")
            # extract the original text
            original_text = doi_doc.data.get('originalText', None)
            if original_text:
                print(f"Extracted original text for doi: {doi}")
                # clean the text
                cleaned_text = clean_elsapy_text(original_text)
                return cleaned_text
            else:
                print("No original text found in the document data.")
                return None
        else:
            print(f"The paper doesn't exist or is not accessible: {doi}")
            return None
    except Exception as e:
        print(f"Error retrieving document: {e}")
        return None

def clean_elsapy_text(text):
    """Cleans the text extracted from Elsevier API."""
    if not text:
        return "Error: No text provided for cleaning."

    url_pattern = re.compile(r'https?://\S+|www\.\S+')
    docs_pattern = re.compile(r'\S+\.(jpg|png|gif|tif|sml|zip|docx|xlsx)', re.IGNORECASE | re.MULTILINE)

    clean_text = text
    # remove URLs
    clean_text = url_pattern.sub('', clean_text)
    # remove document links
    clean_text = docs_pattern.sub('', clean_text)

    # remove everything before 'Introduction'
    intro_pos = clean_text.find('Introduction')
    if intro_pos != -1:
        clean_text = clean_text[intro_pos:]

    # remove everything after 'References'
    ref_pos = clean_text.rfind('References')
    if ref_pos != -1:
        clean_text = clean_text[:ref_pos]

    # remove everything after 'Acknowledgements'
    ack_pos = clean_text.rfind('Acknowledgements')
    if ack_pos != -1:
        clean_text = clean_text[:ack_pos]

    print(
        f"Cleaned the extracted text.\nOriginal length: {len(text)} characters.\nNew length: {len(clean_text)} characters.\nReduction %: {100 - (len(clean_text) / len(text) * 100):.2f}%"
    )
    return clean_text


def fetch_pmc_json(pmc_id=None, save_dir=None):
    """
    Fetches the PMC JSON data for a given PMC ID.

    Args:
        pmc_id (str): The PMC ID of the article.
        save_dir (str): The directory to save the JSON file.

    Returns:
        dict: The JSON data retrieved from the PMC API or an error message.
    """
    if not pmc_id:
        raise ValueError("PMC ID must be provided.")
    
    url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmc_id}/unicode"
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad responses
        
        data = response.json()

        if save_dir:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
                print(f"Created directory: {save_dir}")
            # Save the JSON data to a file
            with open(f"{save_dir}/{pmc_id}.json", 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
                print(f"Saved JSON data to {save_dir}/{pmc_id}.json")
        return data
    except (json.JSONDecodeError, Exception) as e:
        print(f"Paper not in Open Access - Skipping PMID: {pmc_id}")
        return {"error": "Failed to decode JSON response.", "details": str(e)}
    except requests.HTTPError as e:
        print(f"Illegal PMID - Skipping PMC ID {pmc_id}")
        return {"error": "HTTP error occurred.", "details": str(e)}


def extract_pmc_full_text(filepath=None, save_dir=None):
    """
    Extracts the full text from a PMC API JSON response file.

    Args:
        filepath (str): The path to the JSON file.

    Returns:
        str: The concatenated full text from all passages.
    """
    if not filepath or not os.path.isfile(filepath):
        return "Error: A valid file path must be provided."
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # Navigate through the JSON structure to get to the passages
        # The structure is a list containing one dictionary
        documents = data[0].get('documents', [])
        if not documents:
            return "Error: 'documents' array not found or is empty."

        # The 'documents' array contains one dictionary
        passages = documents[0].get('passages', [])
        if not passages:
            return "Error: 'passages' array not found or is empty."

        # Loop through passages and collect relevant data
        full_text_parts = []
        for passage in passages:
            infon = passage.get('infons', {})
            section_type = infon.get('section_type', None)
            content_type = infon.get('type', None)
            text = passage.get('text', '')
            if section_type == 'TITLE':
                full_text_parts.append(f"# {text}\n")
                full_text_parts.append(f"DOI: {infon.get('article-id_doi', 'N/A')}\n")
                
                names = {k: v for k, v in infon.items() if k.startswith('name_')}
                first_author = names["name_0"].split(';')
                first_author = f"{first_author[1].split(':')[1]} {first_author[0].split(':')[1]}"
                last_author = names[f"name_{len(names)-1}"].split(';')
                last_author = f"{last_author[1].split(':')[1]} {last_author[0].split(':')[1]}"
                
                full_text_parts.append(f"First Author: {first_author}\n")
                full_text_parts.append(f"Last Author: {last_author}\n")
                
                full_text_parts.append(f"Year: {infon.get('year', 'N/A')}\n")
            
            else:
                if content_type.startswith('title_1'):
                    full_text_parts.append(f"# {text}\n")
                elif content_type.startswith('title_2'):
                    full_text_parts.append(f"## {text}\n")
                elif content_type.startswith('title_3'):
                    full_text_parts.append(f"### {text}\n")
                elif content_type.startswith('title_4'):
                    full_text_parts.append(f"#### {text}\n")
                elif content_type.startswith('fig'):
                    fig_id = infon.get('id', 'N/A')
                    if content_type == 'fig_title_caption':
                        full_text_parts.append(f"*Figure {fig_id}: {text}*\n")
                    elif content_type == 'fig_caption':
                        full_text_parts.append(f"{text}\n")
                elif content_type.startswith('ref'):
                    continue
                else:
                    full_text_parts.append(f"{text}\n")
                    
        # Join all the text parts together with a newline for readability
        full_text = '\n'.join(full_text_parts)
        
        if save_dir:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
                print(f"Created directory: {save_dir}")
            base_filename = os.path.basename(filepath).replace('.json', '.txt')
            output_path = os.path.join(save_dir, base_filename)
            with open(output_path, 'w', encoding='utf-8') as out_f:
                out_f.write(full_text)
                print(f"Saved extracted text to {output_path}")
                
        return full_text

    except FileNotFoundError:
        return f"Error: The file '{filepath}' was not found."
    except json.JSONDecodeError:
        return f"Error: The file '{filepath}' is not a valid JSON file."
    except (IndexError, KeyError) as e:
        return f"Error: The JSON structure is not as expected. Missing key: {e}"
    except Exception as e:
        return f"An unexpected error occurred: {e}"


def order_result_dict(input_dict_list, schema):
    """Orders the result nested dictionary according to the schema keys."""
    output_dict_list = []
    for input_dict in input_dict_list:
        output_dict = {}
        for key in schema.keys():
            if key in input_dict:
                output_dict[key] = input_dict[key]
            else:
                output_dict[key] = None 
                print(f"Warning: Key '{key}' not found in input dictionary. Setting it to None.")
        output_dict_list.append(output_dict)
    return output_dict_list

def cleanup_response(text):
    """Cleans up the model response text."""
    response_text = text.text
    response_text = response_text.strip().replace("```json", "").replace("```", "").strip()
    response_text = json.loads(response_text)
    ordered_result_dict = order_result_dict(response_text, schema)
    return ordered_result_dict

def extract_data_with_gemini(text_content, schema=schema, model_name=os.getenv("MODEL_NAME")):
    """Sends text to Gemini and asks for structured dictionary extraction."""

    model = GenerativeModel(model_name)
    # This detailed prompt guides the model to return clean JSON that matches your schema.
    prompt = f"""
    You are a an experienced biocurator and data analyst. 
    Your task is to carefully read the provided scientific papers and extract relevant metadata that will be used for a new data repository containing genetic perturbation studies, such as MAVE, CRISPR and Perturb-seq studies. 
    The data needs to be structured as a python dict that adheres to the provided JSON schema.
    You need to be critical and meticulous in your data extraction efforts. 
    When there is no information available in the paper for a particular field - leave the fields empty. 
    When there are more than one experiment reported in the paper, you will create a list of python dict structures - one for each experiment.

    **Instructions:**
    1.  Carefully read the text.
    2.  Extract the data points according to the JSON schema below.
    3.  If a value is not found, use `null` or an empty array `[]` as appropriate.
    4.  Your entire response must be ONLY the list of dict objects, with no introductory text, explanations, or markdown formatting like ```json.

    **JSON schema:**
    
    {schema}
    
    **Scientific Paper Text:**
    ---
    {text_content}
    ---
    """

    print("Constructed the prompt for with custom instructions, schema and scientific paper text.")
    print(f"Prompt length: {len(prompt)} characters.")
    print("Sending the prompt to Gemini model...")
    try:
        response = model.generate_content(prompt)
        clean_response = cleanup_response(response)
        print("Received and cleaned up the response from Gemini model.")
        return clean_response
    except Exception as e:
        print(f"Error calling Gemini model: {e}")
        return None