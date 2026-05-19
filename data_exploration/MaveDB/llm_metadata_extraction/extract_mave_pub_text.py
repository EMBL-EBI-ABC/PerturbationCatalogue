import os
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from threading import Lock
from time import sleep

import pymupdf4llm
import requests
from paperscraper.pdf import save_pdf
from tqdm import tqdm

from xml_parser import xml_to_md

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]

MAVEDB_DUMP_DIR = REPO_ROOT / 'data_exploration' / 'MaveDB' / 'Dump' / 'mavedb-dump.20250612164404' / 'csv'
MAVEDB_METADATA_OUTPUT_DIR = REPO_ROOT / 'data_exploration' / 'MaveDB' / 'llm_metadata_extraction' / 'mavedb_metadata'
MAVEDB_URN_TO_DOIS_OUTPUT_FILE = REPO_ROOT / 'data_exploration' / 'MaveDB' / 'llm_metadata_extraction' / 'mavedb_urn_to_dois.json'
MAVEDB_API_BASE_URL = "https://api.mavedb.org/api/v1"
PAPERSCRAPER_API_KEYS_FILE = SCRIPT_DIR / 'scraper_api_keys.txt'
PAPERSCRAPER_DUMP_DIR = SCRIPT_DIR / 'paperscraper_dumps'
PAPERSCRAPER_FULL_TEXT_RAW_DIR = SCRIPT_DIR / 'pub_full_text_raw'
FULL_TEXT_MD_DIR = SCRIPT_DIR / 'pub_full_text_md'
DOWNLOAD_PROGRESS_LOG_FILE = SCRIPT_DIR / 'pub_full_text_download.log'
DEFAULT_DOWNLOAD_MAX_WORKERS = min(32, (os.cpu_count() or 1) * 4)
PRINT_LOCK = Lock()


def _append_log_block(title: str, *lines: str) -> None:
    timestamp = datetime.now().isoformat(timespec='seconds')
    separator = "=" * 80
    DOWNLOAD_PROGRESS_LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    with DOWNLOAD_PROGRESS_LOG_FILE.open('a', encoding='utf-8') as log_file:
        log_file.write(f"\n[{timestamp}] {separator}\n")
        log_file.write(f"[{timestamp}] {title}\n")
        for line in lines:
            log_file.write(f"[{timestamp}] {line}\n")
        log_file.write(f"[{timestamp}] {separator}\n")


def _append_log_line(message: str) -> None:
    timestamp = datetime.now().isoformat(timespec='seconds')
    DOWNLOAD_PROGRESS_LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    with DOWNLOAD_PROGRESS_LOG_FILE.open('a', encoding='utf-8') as log_file:
        log_file.write(f"[{timestamp}] {message}\n")


def print_status_block(title: str, *lines: str, echo: bool = True) -> None:
    separator = "=" * 80
    with PRINT_LOCK:
        _append_log_block(title, *lines)
        if not echo:
            return
        print(f"\n{separator}")
        print(title)
        for line in lines:
            print(line)
        print(separator)


def remove_references_from_markdown(text: str) -> str:
    """Remove a trailing references section from converted Markdown when it is clearly present.

    The function is intentionally conservative: it only removes a section when it finds
    a Markdown heading for References/Bibliography in the latter half of the document and
    the following lines look like bibliographic entries.
    """
    reference_heading_re = re.compile(r"^\s{0,3}(#{1,6})\s*(references?|bibliography)\s*:??\s*$", re.IGNORECASE)
    markdown_heading_re = re.compile(r"^\s{0,3}(#{1,6})\s+\S")
    numbered_reference_re = re.compile(r"^\s{0,3}(?:\[\d+\]|\d+[.)])\s+")
    doi_re = re.compile(r"(?:doi:\s*|https?://(?:dx\.)?doi\.org/|\b10\.\d{4,9}/)\S*", re.IGNORECASE)
    reference_year_re = re.compile(r"\(\d{4}[a-z]?\)")
    section_terminator_re = re.compile(r"^\s*\*\*(?:==|--)")

    if not text:
        return text

    lines = text.splitlines()
    if not lines:
        return text

    candidate_index = None
    candidate_heading_level = None
    for index, line in enumerate(lines):
        match = reference_heading_re.match(line)
        if match:
            candidate_index = index
            candidate_heading_level = len(match.group(1))

    if candidate_index is None or candidate_index < len(lines) // 2:
        return text

    end_index = len(lines)
    for index in range(candidate_index + 1, len(lines)):
        if section_terminator_re.match(lines[index]):
            end_index = index
            break

        match = markdown_heading_re.match(lines[index])
        if match and len(match.group(1)) <= candidate_heading_level:
            end_index = index
            break

    section_lines = lines[candidate_index + 1:end_index]
    non_empty_lines = [line.strip() for line in section_lines if line.strip()]
    if len(non_empty_lines) < 3:
        return text

    reference_signal_count = 0
    for line in non_empty_lines[:20]:
        if numbered_reference_re.match(line):
            reference_signal_count += 2
            continue
        if doi_re.search(line):
            reference_signal_count += 1
        if reference_year_re.search(line) and "," in line:
            reference_signal_count += 1

    if reference_signal_count < 4:
        return text

    cleaned_lines = lines[:candidate_index] + lines[end_index:]
    return "\n".join(cleaned_lines).strip() + "\n"


def get_unique_mavedb_urns(dump_dir: str | Path) -> list[str]:
    """Extract unique file basenames from CSV files in the dump directory.
    
    Parameters:
    - dump_dir: The directory containing the MaveDB CSVs from the dump.
    Returns:
    - A set of unique MaveDB URNs extracted from the file names.
    """
    dump_dir = Path(dump_dir).resolve()
    unique_files = set()
    for file_path in dump_dir.iterdir():
        if file_path.suffix == '.csv':
            filename = file_path.name
            # Trim .scores.csv or .counts.csv
            cleaned = filename.replace('.scores.csv', '').replace('.counts.csv', '')
            cleaned = cleaned.replace('urn-mavedb-', 'urn:mavedb:')
            unique_files.add(cleaned)
    return list(unique_files)

def get_mavedb_record_type_from_urn(urn: str) -> str:
	"""Infer the MaveDB record type from the URN structure.
 
    Parameters:
    - urn: The MaveDB URN (e.g., 'urn:mavedb:00000034-b-1')
    Returns:
    - The inferred record type as a string (e.g., 'score-sets', 'experiments', 'experiment-sets')
    """
	# MaveDB uses progressively more specific identifiers separated by hyphens.
	# The number of segments tells us which REST collection to query.
	identifier = urn.replace("urn:mavedb:", "")
	segment_count = len(identifier.split("-"))
	if segment_count >= 3:
		return "score-sets"
	if segment_count == 2:
		return "experiments"
	if segment_count == 1:
		return "experiment-sets"

	raise ValueError(f"Unable to infer MaveDB record type from URN: {urn}")

def fetch_mavedb_entry(urn) -> dict:
    """Fetch entry data from MaveDB API for a given URN.
    
    Parameters:
    - urn: The MaveDB URN to fetch (e.g., 'urn:mavedb:00000034-b-1')
    Returns:    
    - A dictionary containing the entry data from MaveDB.
    """
    record_type = get_mavedb_record_type_from_urn(urn)
    url = f"{MAVEDB_API_BASE_URL}/{record_type}/{urn}"
    response = requests.get(url)
    response.raise_for_status()  # Raise an error for bad status codes
    return response.json()

def get_dois_from_mavedb_entry(entry, log: bool = True) -> list[str] | None:
    """Extract DOI from MaveDB entry if available.
    
    Parameters:
    - entry: A dictionary representing the MaveDB entry data.
    - log: Whether to emit status logging for DOI presence/absence.
    Returns:
    - A list of DOIs if found, otherwise None.
    """
    # MaveDB places publication identifiers in a few different fields depending on record type.
    experiment = entry.get("experiment") or {}
    doi_sources = (
        (entry.get("doiIdentifiers") or [], "identifier"),
        (entry.get("primaryPublicationIdentifiers") or [], "doi"),
        (entry.get("secondaryPublicationIdentifiers") or [], "doi"),
        (experiment.get("primaryPublicationIdentifiers") or [], "doi"),
        (experiment.get("secondaryPublicationIdentifiers") or [], "doi"),
    )
    dois = sorted(
        {
            doi
            for identifiers, field_name in doi_sources
            for identifier in identifiers
            if isinstance(identifier, dict)
            if (doi := identifier.get(field_name))
        }
    )
    if dois:
        if log:
            print_status_block(
                "MaveDB publication identifiers found",
                f"Count: {len(dois)}",
                f"DOIs: {', '.join(dois)}",
            )
        return dois
    urn = entry.get("urn", "<unknown URN>")
    if log:
        print_status_block(
            "No publication identifiers found",
            f"Entry: {urn}",
        )
    return None


def export_mavedb_urn_to_dois_json(
    mavedb_entries_dir: Path | str = MAVEDB_METADATA_OUTPUT_DIR,
    output_file: Path | str = MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
) -> dict[str, list[str]]:
    """Write a JSON mapping from MaveDB URNs to publication DOIs.

    Parameters:
    - mavedb_entries_dir: Directory containing cached MaveDB entry JSON files.
    - output_file: Destination JSON file for the URN to DOI mapping.
    Returns:
    - A dictionary mapping each URN with at least one DOI to a sorted DOI list.
    """
    if isinstance(mavedb_entries_dir, (str, Path)):
        mavedb_entries_dir = Path(mavedb_entries_dir).resolve()
    output_file = Path(output_file).resolve()

    if not mavedb_entries_dir.is_dir():
        raise ValueError(f"MaveDB entries directory not found: {mavedb_entries_dir}")

    urn_to_dois: dict[str, list[str]] = {}
    entry_files = sorted(mavedb_entries_dir.glob("*.json"))
    for entry_file in tqdm(entry_files, desc="Building URN to DOI mapping", unit="entry"):
        try:
            entry = json.loads(entry_file.read_text(encoding='utf-8'))
            urn = entry.get("urn")
            if not urn:
                continue
            dois = get_dois_from_mavedb_entry(entry, log=False)
            if dois:
                urn_to_dois[urn] = dois
        except Exception as exc:
            print_status_block(
                "Error exporting URN to DOI mapping",
                f"File: {entry_file}",
                f"Error: {exc}",
            )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(json.dumps(urn_to_dois, indent=2, sort_keys=True), encoding='utf-8')

    print_status_block(
        "MaveDB URN to DOI mapping written",
        f"Entries scanned: {len(entry_files)}",
        f"URNs with DOIs: {len(urn_to_dois)}",
        f"Output file: {output_file}",
    )
    return urn_to_dois

def retrieve_pub_full_text(
    doi: str,
    api_keys_file: str | Path = PAPERSCRAPER_API_KEYS_FILE,
    output_dir: str | Path = PAPERSCRAPER_FULL_TEXT_RAW_DIR,
    overwrite: bool = False,
) -> Path | None:
    """Retrieve the full text of a publication given its DOI using paperscraper.
    
    Parameters:
    - doi: The DOI of the publication to retrieve.
    - api_keys_file: Path to a file containing API keys for paperscraper.
    - output_dir: Directory where the retrieved PDF should be saved.
    - overwrite: Whether to overwrite an existing downloaded file.
    
    Returns:
    - The file path to the saved PDF of the publication.
    """
    api_keys_file = Path(api_keys_file).resolve()
    output_dir = Path(output_dir).resolve()

    if not api_keys_file.is_file():
        raise ValueError(f"API keys file not found: {api_keys_file}")
    output_dir.mkdir(parents=True, exist_ok=True)
    paper_data = {'doi': doi}
    # Replace / and . in DOI with _ for the filename
    filename_base = re.sub(r'[/.]', '_', doi)
    filepath_pdf = output_dir / f"{filename_base}.pdf"
    filepath_xml = output_dir / f"{filename_base}.xml"

    if not overwrite:
        if filepath_pdf.is_file():
            print_status_block(
                "Publication full text already present",
                f"DOI: {doi}",
                "Format: PDF",
                f"Using existing file: {filepath_pdf}",
            )
            return filepath_pdf
        if filepath_xml.is_file():
            print_status_block(
                "Publication full text already present",
                f"DOI: {doi}",
                "Format: XML",
                f"Using existing file: {filepath_xml}",
            )
            return filepath_xml
    
    out = save_pdf(paper_data, api_keys=str(api_keys_file), filepath=str(filepath_pdf)) # returns True if successful, False otherwise 
    # despite name save_pdf, the function might save xml, if PDF is not available. We should check the file type and return correct file path.
    if out:
        # Check if the saved file is a PDF
        if filepath_pdf.is_file():
            print_status_block(
                "Publication full text successfully retrieved",
                f"DOI: {doi}",
                f"Format: PDF",
                f"Saved to: {filepath_pdf}",
            )
            return filepath_pdf
        else:
            if filepath_xml.is_file():
                print_status_block(
                    "Publication full text successfully retrieved",
                    f"DOI: {doi}",
                    f"Format: XML",
                    f"Saved to: {filepath_xml}",
                )
                return filepath_xml
    print_status_block(
        "Publication full text retrieval failed",
        f"DOI: {doi}",
        "No file was saved.",
    )
    return None

def pdf_to_md(
    filepath_pdf: str | Path,
    output_dir: str | Path = FULL_TEXT_MD_DIR,
    remove_references: bool = True,
) -> Path:
    """Extract text from a PDF file and save it as Markdown.
    Parameters:
    - filepath_pdf: The path to the PDF file to extract text from.
    - output_dir: The directory where the Markdown file should be saved.
    - remove_references: Whether to remove a trailing references section when confidently detected.
    Returns:
    - The file path to the saved Markdown file.
    """
    filepath_pdf = Path(filepath_pdf).resolve()
    output_dir = Path(output_dir).resolve()
    if not filepath_pdf.is_file():
        raise ValueError(f"PDF file not found: {filepath_pdf}")
    output_dir.mkdir(parents=True, exist_ok=True)
    filename_base = filepath_pdf.stem
    filepath_md = output_dir / f"{filename_base}.md"
    
    print_status_block(
        "Converting PDF to Markdown",
        f"Input PDF: {filepath_pdf}"
    )
    
    md = pymupdf4llm.to_markdown(str(filepath_pdf), footer=False, header=False)
    if remove_references:
        md = remove_references_from_markdown(md)
    
    filepath_md.write_text(md, encoding='utf-8')
        
    print_status_block(
        "PDF successfully converted to Markdown",
        f"Output Markdown: {filepath_md}"
    )
    
    return filepath_md

def convert_pub_full_text_to_md(
    input_filepath: str | Path,
    output_dir: str | Path = FULL_TEXT_MD_DIR,
    remove_references: bool = True,
) -> Path:
    """Convert a publication full text file (PDF or XML) to Markdown format.
    Parameters:
    - input_filepath: The path to the publication full text file (PDF or XML).
    - output_dir: The directory where the Markdown file should be saved.
    - remove_references: Whether to remove a trailing references section when confidently detected (only applicable to PDF input).
    Returns:
    - The file path to the saved Markdown file.
    """
    input_filepath = Path(input_filepath).resolve()
    output_dir = Path(output_dir).resolve()
    if not input_filepath.is_file():
        raise ValueError(f"File not found: {input_filepath}")
    output_dir.mkdir(parents=True, exist_ok=True)

    if input_filepath.suffix.lower() == ".pdf":
        return pdf_to_md(input_filepath, output_dir=output_dir, remove_references=remove_references)
    if input_filepath.suffix.lower() == ".xml":
        return xml_to_md(input_filepath, output_dir=output_dir, remove_references=remove_references)
    raise ValueError(f"Unsupported file type: {input_filepath.suffix}")

def bulk_fetch_mavedb_entries(
    input_files: list[str],
    output_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
    overwrite: bool = False,
    sleep_time: float = 0.1
) -> list[dict]:
    """
    Fetch MaveDB entries for a list of unique URNs, with optional sleep between requests.
    Optionally saves each entry as a JSON file in the specified output directory.
    If the file is already present, it will be skipped to avoid redundant API calls.
    
    Parameters:    
    - input_files: A list of unique MaveDB URNs to fetch entries for.
    - output_dir: The directory where the fetched entries should be saved (optional).
    - overwrite: Whether to overwrite existing files in the output directory (default: False). If False, existing files will be used and API calls will be skipped for those entries.
    - sleep_time: The number of seconds to sleep between requests to avoid overwhelming the API.
    Returns:    
    - A list of dictionaries containing the fetched MaveDB entries.
    """
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    for urn in tqdm(input_files, desc="Fetching MaveDB entries", unit="entry"):
        entry_output_path = output_dir / f"{urn.replace(':', '_')}.json"
        try:
            if entry_output_path.is_file() and not overwrite:
                print_status_block(
                    "MaveDB entry already present",
                    f"URN: {urn}",
                    f"Skipping download for file: {entry_output_path}",
                )
                entries.append(json.loads(entry_output_path.read_text(encoding='utf-8')))
                continue

            entry = fetch_mavedb_entry(urn)
            entries.append(entry)
            sleep(sleep_time)  # Sleep between requests
            entry_output_path.write_text(json.dumps(entry, indent=2), encoding='utf-8')
        except Exception as exc:
            print_status_block(
                "Error fetching MaveDB entry",
                f"URN: {urn}",
                f"Error: {exc}",
            )
    return entries

def bulk_download_pub_full_texts(
    mavedb_entries_dir: Path | str = MAVEDB_METADATA_OUTPUT_DIR,
    output_dir: str | Path = PAPERSCRAPER_FULL_TEXT_RAW_DIR,
    overwrite: bool = False,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
):
    """Bulk download publication full texts for a list of MaveDB entries.
    Parameters:
    - mavedb_entries_dir: A path to a directory containing JSON files of MaveDB entries.
    - output_dir: The directory where the downloaded full texts should be saved (optional).
    - overwrite: Whether to overwrite existing files in the output directory (default: False). If False, existing files will be used and downloads will be skipped for those entries.
    - max_workers: Number of worker threads used to download unique publication DOIs concurrently.
    Returns:
    - A list of file paths to the downloaded publication full texts.
    """
    if isinstance(mavedb_entries_dir, (str, Path)):
        mavedb_entries_dir = Path(mavedb_entries_dir).resolve()
    output_dir = Path(output_dir).resolve()
    if not mavedb_entries_dir.is_dir():
        raise ValueError(f"MaveDB entries directory not found: {mavedb_entries_dir}")
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")
    output_dir.mkdir(parents=True, exist_ok=True)

    doi_to_urns: dict[str, set[str]] = {}
    doi_order: list[str] = []
    doi_reference_count = 0
    entries_with_dois = 0
    entry_files = list(mavedb_entries_dir.glob("*.json"))
    for entry_file in tqdm(entry_files, desc="Collecting publication DOIs", unit="entry"):
        try:
            entry = json.loads(entry_file.read_text(encoding='utf-8'))
            urn = entry.get("urn", "<unknown URN>")
            dois = get_dois_from_mavedb_entry(entry, log=False)
            if dois:
                entries_with_dois += 1
                doi_reference_count += len(dois)
                for doi in dois:
                    if doi not in doi_to_urns:
                        doi_order.append(doi)
                    doi_to_urns.setdefault(doi, set()).add(urn)
        except Exception as exc:
            print_status_block(
                "Error processing MaveDB entry file",
                f"File: {entry_file}",
                f"Error: {exc}",
            )

    print_status_block(
        "Publication DOI collection complete",
        f"Entries scanned: {len(entry_files)}",
        f"Entries with DOIs: {entries_with_dois}",
        f"DOI references found: {doi_reference_count}",
        f"Unique DOIs found: {len(doi_order)}",
        f"Progress log: {DOWNLOAD_PROGRESS_LOG_FILE}",
    )

    results_by_doi: dict[str, Path] = {}
    total_dois = len(doi_order)
    completed_dois = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_doi = {}
        for doi in doi_order:
            urns = doi_to_urns[doi]
            urns_str = ", ".join(sorted(urns))
            print_status_block(
                "Queueing publication full text retrieval",
                f"DOI: {doi}",
                f"MaveDB URNs: {urns_str}",
            )
            future = executor.submit(
                retrieve_pub_full_text,
                doi=doi,
                output_dir=output_dir,
                overwrite=overwrite,
            )
            future_to_doi[future] = doi

        for future in tqdm(as_completed(future_to_doi), total=len(future_to_doi), desc="Downloading publication full texts", unit="doi"):
            doi = future_to_doi[future]
            try:
                full_text_path = future.result()
                if full_text_path:
                    results_by_doi[doi] = full_text_path
                completed_dois += 1
                _append_log_line(
                    f"Download progress: {completed_dois}/{total_dois} DOIs processed; latest DOI: {doi}; status: {'ok' if full_text_path else 'no-file'}"
                )
            except Exception as exc:
                completed_dois += 1
                _append_log_line(
                    f"Download progress: {completed_dois}/{total_dois} DOIs processed; latest DOI: {doi}; status: error; error: {exc}"
                )
                print_status_block(
                    "Error retrieving publication full text",
                    f"DOI: {doi}",
                    f"Error: {exc}",
                )
    return [results_by_doi[doi] for doi in doi_order if doi in results_by_doi]

def bulk_convert_full_texts_to_md(
    input_files: list[str | Path],
    output_dir: Path,
    remove_references: bool = True,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
) -> list[Path]:
    """Bulk convert publication full text files (PDF or XML) to Markdown format.
    Parameters:
    - input_files: A list of file paths to the publication full text files (PDF or XML) to convert.
    - output_dir: The directory where the converted Markdown files should be saved.
    - remove_references: Whether to remove a trailing references section when confidently detected (only applicable to PDF input).
    - max_workers: Number of worker threads used to convert files concurrently.
    Returns:
    - A list of file paths to the converted Markdown files, in the same order as the input files (files that failed to convert will be skipped).
    """
    
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")

    files = [Path(file_path).resolve() for file_path in input_files]
    print_status_block(
        "Starting publication full text Markdown conversion",
        f"Files queued: {len(files)}",
        f"Output directory: {output_dir}",
        f"Remove references: {remove_references}",
        f"Max workers: {max_workers}",
    )

    converted_files_by_index: dict[int, Path] = {}
    completed_files = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {
            executor.submit(
                convert_pub_full_text_to_md,
                file_path,
                output_dir=output_dir,
                remove_references=remove_references,
            ): index
            for index, file_path in enumerate(files)
        }

        for future in tqdm(
            as_completed(future_to_index),
            total=len(future_to_index),
            desc="Converting publication full texts to Markdown",
            unit="file",
        ):
            index = future_to_index[future]
            file_path = files[index]
            try:
                md_file = future.result()
                completed_files += 1
                if md_file:
                    converted_files_by_index[index] = md_file
                _append_log_line(
                    f"Markdown conversion progress: {completed_files}/{len(files)} files processed; latest file: {file_path}; status: {'ok' if md_file else 'no-output'}"
                )
            except Exception as exc:
                completed_files += 1
                _append_log_line(
                    f"Markdown conversion progress: {completed_files}/{len(files)} files processed; latest file: {file_path}; status: error; error: {exc}"
                )
                print_status_block(
                    "Error converting publication full text to Markdown",
                    f"File: {file_path}",
                    f"Error: {exc}",
                )

    converted_files = [converted_files_by_index[index] for index in sorted(converted_files_by_index)]
    print_status_block(
        "Publication full text Markdown conversion complete",
        f"Files queued: {len(files)}",
        f"Files converted: {len(converted_files)}",
        f"Files failed: {len(files) - len(converted_files)}",
        f"Output directory: {output_dir}",
        f"Progress log: {DOWNLOAD_PROGRESS_LOG_FILE}",
    )
    return converted_files


def main() -> None:
    unique_files = get_unique_mavedb_urns(MAVEDB_DUMP_DIR)

    bulk_fetch_mavedb_entries(unique_files, output_dir=MAVEDB_METADATA_OUTPUT_DIR, overwrite=False)
    export_mavedb_urn_to_dois_json(mavedb_entries_dir=MAVEDB_METADATA_OUTPUT_DIR, output_file=MAVEDB_URN_TO_DOIS_OUTPUT_FILE)
    full_text_paths = bulk_download_pub_full_texts(
        mavedb_entries_dir=MAVEDB_METADATA_OUTPUT_DIR,
        output_dir=PAPERSCRAPER_FULL_TEXT_RAW_DIR,
        overwrite=False,
    )
    bulk_convert_full_texts_to_md(full_text_paths, output_dir=FULL_TEXT_MD_DIR, remove_references=True)


if __name__ == "__main__":
    main()


