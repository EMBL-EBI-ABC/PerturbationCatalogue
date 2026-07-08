import json
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pymupdf4llm
from paperscraper.pdf import save_pdf
from tqdm import tqdm

from curation_tools.llm_curation.logging_utils import (
    append_log_line,
    print_status_block,
)
from curation_tools.llm_curation.xml_parser import xml_to_md

LLM_CURATION_DIR = Path(__file__).resolve().parent
PAPERSCRAPER_DUMP_DIR = LLM_CURATION_DIR / "paperscraper_dumps"
PAPERSCRAPER_FULL_TEXT_RAW_DIR = LLM_CURATION_DIR / "mavedb" / "pub_full_text_raw"
FULL_TEXT_MD_DIR = LLM_CURATION_DIR / "mavedb" / "pub_full_text_md"
DOWNLOAD_PROGRESS_LOG_FILE = LLM_CURATION_DIR / "pub_full_text_download.log"
DEFAULT_DOWNLOAD_MAX_WORKERS = min(32, (os.cpu_count() or 1) * 4)


def retrieve_pub_full_text(
    doi: str,
    output_dir: str | Path = PAPERSCRAPER_FULL_TEXT_RAW_DIR,
    overwrite: bool = False,
) -> Path | None:
    """Retrieve the full text of a publication given its DOI using paperscraper."""
    output_dir = Path(output_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    paper_data = {"doi": doi}
    filename_base = re.sub(r"[/.]", "_", doi)
    filepath_pdf = output_dir / f"{filename_base}.pdf"
    filepath_xml = output_dir / f"{filename_base}.xml"

    if not overwrite:
        if filepath_pdf.is_file():
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "Publication full text already present",
                f"DOI: {doi}",
                "Format: PDF",
                f"Using existing file: {filepath_pdf}",
            )
            return filepath_pdf
        if filepath_xml.is_file():
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "Publication full text already present",
                f"DOI: {doi}",
                "Format: XML",
                f"Using existing file: {filepath_xml}",
            )
            return filepath_xml

    out = save_pdf(paper_data, filepath=str(filepath_pdf))
    if out:
        if filepath_pdf.is_file():
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "Publication full text successfully retrieved",
                f"DOI: {doi}",
                "Format: PDF",
                f"Saved to: {filepath_pdf}",
            )
            return filepath_pdf
        if filepath_xml.is_file():
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "Publication full text successfully retrieved",
                f"DOI: {doi}",
                "Format: XML",
                f"Saved to: {filepath_xml}",
            )
            return filepath_xml
    print_status_block(
        DOWNLOAD_PROGRESS_LOG_FILE,
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
    """Extract text from a PDF file and save it as Markdown."""
    filepath_pdf = Path(filepath_pdf).resolve()
    output_dir = Path(output_dir).resolve()
    if not filepath_pdf.is_file():
        raise ValueError(f"PDF file not found: {filepath_pdf}")
    output_dir.mkdir(parents=True, exist_ok=True)
    filepath_md = output_dir / f"{filepath_pdf.stem}.md"

    print_status_block(
        DOWNLOAD_PROGRESS_LOG_FILE,
        "Converting PDF to Markdown",
        f"Input PDF: {filepath_pdf}",
    )

    markdown_text = pymupdf4llm.to_markdown(
        str(filepath_pdf), footer=False, header=False
    )
    if remove_references:
        markdown_text = remove_references_from_markdown(markdown_text)

    filepath_md.write_text(markdown_text, encoding="utf-8")

    print_status_block(
        DOWNLOAD_PROGRESS_LOG_FILE,
        "PDF successfully converted to Markdown",
        f"Output Markdown: {filepath_md}",
    )
    return filepath_md


def convert_pub_full_text_to_md(
    input_filepath: str | Path,
    output_dir: str | Path = FULL_TEXT_MD_DIR,
    remove_references: bool = True,
) -> Path:
    """Convert a publication full text file (PDF or XML) to Markdown format."""
    input_filepath = Path(input_filepath).resolve()
    output_dir = Path(output_dir).resolve()
    if not input_filepath.is_file():
        raise ValueError(f"File not found: {input_filepath}")
    output_dir.mkdir(parents=True, exist_ok=True)

    if input_filepath.suffix.lower() == ".pdf":
        return pdf_to_md(
            input_filepath, output_dir=output_dir, remove_references=remove_references
        )
    if input_filepath.suffix.lower() == ".xml":
        return xml_to_md(
            input_filepath, output_dir=output_dir, remove_references=remove_references
        )
    raise ValueError(f"Unsupported file type: {input_filepath.suffix}")


def remove_references_from_markdown(text: str) -> str:
    """Remove a trailing references section from converted Markdown when clearly present."""
    reference_heading_re = re.compile(
        r"^\s{0,3}(#{1,6})\s*(references?|bibliography)\s*:??\s*$", re.IGNORECASE
    )
    markdown_heading_re = re.compile(r"^\s{0,3}(#{1,6})\s+\S")
    numbered_reference_re = re.compile(r"^\s{0,3}(?:\[\d+\]|\d+[.)])\s+")
    doi_re = re.compile(
        r"(?:doi:\s*|https?://(?:dx\.)?doi\.org/|\b10\.\d{4,9}/)\S*", re.IGNORECASE
    )
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

    section_lines = lines[candidate_index + 1 : end_index]
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


def bulk_download_pub_full_texts(
    doi_to_urns: dict[str, set[str]],
    output_dir: str | Path = PAPERSCRAPER_FULL_TEXT_RAW_DIR,
    doi_to_fulltext_output_file: str | Path | None = None,
    overwrite: bool = False,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
) -> list[Path]:
    """Bulk download publication full texts for a mapping of DOIs to MaveDB URNs."""
    output_dir = Path(output_dir).resolve()
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")
    output_dir.mkdir(parents=True, exist_ok=True)

    doi_order = list(doi_to_urns)
    results_by_doi: dict[str, Path | None] = {doi: None for doi in doi_order}
    total_dois = len(doi_order)
    completed_dois = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_doi = {}
        for doi in doi_order:
            urns = doi_to_urns[doi]
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "Queueing publication full text retrieval",
                f"DOI: {doi}",
                f"MaveDB URNs: {', '.join(sorted(urns))}",
            )
            future = executor.submit(
                retrieve_pub_full_text,
                doi=doi,
                output_dir=output_dir,
                overwrite=overwrite,
            )
            future_to_doi[future] = doi

        for future in tqdm(
            as_completed(future_to_doi),
            total=len(future_to_doi),
            desc="Downloading publication full texts",
            unit="doi",
        ):
            doi = future_to_doi[future]
            try:
                full_text_path = future.result()
                results_by_doi[doi] = full_text_path
                completed_dois += 1
                append_log_line(
                    DOWNLOAD_PROGRESS_LOG_FILE,
                    f"Download progress: {completed_dois}/{total_dois} DOIs processed; latest DOI: {doi}; status: {'ok' if full_text_path else 'no-file'}",
                )
            except Exception as exc:
                completed_dois += 1
                append_log_line(
                    DOWNLOAD_PROGRESS_LOG_FILE,
                    f"Download progress: {completed_dois}/{total_dois} DOIs processed; latest DOI: {doi}; status: error; error: {exc}",
                )
                print_status_block(
                    DOWNLOAD_PROGRESS_LOG_FILE,
                    "Error retrieving publication full text",
                    f"DOI: {doi}",
                    f"Error: {exc}",
                )
    if doi_to_fulltext_output_file is not None:
        doi_to_fulltext_output_file = Path(doi_to_fulltext_output_file).resolve()
        doi_to_fulltext_output_file.parent.mkdir(parents=True, exist_ok=True)
        doi_to_fulltext = {
            doi: str(full_text_path) if full_text_path is not None else None
            for doi, full_text_path in results_by_doi.items()
        }
        doi_to_fulltext_output_file.write_text(
            json.dumps(doi_to_fulltext, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        print_status_block(
            DOWNLOAD_PROGRESS_LOG_FILE,
            "DOI to publication full text mapping written",
            f"DOIs attempted: {len(doi_order)}",
            f"Full texts retrieved: {sum(path is not None for path in results_by_doi.values())}",
            f"Output file: {doi_to_fulltext_output_file}",
        )

    return [
        full_text_path
        for doi in doi_order
        if (full_text_path := results_by_doi[doi]) is not None
    ]


def bulk_convert_full_texts_to_md(
    input_files: list[str | Path],
    output_dir: Path | str = FULL_TEXT_MD_DIR,
    remove_references: bool = True,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
) -> list[Path]:
    """Bulk convert publication full text files (PDF or XML) to Markdown."""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")

    files = [Path(file_path).resolve() for file_path in input_files]
    print_status_block(
        DOWNLOAD_PROGRESS_LOG_FILE,
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
                append_log_line(
                    DOWNLOAD_PROGRESS_LOG_FILE,
                    f"Markdown conversion progress: {completed_files}/{len(files)} files processed; latest file: {file_path}; status: {'ok' if md_file else 'no-output'}",
                )
            except Exception as exc:
                completed_files += 1
                append_log_line(
                    DOWNLOAD_PROGRESS_LOG_FILE,
                    f"Markdown conversion progress: {completed_files}/{len(files)} files processed; latest file: {file_path}; status: error; error: {exc}",
                )
                print_status_block(
                    DOWNLOAD_PROGRESS_LOG_FILE,
                    "Error converting publication full text to Markdown",
                    f"File: {file_path}",
                    f"Error: {exc}",
                )

    converted_files = [
        converted_files_by_index[index] for index in sorted(converted_files_by_index)
    ]
    print_status_block(
        DOWNLOAD_PROGRESS_LOG_FILE,
        "Publication full text Markdown conversion complete",
        f"Files queued: {len(files)}",
        f"Files converted: {len(converted_files)}",
        f"Files failed: {len(files) - len(converted_files)}",
        f"Output directory: {output_dir}",
        f"Progress log: {DOWNLOAD_PROGRESS_LOG_FILE}",
    )
    return converted_files
