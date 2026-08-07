"""Pipeline for MaveDB publication text collection and prompt context preparation."""

import argparse
import json
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from time import sleep
from typing import Type

import instructor
import requests
from pydantic import BaseModel
from tqdm import tqdm

from curation_tools.llm_curation.logging_utils import (
    append_log_line,
    print_status_block,
)
from curation_tools.llm_curation.metadata_extraction import (
    DEFAULT_LLM_MODEL_NAME,
    _ensure_log_file,
    create_csv_from_curated_metadata_json,
    dataset_id_from_source_urns,
    render_metadata_extraction_prompt,
)
from curation_tools.llm_curation.publication_text import (
    DEFAULT_DOWNLOAD_MAX_WORKERS,
    DOWNLOAD_PROGRESS_LOG_FILE,
    FULL_TEXT_MD_DIR,
    PAPERSCRAPER_FULL_TEXT_RAW_DIR,
    bulk_convert_full_texts_to_md,
    bulk_download_pub_full_texts,
)

REPO_ROOT = Path(__file__).resolve().parents[4]
MAVEDB_DIR = REPO_ROOT / "data_exploration" / "MaveDB"
LLM_CURATION_DIR = REPO_ROOT / "data_exploration" / "curation_tools" / "llm_curation"

MAVEDB_DUMP_DIR = MAVEDB_DIR / "Dump" / "mavedb-dump.20250612164404" / "csv"
MAVEDB_METADATA_OUTPUT_DIR = MAVEDB_DIR / "llm_metadata_extraction" / "mavedb_metadata"
MAVEDB_URN_TO_DOIS_OUTPUT_FILE = (
    MAVEDB_DIR / "llm_metadata_extraction" / "mavedb_urn_to_dois.json"
)
MAVEDB_DOI_TO_FULLTEXT_OUTPUT_FILE = (
    MAVEDB_DIR / "llm_metadata_extraction" / "mavedb_doi_to_fulltext.json"
)
MAVEDB_API_BASE_URL = "https://api.mavedb.org/api/v1"
DEFAULT_FETCH_SLEEP_TIME = 0.1
JSON_INDENT = 2

DEFAULT_EXCLUDED_DOIS: tuple[str, ...] = (
    "10.1186/s13059-017-1272-5",  # Enrich2 software paper
    "10.1038/nmeth.1492",  # Fowler 2010 DMS method paper
    "10.1038/s41588-018-0122-z",  # VAMP-seq method paper
)

_MAVEDB_URN_TO_DOIS_CACHE: dict[str, list[str]] | None = None


def parse_excluded_dois(
    value: str | list[str] | tuple[str, ...] | set[str] | Path | None,
) -> tuple[str, ...]:
    """Parse user-configured excluded DOIs from a string, list, set, or file path."""
    if value is None:
        return DEFAULT_EXCLUDED_DOIS
    if isinstance(value, (list, tuple, set)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    if isinstance(value, Path) or (
        isinstance(value, str)
        and (Path(value).is_file() or "/" in value or "\\" in value)
    ):
        path = Path(value)
        if path.is_file():
            lines = path.read_text(encoding="utf-8").splitlines()
            parsed = [
                line.strip()
                for line in lines
                if line.strip() and not line.strip().startswith("#")
            ]
            return tuple(parsed)
    if isinstance(value, str):
        items = [item.strip() for item in re.split(r"[,; \s]+", value) if item.strip()]
        return tuple(items)
    return DEFAULT_EXCLUDED_DOIS


def _is_doi_excluded(
    doi: str,
    excluded_dois: set[str] | list[str] | tuple[str, ...] | None,
) -> bool:
    """Check if a DOI matches any entry in the user-configured exclusion list."""
    if not excluded_dois:
        return False
    doi_norm = doi.strip().lower()
    doi_lookup_norm = format_identifier_for_lookup(doi).lower()
    doi_alnum = re.sub(r"[^a-z0-9]", "", doi_norm)
    for exc in excluded_dois:
        exc_norm = str(exc).strip().lower()
        exc_lookup_norm = format_identifier_for_lookup(str(exc)).lower()
        exc_alnum = re.sub(r"[^a-z0-9]", "", exc_norm)
        if (
            doi_norm == exc_norm
            or doi_lookup_norm == exc_lookup_norm
            or doi_alnum == exc_alnum
        ):
            return True
    return False


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for the MaveDB text collection pipeline."""
    parser = argparse.ArgumentParser(
        description="Collect MaveDB publication full text and convert it to Markdown.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--dump-dir",
        type=Path,
        default=MAVEDB_DUMP_DIR,
        help="Directory containing the MaveDB CSV dump.",
    )
    parser.add_argument(
        "--metadata-output-dir",
        type=Path,
        default=MAVEDB_METADATA_OUTPUT_DIR,
        help="Directory used to cache MaveDB entry metadata JSON files.",
    )
    parser.add_argument(
        "--urn-to-dois-output-file",
        type=Path,
        default=MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
        help="Output file for the URN-to-DOI mapping JSON.",
    )
    parser.add_argument(
        "--doi-to-fulltext-output-file",
        type=Path,
        default=MAVEDB_DOI_TO_FULLTEXT_OUTPUT_FILE,
        help="Output file for the DOI-to-full-text mapping JSON.",
    )
    parser.add_argument(
        "--raw-output-dir",
        type=Path,
        default=PAPERSCRAPER_FULL_TEXT_RAW_DIR,
        help="Directory used to cache downloaded publication full text files.",
    )
    parser.add_argument(
        "--markdown-output-dir",
        type=Path,
        default=FULL_TEXT_MD_DIR,
        help="Directory used to write converted Markdown files.",
    )
    parser.add_argument(
        "--include-secondary",
        action="store_true",
        help="Include secondary publication DOIs even when primary publication DOIs exist.",
    )
    parser.add_argument(
        "--excluded-dois",
        type=str,
        default=None,
        help=(
            "User-configurable list of DOIs to exclude (e.g. software/method papers). "
            "Can be a comma-separated string, space-separated string, or path to a text file containing DOIs."
        ),
    )
    return parser


def get_unique_mavedb_urns(dump_dir: str | Path = MAVEDB_DUMP_DIR) -> list[str]:
    """Extract unique file basenames from CSV files in the dump directory."""
    dump_dir = Path(dump_dir).resolve()
    unique_files = set()
    for file_path in dump_dir.iterdir():
        if file_path.suffix == ".csv":
            cleaned = file_path.name.replace(".scores.csv", "").replace(
                ".counts.csv", ""
            )
            cleaned = cleaned.replace("urn-mavedb-", "urn:mavedb:")
            unique_files.add(cleaned)
    return list(unique_files)


def get_mavedb_record_type_from_urn(urn: str) -> str:
    """Infer the MaveDB record type from the URN structure."""
    identifier = urn.replace("urn:mavedb:", "")
    segment_count = len(identifier.split("-"))
    if segment_count >= 3:
        return "score-sets"
    if segment_count == 2:
        return "experiments"
    if segment_count == 1:
        return "experiment-sets"
    raise ValueError(f"Unable to infer MaveDB record type from URN: {urn}")


def fetch_mavedb_entry(urn: str) -> dict:
    """Fetch entry data from the MaveDB API for a given URN."""
    record_type = get_mavedb_record_type_from_urn(urn)
    response = requests.get(f"{MAVEDB_API_BASE_URL}/{record_type}/{urn}")
    response.raise_for_status()
    return response.json()


def get_dois_from_mavedb_entry(
    entry: dict,
    primary_only: bool = True,
    excluded_dois: (
        set[str] | list[str] | tuple[str, ...] | str | Path | None
    ) = DEFAULT_EXCLUDED_DOIS,
    log: bool = True,
) -> list[str] | None:
    """Extract DOI values from a MaveDB entry, prioritizing primary publications."""
    experiment = entry.get("experiment") or {}
    parsed_excluded = parse_excluded_dois(excluded_dois)

    primary_sources = (
        (entry.get("primaryPublicationIdentifiers") or [], "doi"),
        (entry.get("doiIdentifiers") or [], "identifier"),
        (experiment.get("primaryPublicationIdentifiers") or [], "doi"),
    )
    secondary_sources = (
        (entry.get("secondaryPublicationIdentifiers") or [], "doi"),
        (experiment.get("secondaryPublicationIdentifiers") or [], "doi"),
    )

    raw_primary_dois = sorted(
        {
            doi
            for identifiers, field_name in primary_sources
            for identifier in identifiers
            if isinstance(identifier, dict)
            if (doi := identifier.get(field_name))
        }
    )
    raw_secondary_dois = sorted(
        {
            doi
            for identifiers, field_name in secondary_sources
            for identifier in identifiers
            if isinstance(identifier, dict)
            if (doi := identifier.get(field_name))
        }
    )

    primary_dois = [
        doi for doi in raw_primary_dois if not _is_doi_excluded(doi, parsed_excluded)
    ]
    secondary_dois = [
        doi for doi in raw_secondary_dois if not _is_doi_excluded(doi, parsed_excluded)
    ]

    if primary_dois:
        dois = primary_dois
    elif secondary_dois and not primary_only:
        dois = secondary_dois
    elif secondary_dois:
        dois = secondary_dois
    elif raw_primary_dois:
        # Fallback: if all primary DOIs were excluded and no other non-excluded DOIs exist, retain raw primary DOIs
        dois = raw_primary_dois
    elif raw_secondary_dois:
        # Fallback: if all secondary DOIs were excluded and no other non-excluded DOIs exist, retain raw secondary DOIs
        dois = raw_secondary_dois
    else:
        dois = None

    if dois:
        if log:
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "MaveDB publication identifiers found",
                f"Count: {len(dois)}",
                f"DOIs: {', '.join(dois)}",
            )
        return dois
    if log:
        print_status_block(
            DOWNLOAD_PROGRESS_LOG_FILE,
            "No publication identifiers found",
            f"Entry: {entry.get('urn', '<unknown URN>')}",
        )
    return None


def export_mavedb_urn_to_dois_json(
    mavedb_entries_dir: Path | str = MAVEDB_METADATA_OUTPUT_DIR,
    output_file: Path | str = MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
    primary_only: bool = True,
    excluded_dois: (
        set[str] | list[str] | tuple[str, ...] | str | Path | None
    ) = DEFAULT_EXCLUDED_DOIS,
) -> dict[str, list[str]]:
    """Write a JSON mapping from MaveDB URNs to publication DOIs."""
    mavedb_entries_dir = Path(mavedb_entries_dir).resolve()
    output_file = Path(output_file).resolve()

    if not mavedb_entries_dir.is_dir():
        raise ValueError(f"MaveDB entries directory not found: {mavedb_entries_dir}")

    urn_to_dois: dict[str, list[str]] = {}
    entry_files = sorted(mavedb_entries_dir.glob("*.json"))
    for entry_file in tqdm(
        entry_files, desc="Building URN to DOI mapping", unit="entry"
    ):
        try:
            entry = json.loads(entry_file.read_text(encoding="utf-8"))
            urn = entry.get("urn")
            if not urn:
                continue
            dois = get_dois_from_mavedb_entry(
                entry,
                primary_only=primary_only,
                excluded_dois=excluded_dois,
                log=False,
            )
            if dois:
                urn_to_dois[urn] = dois
        except Exception as exc:
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "Error exporting URN to DOI mapping",
                f"File: {entry_file}",
                f"Error: {exc}",
            )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(
        json.dumps(urn_to_dois, indent=2, sort_keys=True), encoding="utf-8"
    )

    print_status_block(
        DOWNLOAD_PROGRESS_LOG_FILE,
        "MaveDB URN to DOI mapping written",
        f"Entries scanned: {len(entry_files)}",
        f"URNs with DOIs: {len(urn_to_dois)}",
        f"Output file: {output_file}",
    )
    return urn_to_dois


def bulk_fetch_mavedb_entries(
    input_files: list[str],
    output_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
    overwrite: bool = False,
    sleep_time: float = DEFAULT_FETCH_SLEEP_TIME,
) -> list[dict]:
    """Fetch MaveDB entries for a list of URNs and cache them as JSON."""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    for urn in tqdm(input_files, desc="Fetching MaveDB entries", unit="entry"):
        entry_output_path = output_dir / f"{urn.replace(':', '_')}.json"
        try:
            if entry_output_path.is_file() and not overwrite:
                print_status_block(
                    DOWNLOAD_PROGRESS_LOG_FILE,
                    "MaveDB entry already present",
                    f"URN: {urn}",
                    f"Skipping download for file: {entry_output_path}",
                )
                entries.append(
                    json.loads(entry_output_path.read_text(encoding="utf-8"))
                )
                continue

            entry = fetch_mavedb_entry(urn)
            entries.append(entry)
            sleep(sleep_time)
            entry_output_path.write_text(json.dumps(entry, indent=2), encoding="utf-8")
        except Exception as exc:
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "Error fetching MaveDB entry",
                f"URN: {urn}",
                f"Error: {exc}",
            )
    return entries


def collect_publication_dois(
    mavedb_entries_dir: Path | str = MAVEDB_METADATA_OUTPUT_DIR,
    primary_only: bool = True,
    excluded_dois: (
        set[str] | list[str] | tuple[str, ...] | str | Path | None
    ) = DEFAULT_EXCLUDED_DOIS,
) -> dict[str, set[str]]:
    """Collect a DOI-to-URN mapping from cached MaveDB entry JSON files."""
    mavedb_entries_dir = Path(mavedb_entries_dir).resolve()
    if not mavedb_entries_dir.is_dir():
        raise ValueError(f"MaveDB entries directory not found: {mavedb_entries_dir}")

    doi_to_urns: dict[str, set[str]] = {}
    doi_reference_count = 0
    entries_with_dois = 0
    entry_files = list(mavedb_entries_dir.glob("*.json"))
    for entry_file in tqdm(
        entry_files, desc="Collecting publication DOIs", unit="entry"
    ):
        try:
            entry = json.loads(entry_file.read_text(encoding="utf-8"))
            urn = entry.get("urn", "<unknown URN>")
            dois = get_dois_from_mavedb_entry(
                entry,
                primary_only=primary_only,
                excluded_dois=excluded_dois,
                log=False,
            )
            if dois:
                entries_with_dois += 1
                doi_reference_count += len(dois)
                for doi in dois:
                    doi_to_urns.setdefault(doi, set()).add(urn)
        except Exception as exc:
            print_status_block(
                DOWNLOAD_PROGRESS_LOG_FILE,
                "Error processing MaveDB entry file",
                f"File: {entry_file}",
                f"Error: {exc}",
            )

    print_status_block(
        DOWNLOAD_PROGRESS_LOG_FILE,
        "Publication DOI collection complete",
        f"Entries scanned: {len(entry_files)}",
        f"Entries with DOIs: {entries_with_dois}",
        f"DOI references found: {doi_reference_count}",
        f"Unique DOIs found: {len(doi_to_urns)}",
        f"Progress log: {DOWNLOAD_PROGRESS_LOG_FILE}",
    )
    return doi_to_urns


def normalize_prompt_text(value: str | None) -> str | None:
    """Normalize prompt text by collapsing whitespace and dropping empty values."""
    if not value:
        return None
    normalized_value = " ".join(str(value).split())
    return normalized_value or None


def format_identifier_for_lookup(identifier: str) -> str:
    """Normalize a publication identifier into the filename-safe lookup form."""
    return str(identifier).replace(".", "_").replace("/", "_")


def format_urn_for_filename(urn: str) -> str:
    """Convert a MaveDB URN into the filename-safe form used on disk."""
    return str(urn).replace(":", "_")


def load_mavedb_urn_to_dois(
    mapping_file: str | Path = MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
) -> dict[str, list[str]]:
    """Load and cache the mapping from MaveDB URNs to publication identifiers."""
    global _MAVEDB_URN_TO_DOIS_CACHE
    mapping_file = Path(mapping_file).resolve()
    if _MAVEDB_URN_TO_DOIS_CACHE is None:
        _MAVEDB_URN_TO_DOIS_CACHE = json.loads(mapping_file.read_text(encoding="utf-8"))
    return _MAVEDB_URN_TO_DOIS_CACHE


def find_matching_mavedb_entry_paths(
    publication_full_text_path: str | Path,
    metadata_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
    mapping_file: str | Path = MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
) -> list[tuple[str, Path]]:
    """Return MaveDB metadata files whose identifiers match a publication file."""
    publication_full_text_path = Path(publication_full_text_path).resolve()
    publication_identifier = publication_full_text_path.stem
    metadata_dir = Path(metadata_dir).resolve()
    urn_to_dois = load_mavedb_urn_to_dois(mapping_file)
    matching_entries: list[tuple[str, Path]] = []
    for urn, identifiers in sorted(urn_to_dois.items()):
        if any(
            format_identifier_for_lookup(identifier) == publication_identifier
            for identifier in identifiers
        ):
            entry_path = metadata_dir / f"{format_urn_for_filename(urn)}.json"
            if entry_path.is_file():
                matching_entries.append((urn, entry_path))
    return matching_entries


def extract_curated_mavedb_prompt_metadata(entry_payload: dict) -> dict[str, object]:
    """Extract the subset of MaveDB entry fields used to condition the prompt."""
    experiment_payload = entry_payload.get("experiment") or {}
    target_genes = sorted(
        {
            gene.get("name")
            for gene in entry_payload.get("targetGenes", [])
            if gene.get("name")
        }
    )
    primary_publications = entry_payload.get("primaryPublicationIdentifiers", [])
    publication_titles = sorted(
        {
            publication.get("title")
            for publication in primary_publications
            if publication.get("title")
        }
    )
    publication_years = sorted(
        {
            publication.get("publicationYear")
            for publication in primary_publications
            if publication.get("publicationYear")
        }
    )
    publication_dois = sorted(
        {
            publication.get("doi")
            for publication in primary_publications
            if publication.get("doi")
        }
    )
    license_info = entry_payload.get("license") or {}
    license_short_name = license_info.get("shortName")

    curated_metadata: dict[str, object] = {
        "score_set_title": normalize_prompt_text(entry_payload.get("title")),
        "score_set_short_description": normalize_prompt_text(
            entry_payload.get("shortDescription")
        ),
        "score_set_abstract": normalize_prompt_text(entry_payload.get("abstractText")),
        "score_set_method": normalize_prompt_text(entry_payload.get("methodText")),
        "experiment_title": normalize_prompt_text(experiment_payload.get("title")),
        "experiment_short_description": normalize_prompt_text(
            experiment_payload.get("shortDescription")
        ),
        "experiment_abstract": normalize_prompt_text(
            experiment_payload.get("abstractText")
        ),
        "experiment_method": normalize_prompt_text(
            experiment_payload.get("methodText")
        ),
        "target_genes": target_genes or None,
        "score_columns": entry_payload.get("datasetColumns", {}).get("scoreColumns")
        or None,
        "primary_publication_titles": publication_titles or None,
        "primary_publication_years": publication_years or None,
        "total_variants": entry_payload.get("numVariants"),
        "primary_publication_dois": publication_dois or None,
        "license": license_short_name,
    }
    return {
        field_name: field_value
        for field_name, field_value in curated_metadata.items()
        if field_value not in (None, [], {})
    }


def build_mavedb_context_signature(curated_metadata: dict[str, object]) -> str:
    """Build a stable signature for grouping equivalent prompt contexts."""
    signature_payload = {
        "score_set_short_description": curated_metadata.get(
            "score_set_short_description"
        ),
        "score_set_method": curated_metadata.get("score_set_method"),
        "experiment_short_description": curated_metadata.get(
            "experiment_short_description"
        ),
        "experiment_method": curated_metadata.get("experiment_method"),
    }
    if not any(signature_payload.values()):
        signature_payload = {
            "score_set_title": curated_metadata.get("score_set_title"),
            "experiment_title": curated_metadata.get("experiment_title"),
        }
    return json.dumps(signature_payload, sort_keys=True, ensure_ascii=True)


def merge_prompt_metadata_value(existing_value: object, new_value: object) -> object:
    """Merge prompt metadata values while preserving scalar and list semantics."""
    if existing_value == new_value or new_value in (None, [], {}):
        return existing_value
    if existing_value in (None, [], {}):
        return new_value

    if isinstance(existing_value, list):
        merged_values = list(existing_value)
    else:
        merged_values = [existing_value]

    if isinstance(new_value, list):
        for value in new_value:
            if value not in merged_values:
                merged_values.append(value)
    elif new_value not in merged_values:
        merged_values.append(new_value)
    return merged_values


def build_mavedb_prompt_contexts(
    publication_full_text_path: str | Path,
    metadata_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
    mapping_file: str | Path = MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
) -> list[dict[str, object]]:
    """Build deduplicated prompt contexts for a publication."""
    matching_entries = find_matching_mavedb_entry_paths(
        publication_full_text_path=publication_full_text_path,
        metadata_dir=metadata_dir,
        mapping_file=mapping_file,
    )
    if not matching_entries:
        return []

    grouped_contexts: dict[str, dict[str, object]] = {}
    for urn, entry_path in matching_entries:
        entry_payload = json.loads(entry_path.read_text(encoding="utf-8"))
        curated_metadata = extract_curated_mavedb_prompt_metadata(entry_payload)
        context_signature = build_mavedb_context_signature(curated_metadata)
        if context_signature not in grouped_contexts:
            grouped_contexts[context_signature] = {
                "source_urns": [urn],
                "source_files": [entry_path.name],
                "metadata": curated_metadata,
            }
            continue

        grouped_contexts[context_signature]["source_urns"].append(urn)
        grouped_contexts[context_signature]["source_files"].append(entry_path.name)
        merged_metadata = grouped_contexts[context_signature]["metadata"]
        for field_name, field_value in curated_metadata.items():
            merged_metadata[field_name] = merge_prompt_metadata_value(
                merged_metadata.get(field_name), field_value
            )

    return sorted(
        grouped_contexts.values(),
        key=lambda context: tuple(context["source_urns"]),
    )


def prompt_context_builder(publication_full_text_path: Path) -> list[dict[str, object]]:
    """Build supplementary MaveDB prompt contexts for a publication text file."""
    return build_mavedb_prompt_contexts(
        publication_full_text_path,
        metadata_dir=MAVEDB_METADATA_OUTPUT_DIR,
        mapping_file=MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
    )


def context_output_suffix_builder(
    prompt_context: dict[str, object] | None,
    context_index: int,
    total_contexts: int,
) -> str:
    """Build a stable output filename suffix for a MaveDB prompt context."""
    base_suffix = "" if total_contexts == 1 else f"__ctx_{context_index:02d}"
    if not prompt_context or not prompt_context.get("source_urns"):
        return base_suffix

    source_urns = prompt_context.get("source_urns", [])
    urn_suffix = "__" + "__".join(format_urn_for_filename(urn) for urn in source_urns)
    return f"{base_suffix}{urn_suffix}" if base_suffix else urn_suffix


def output_metadata_builder(
    prompt_context: dict[str, object] | None,
) -> dict[str, object]:
    """Preserve MaveDB source provenance in each extraction output payload."""
    if not prompt_context:
        return {}
    return {
        "__source_urns": list(prompt_context.get("source_urns", [])),
        "__source_files": list(prompt_context.get("source_files", [])),
    }


def format_supplementary_mavedb_metadata(
    prompt_context: dict[str, object] | None,
) -> str:
    """Render supplementary MaveDB prompt context as JSON or a fallback message."""
    if not prompt_context:
        return "No supplementary MaveDB metadata was available for this publication."

    return json.dumps(
        {
            "source_urns": prompt_context["source_urns"],
            "source_files": prompt_context["source_files"],
            "curated_mavedb_metadata": prompt_context["metadata"],
        },
        indent=JSON_INDENT,
        ensure_ascii=True,
    )


def build_mavedb_publication_full_text(
    publication_full_text_path: Path,
    prompt_context: dict[str, object] | None = None,
    publication_full_text_dir: Path | str = FULL_TEXT_MD_DIR,
    mapping_file: Path | str = MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
) -> str:
    """Build publication full text for a prompt, concatenating multiple primary papers if linked to the same URN."""
    publication_full_text_path = Path(publication_full_text_path).resolve()
    publication_full_text_dir = Path(publication_full_text_dir).resolve()

    if not prompt_context or not prompt_context.get("source_urns"):
        return publication_full_text_path.read_text(encoding="utf-8")

    urn_to_dois = load_mavedb_urn_to_dois(mapping_file)
    source_urns = prompt_context.get("source_urns", [])

    all_dois_for_urns: list[str] = []
    for urn in source_urns:
        for doi in urn_to_dois.get(urn, []):
            if doi not in all_dois_for_urns:
                all_dois_for_urns.append(doi)

    matching_md_files: list[tuple[str, Path]] = []
    for doi in all_dois_for_urns:
        stem = format_identifier_for_lookup(doi)
        md_path = publication_full_text_dir / f"{stem}.md"
        if md_path.is_file():
            matching_md_files.append((doi, md_path))

    if not matching_md_files:
        return publication_full_text_path.read_text(encoding="utf-8")

    if len(matching_md_files) == 1:
        return matching_md_files[0][1].read_text(encoding="utf-8")

    sections = []
    for idx, (doi, md_path) in enumerate(matching_md_files, start=1):
        content = md_path.read_text(encoding="utf-8").strip()
        sections.append(f"# Primary Publication {idx} (DOI: {doi})\n\n{content}")

    return "\n\n---\n\n".join(sections)


def extract_evidence_for_mavedb_urn(
    urn: str,
    dois: list[str],
    extraction_schema: Type[BaseModel],
    output_dir: str | Path,
    log_file: str | Path,
    prompt_template_file: str | Path,
    metadata_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
    publication_full_text_dir: str | Path = FULL_TEXT_MD_DIR,
    overwrite: bool = False,
    model_name: str = DEFAULT_LLM_MODEL_NAME,
    prompt_context_formatter=format_supplementary_mavedb_metadata,
    output_metadata_builder=output_metadata_builder,
    verbose: bool = False,
) -> Path | None:
    """Extract evidence for a single MaveDB URN, loading its metadata and concatenating its primary papers."""
    output_dir = Path(output_dir).resolve()
    metadata_dir = Path(metadata_dir).resolve()
    publication_full_text_dir = Path(publication_full_text_dir).resolve()
    prompt_template_file = Path(prompt_template_file).resolve()
    log_file = _ensure_log_file(log_file)

    urn_filename_stem = format_urn_for_filename(urn)
    output_path = output_dir / f"{urn_filename_stem}.json"

    if not overwrite and output_path.is_file():
        print_status_block(
            log_file,
            "MaveDB URN evidence extraction skipped - output already exists",
            f"URN: {urn}",
            f"Evidence output: {output_path}",
        )
        return output_path

    entry_metadata_path = metadata_dir / f"{urn_filename_stem}.json"
    if entry_metadata_path.is_file():
        entry_payload = json.loads(entry_metadata_path.read_text(encoding="utf-8"))
        curated_metadata = extract_curated_mavedb_prompt_metadata(entry_payload)
        source_files = [entry_metadata_path.name]
    else:
        curated_metadata = {}
        source_files = []

    prompt_context = {
        "source_urns": [urn],
        "source_files": source_files,
        "metadata": curated_metadata,
    }

    matching_md_files: list[tuple[str, Path]] = []
    for doi in dois:
        stem = format_identifier_for_lookup(doi)
        md_path = publication_full_text_dir / f"{stem}.md"
        if md_path.is_file():
            matching_md_files.append((doi, md_path))

    if matching_md_files:
        if len(matching_md_files) == 1:
            publication_full_text = matching_md_files[0][1].read_text(encoding="utf-8")
        else:
            sections = []
            for idx, (doi, md_path) in enumerate(matching_md_files, start=1):
                content = md_path.read_text(encoding="utf-8").strip()
                sections.append(
                    f"# Primary Publication {idx} (DOI: {doi})\n\n{content}"
                )
            publication_full_text = "\n\n---\n\n".join(sections)
    else:
        publication_full_text = (
            "No publication full text available for this MaveDB dataset."
        )

    prompt_template = prompt_template_file.read_text(encoding="utf-8")
    supplementary_metadata = prompt_context_formatter(prompt_context)
    prompt = render_metadata_extraction_prompt(
        prompt_template=prompt_template,
        publication_full_text=publication_full_text,
        supplementary_metadata=supplementary_metadata,
    )

    if verbose:
        print_status_block(
            log_file,
            "[VERBOSE] MaveDB URN full prompt, metadata, and publication text",
            f"URN: {urn}",
            f"DOIs: {dois}",
            "----- PROMPT START -----",
            prompt,
            "----- PROMPT END -----",
        )

    client = instructor.from_provider(
        model_name,
        location="global",
        vertexai=True,
    )
    extraction_response = client.create(
        response_model=extraction_schema,
        messages=[{"role": "user", "content": prompt}],
        thinking_config={"thinking_level": "high"},
        generation_config={"temperature": 0.2},
    )

    output_metadata = output_metadata_builder(prompt_context)
    evidence_payload = extraction_response.model_dump()
    evidence_payload.update(output_metadata)
    evidence_payload["curation_agent_type"] = "LLM"
    evidence_payload["curation_agent_name"] = model_name

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(evidence_payload, indent=JSON_INDENT),
        encoding="utf-8",
    )

    print_status_block(
        log_file,
        "MaveDB URN evidence extraction complete",
        f"URN: {urn}",
        f"Primary DOIs: {dois}",
        f"Output: {output_path}",
    )
    return output_path


def bulk_extract_evidence_for_mavedb_urns(
    urn_to_dois: dict[str, list[str]],
    extraction_schema: Type[BaseModel],
    output_dir: str | Path,
    log_file: str | Path,
    prompt_template_file: str | Path,
    metadata_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
    publication_full_text_dir: str | Path = FULL_TEXT_MD_DIR,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
    overwrite: bool = False,
    model_name: str = DEFAULT_LLM_MODEL_NAME,
    create_csv: bool = True,
    verbose: bool = False,
) -> list[Path]:
    """Bulk extract evidence for MaveDB URNs in parallel."""
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    prompt_template_file = Path(prompt_template_file).resolve()

    urn_list = sorted(urn_to_dois.keys())
    print_status_block(
        log_file,
        "Starting bulk URN-centric evidence extraction",
        f"MaveDB URNs queued: {len(urn_list)}",
        f"Output directory: {output_dir}",
        f"Model: {model_name}",
        f"Max workers: {max_workers}",
        f"Overwrite: {overwrite}",
        f"Log file: {log_file}",
    )

    extracted_outputs: list[Path] = []
    completed_urns = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_urn = {
            executor.submit(
                extract_evidence_for_mavedb_urn,
                urn=urn,
                dois=urn_to_dois[urn],
                extraction_schema=extraction_schema,
                output_dir=output_dir,
                log_file=log_file,
                prompt_template_file=prompt_template_file,
                metadata_dir=metadata_dir,
                publication_full_text_dir=publication_full_text_dir,
                overwrite=overwrite,
                model_name=model_name,
                verbose=verbose,
            ): urn
            for urn in urn_list
        }

        for future in tqdm(
            as_completed(future_to_urn),
            total=len(future_to_urn),
            desc="Extracting MaveDB URN evidence",
            unit="urn",
        ):
            urn = future_to_urn[future]
            try:
                out_path = future.result()
                if out_path:
                    extracted_outputs.append(out_path)
                completed_urns += 1
                append_log_line(
                    log_file,
                    f"Bulk URN evidence extraction progress: {completed_urns}/{len(urn_list)} processed; URN: {urn}; status: ok",
                )
            except Exception as exc:
                completed_urns += 1
                append_log_line(
                    log_file,
                    f"Bulk URN evidence extraction progress: {completed_urns}/{len(urn_list)} processed; URN: {urn}; status: error; error: {exc}",
                )
                print_status_block(
                    log_file,
                    "Bulk URN evidence extraction failed",
                    f"URN: {urn}",
                    f"Error: {exc}",
                )

    if create_csv and extracted_outputs:
        create_csv_from_curated_metadata_json(
            input_dir=output_dir,
            output_csv_path=output_dir / "clean_metadata.csv",
            log_file=log_file,
        )

    return extracted_outputs


def run_full_text_collection_pipeline(
    dump_dir: str | Path = MAVEDB_DUMP_DIR,
    metadata_output_dir: str | Path = MAVEDB_METADATA_OUTPUT_DIR,
    urn_to_dois_output_file: str | Path = MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
    doi_to_fulltext_output_file: str | Path = MAVEDB_DOI_TO_FULLTEXT_OUTPUT_FILE,
    raw_output_dir: str | Path = PAPERSCRAPER_FULL_TEXT_RAW_DIR,
    markdown_output_dir: str | Path = FULL_TEXT_MD_DIR,
    overwrite: bool = False,
    max_workers: int = DEFAULT_DOWNLOAD_MAX_WORKERS,
    primary_only: bool = True,
    excluded_dois: (
        set[str] | list[str] | tuple[str, ...] | str | Path | None
    ) = DEFAULT_EXCLUDED_DOIS,
) -> None:
    """Run the end-to-end MaveDB publication text collection pipeline."""
    unique_files = get_unique_mavedb_urns(dump_dir)
    bulk_fetch_mavedb_entries(
        unique_files,
        output_dir=metadata_output_dir,
        overwrite=overwrite,
    )
    export_mavedb_urn_to_dois_json(
        mavedb_entries_dir=metadata_output_dir,
        output_file=urn_to_dois_output_file,
        primary_only=primary_only,
        excluded_dois=excluded_dois,
    )
    doi_to_urns = collect_publication_dois(
        mavedb_entries_dir=metadata_output_dir,
        primary_only=primary_only,
        excluded_dois=excluded_dois,
    )
    full_text_paths = bulk_download_pub_full_texts(
        doi_to_urns=doi_to_urns,
        output_dir=raw_output_dir,
        doi_to_fulltext_output_file=doi_to_fulltext_output_file,
        overwrite=overwrite,
        max_workers=max_workers,
    )
    bulk_convert_full_texts_to_md(
        full_text_paths,
        output_dir=markdown_output_dir,
        remove_references=True,
        max_workers=max_workers,
    )


def main() -> None:
    args = build_parser().parse_args()
    excluded_dois = (
        parse_excluded_dois(args.excluded_dois)
        if args.excluded_dois
        else DEFAULT_EXCLUDED_DOIS
    )
    run_full_text_collection_pipeline(
        dump_dir=args.dump_dir,
        metadata_output_dir=args.metadata_output_dir,
        urn_to_dois_output_file=args.urn_to_dois_output_file,
        doi_to_fulltext_output_file=args.doi_to_fulltext_output_file,
        raw_output_dir=args.raw_output_dir,
        markdown_output_dir=args.markdown_output_dir,
        primary_only=not args.include_secondary,
        excluded_dois=excluded_dois,
    )


if __name__ == "__main__":
    main()
