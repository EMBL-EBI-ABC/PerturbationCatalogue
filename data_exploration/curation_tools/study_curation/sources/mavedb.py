"""Pipeline for MaveDB publication text collection and prompt context preparation."""

import argparse
import json
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from time import sleep

import requests

from curation_tools.study_curation.logging_utils import (
    append_log_line,
    print_status_block,
)
from curation_tools.study_curation.paths import (
    MAVEDB_DOI_TO_FULLTEXT_OUTPUT_FILE,
    MAVEDB_DUMP_DIR,
    MAVEDB_METADATA_OUTPUT_DIR,
    MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
)
from curation_tools.study_curation.sources.publication_text import (
    DEFAULT_DOWNLOAD_MAX_WORKERS,
    DOWNLOAD_PROGRESS_LOG_FILE,
    FULL_TEXT_MD_DIR,
    PAPERSCRAPER_FULL_TEXT_RAW_DIR,
    bulk_convert_full_texts_to_md,
    bulk_download_pub_full_texts,
)

MAVEDB_API_BASE_URL = "https://api.mavedb.org/api/v1"
DEFAULT_FETCH_SLEEP_TIME = 0.1
JSON_INDENT = 2

DEFAULT_EXCLUDED_DOIS: tuple[str, ...] = (
    "10.1186/s13059-017-1272-5",  # Enrich2 software paper
    "10.1038/nmeth.1492",  # Fowler 2010 DMS method paper
    "10.1038/s41588-018-0122-z",  # VAMP-seq method paper
    "10.1101/2024.04.26.591310",  # domainome paper
)


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


def filter_excluded_dois(
    dois: list[str] | tuple[str, ...] | None,
    excluded_dois: (
        set[str] | list[str] | tuple[str, ...] | str | Path | None
    ) = DEFAULT_EXCLUDED_DOIS,
) -> list[str]:
    """Remove configured excluded DOIs while preserving input order."""
    parsed_excluded = parse_excluded_dois(excluded_dois)
    return [doi for doi in dois or [] if not _is_doi_excluded(doi, parsed_excluded)]


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
        "perturbed_target_symbol": target_genes or None,
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
    """Parse command-line options and run the publication text pipeline."""
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
