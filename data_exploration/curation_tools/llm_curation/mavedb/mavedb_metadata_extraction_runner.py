"""CLI runner for extracting structured metadata from MaveDB-linked papers."""

import argparse
import json
import re
from pathlib import Path

from curation_tools.llm_curation.metadata_extraction import (
    DEFAULT_DOWNLOAD_MAX_WORKERS,
    DEFAULT_LLM_MODEL_NAME,
    bulk_extract_evidence_from_publications,
)
from curation_tools.llm_curation.schema_loading import load_extraction_schema
from curation_tools.llm_curation.mavedb.processing import (
    MAVEDB_METADATA_OUTPUT_DIR,
    MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
    FULL_TEXT_MD_DIR,
    build_mavedb_publication_full_text,
    bulk_extract_evidence_for_mavedb_urns,
    context_output_suffix_builder,
    format_supplementary_mavedb_metadata,
    format_urn_for_filename,
    get_dois_from_mavedb_entry,
    load_mavedb_urn_to_dois,
    output_metadata_builder,
    prompt_context_builder,
)

SCRIPT_DIR = Path(__file__).resolve().parents[1]
MAVEDB_LLM_METADATA_DIR = MAVEDB_METADATA_OUTPUT_DIR.parent
PROMPT_TEMPLATE_FILE = SCRIPT_DIR / "step1_evidence_extraction_prompt.md"
EXTRACTED_METADATA_OUTPUT_DIR = MAVEDB_LLM_METADATA_DIR / "extracted_metadata"
METADATA_EXTRACTION_LOG_FILE = (
    MAVEDB_LLM_METADATA_DIR / "mavedb_metadata_extraction.log"
)
DEFAULT_BULK_EXCLUDED_PUBLICATION_FILES = frozenset({"10_1101_2024_04_26_591310.md"})


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for bulk MaveDB metadata extraction."""
    parser = argparse.ArgumentParser(
        description="Extract structured metadata from MaveDB-linked publication Markdown.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--extraction-schema",
        required=True,
        help=(
            "Pydantic schema to use for metadata extraction, as "
            "'package.module:SchemaClass' or '/path/to/schema.py:SchemaClass'."
        ),
    )
    parser.add_argument(
        "--publication-full-text-dir",
        type=Path,
        default=FULL_TEXT_MD_DIR,
        help="Directory containing publication Markdown files to process.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=EXTRACTED_METADATA_OUTPUT_DIR,
        help="Directory used to write extracted metadata outputs.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=METADATA_EXTRACTION_LOG_FILE,
        help="Log file used for extraction progress and errors.",
    )
    parser.add_argument(
        "--prompt-template-file",
        type=Path,
        default=PROMPT_TEMPLATE_FILE,
        help="Prompt template file used to build extraction prompts.",
    )
    parser.add_argument(
        "--llm-model",
        default=DEFAULT_LLM_MODEL_NAME,
        help="LLM model ID to use for metadata extraction.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_DOWNLOAD_MAX_WORKERS,
        help="Number of worker threads to use for bulk extraction.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing extracted metadata outputs.",
    )
    parser.add_argument(
        "--create-csv",
        action="store_true",
        help="Create a single CSV from the curated clean metadata JSON outputs.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print and log the full prompt, MaveDB metadata, and publication full text for each item.",
    )
    parser.add_argument(
        "--urn-to-dois-file",
        type=Path,
        default=MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
        help="JSON file mapping MaveDB URNs to primary DOIs.",
    )
    parser.add_argument(
        "--urns",
        nargs="+",
        help=(
            "Filter execution to specific MaveDB URNs (e.g. --urns mavedb:00000093-a-1 "
            "or --urns urn:mavedb:00000093-a-1,urn:mavedb:00000002-a-1)."
        ),
    )
    return parser


def parse_target_urns(
    values: list[str] | tuple[str, ...] | set[str] | str | None,
) -> set[str] | None:
    """Parse user-provided URN string(s) into standardized 'urn:mavedb:...' strings."""
    if not values:
        return None
    raw_list: list[str] = []
    if isinstance(values, str):
        raw_list = [
            item.strip() for item in re.split(r"[,; \s]+", values) if item.strip()
        ]
    elif isinstance(values, (list, tuple, set)):
        for item in values:
            if isinstance(item, str):
                raw_list.extend(
                    [sub.strip() for sub in re.split(r"[,; \s]+", item) if sub.strip()]
                )

    standardized: set[str] = set()
    for raw in raw_list:
        clean = raw.strip()
        if not clean:
            continue
        if clean.startswith("urn_mavedb_"):
            clean = clean.replace("urn_mavedb_", "urn:mavedb:")
        elif clean.startswith("mavedb:"):
            clean = f"urn:{clean}"
        elif not clean.startswith("urn:mavedb:"):
            clean = f"urn:mavedb:{clean}"
        standardized.add(clean)
    return standardized or None


def main() -> None:
    """Run bulk metadata extraction for cached MaveDB URN datasets."""
    args = build_parser().parse_args()
    if args.max_workers < 1:
        raise ValueError("--max-workers must be at least 1")

    extraction_schema = load_extraction_schema(args.extraction_schema)

    urn_to_dois = load_mavedb_urn_to_dois(args.urn_to_dois_file)
    target_urns = parse_target_urns(args.urns)
    if target_urns:
        filtered_urn_to_dois: dict[str, list[str]] = {}
        for target_urn in sorted(target_urns):
            if target_urn in urn_to_dois:
                filtered_urn_to_dois[target_urn] = urn_to_dois[target_urn]
            else:
                urn_stem = format_urn_for_filename(target_urn)
                meta_file = MAVEDB_METADATA_OUTPUT_DIR / f"{urn_stem}.json"
                if meta_file.is_file():
                    entry = json.loads(meta_file.read_text(encoding="utf-8"))
                    dois = get_dois_from_mavedb_entry(entry, log=False) or []
                    filtered_urn_to_dois[target_urn] = dois
                else:
                    filtered_urn_to_dois[target_urn] = []
        urn_to_dois = filtered_urn_to_dois

    bulk_extract_evidence_for_mavedb_urns(
        urn_to_dois=urn_to_dois,
        extraction_schema=extraction_schema,
        output_dir=args.output_dir,
        log_file=args.log_file,
        prompt_template_file=args.prompt_template_file,
        publication_full_text_dir=args.publication_full_text_dir,
        max_workers=args.max_workers,
        overwrite=args.overwrite,
        model_name=args.llm_model,
        create_csv=args.create_csv,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
