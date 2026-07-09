"""CLI runner for extracting structured metadata from MaveDB-linked papers."""

import argparse
from pathlib import Path

from curation_tools.llm_curation.metadata_extraction import (
    DEFAULT_DOWNLOAD_MAX_WORKERS,
    DEFAULT_LLM_MODEL_NAME,
    bulk_extract_metadata_from_publications,
)
from curation_tools.llm_curation.schema_loading import load_extraction_schema
from curation_tools.llm_curation.mavedb.processing import (
    MAVEDB_METADATA_OUTPUT_DIR,
    FULL_TEXT_MD_DIR,
    context_output_suffix_builder,
    format_supplementary_mavedb_metadata,
    output_metadata_builder,
    prompt_context_builder,
)

SCRIPT_DIR = Path(__file__).resolve().parents[1]
MAVEDB_LLM_METADATA_DIR = MAVEDB_METADATA_OUTPUT_DIR.parent
PROMPT_TEMPLATE_FILE = SCRIPT_DIR / "mavedb_metadata_extraction_prompt_template.md"
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
    return parser


def main() -> None:
    """Run bulk metadata extraction for all cached MaveDB publication Markdown files."""
    args = build_parser().parse_args()
    if args.max_workers < 1:
        raise ValueError("--max-workers must be at least 1")

    extraction_schema = load_extraction_schema(args.extraction_schema)

    bulk_extract_metadata_from_publications(
        publication_full_text_paths=[
            str(path) for path in args.publication_full_text_dir.glob("*.md")
        ],
        extraction_schema=extraction_schema,
        output_dir=args.output_dir,
        log_file=args.log_file,
        prompt_template_file=args.prompt_template_file,
        max_workers=args.max_workers,
        overwrite=args.overwrite,
        model_name=args.llm_model,
        create_csv=args.create_csv,
        excluded_publication_files=DEFAULT_BULK_EXCLUDED_PUBLICATION_FILES,
        prompt_context_builder=prompt_context_builder,
        prompt_context_formatter=format_supplementary_mavedb_metadata,
        context_output_suffix_builder=context_output_suffix_builder,
        output_metadata_builder=output_metadata_builder,
    )


if __name__ == "__main__":
    main()
