"""Shared locations for study curation code, source caches, and run history."""

from pathlib import Path

STUDY_CURATION_DIR = Path(__file__).resolve().parent
REPO_ROOT = STUDY_CURATION_DIR.parents[2]
LLM_DIR = STUDY_CURATION_DIR / "llm"
PROMPTS_DIR = LLM_DIR / "prompts"
DEFAULT_SCHEMA_PATH = LLM_DIR / "llm_curation_schema.py"
DEFAULT_STEP1_PROMPT_PATH = PROMPTS_DIR / "step1_evidence_extraction_prompt.md"
DEFAULT_STEP2_PROMPT_PATH = PROMPTS_DIR / "step2_specific_term_extraction.md"
DEFAULT_STEP3_PROMPT_PATH = PROMPTS_DIR / "step3_candidate_discovery_prompt.md"

MAVEDB_DIR = REPO_ROOT / "data_exploration" / "MaveDB"
MAVEDB_DUMP_DIR = MAVEDB_DIR / "Dump" / "mavedb-dump.20250612164404" / "csv"
MAVEDB_RESOURCE_DIR = STUDY_CURATION_DIR / "mavedb" / "resources"
DEFAULT_MAPPING_PATH = MAVEDB_RESOURCE_DIR / "metadata_mappings.csv"
DEFAULT_OVERRIDE_PATH = MAVEDB_RESOURCE_DIR / "metadata_overrides.csv"

CURATION_CACHE_DIR = REPO_ROOT / "data_exploration" / "curation_cache"
PAPERSCRAPER_FULL_TEXT_RAW_DIR = CURATION_CACHE_DIR / "publications" / "raw"
FULL_TEXT_MD_DIR = CURATION_CACHE_DIR / "publications" / "markdown"
DOWNLOAD_PROGRESS_LOG_FILE = CURATION_CACHE_DIR / "publications" / "download.log"
MAVEDB_METADATA_OUTPUT_DIR = CURATION_CACHE_DIR / "mavedb" / "metadata"
MAVEDB_URN_TO_DOIS_OUTPUT_FILE = CURATION_CACHE_DIR / "mavedb" / "urn_to_dois.json"
MAVEDB_DOI_TO_FULLTEXT_OUTPUT_FILE = (
    CURATION_CACHE_DIR / "mavedb" / "doi_to_fulltext.json"
)
CURATION_RUNS_DIR = REPO_ROOT / "curation_runs"
