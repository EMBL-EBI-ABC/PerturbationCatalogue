"""Elasticsearch index alias names used by the backend."""

import os
import re


ES_INDEX_SUFFIX = os.getenv("ES_INDEX_SUFFIX", "")


def validate_es_index_suffix(suffix: str) -> str:
    if suffix and not re.fullmatch(r"[A-Za-z0-9_-]+", suffix):
        raise ValueError(
            "ES_INDEX_SUFFIX may only contain letters, numbers, underscores, and hyphens"
        )
    return suffix


def apply_es_index_suffix(index_base: str) -> str:
    suffix = validate_es_index_suffix(ES_INDEX_SUFFIX)
    return f"{index_base}{suffix}"


ES_LANDING_PAGE_SUMMARY = apply_es_index_suffix("landing-page-summary")
ES_TARGET_SUMMARY = apply_es_index_suffix("target-summary")
ES_DATASET_SUMMARY = apply_es_index_suffix("dataset-summary")
