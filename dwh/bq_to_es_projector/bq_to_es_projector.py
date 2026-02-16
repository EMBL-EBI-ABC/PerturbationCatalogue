#!/usr/bin/env python3
"""
Projector: BigQuery -> Elasticsearch (Python client).

Reads rows from:
  <BQ_PROJECT>.<BQ_DATASET>.<dataset_summary|target_summary|landing_page_summary>

Writes to ES index:
  <dataset-summary|target-summary|landing-page-summary>
"""

import os
import sys
import json
import logging
import datetime
from typing import Any, Dict, Iterable, Tuple

from google.cloud import bigquery
from elasticsearch import Elasticsearch, helpers, ApiError
from tqdm import tqdm

# ---------------- Config ----------------
BQ_PROJECT = os.getenv("GCLOUD_PROJECT")
BQ_DATASET = os.getenv("BQ_DATASET")

TABLE_CONFIG = {
    "dataset_summary": {
        "index_base": "dataset-summary",
        "key_field": "dataset_id",
        "prefix": "dataset",
    },
    "target_summary": {
        "index_base": "target-summary",
        "key_field": "perturbed_target_symbol",
        "prefix": "target",
    },
    "landing_page_summary": {
        "index_base": "landing-page-summary",
        "key_field": "summary",
        "prefix": "landing-page",
    },
}

ES_URL = (os.getenv("ES_URL") or "").rstrip("/")

ES_USER = os.getenv("ES_USERNAME")
ES_PASS = os.getenv("ES_PASSWORD")

BULK_CHUNK_SIZE = int(os.getenv("BULK_CHUNK_SIZE", "2000"))
BULK_MAX_RETRIES = int(os.getenv("BULK_MAX_RETRIES", "5"))
BULK_TIMEOUT = int(os.getenv("BULK_TIMEOUT", "120"))

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

# Silence noisy dependency logs
logging.getLogger("elastic_transport").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


# ------------- ES client ----------------
def make_es_client() -> Elasticsearch:
    if not ES_URL or not ES_USER or not ES_PASS:
        raise RuntimeError("ES credentials are not set")
    es = Elasticsearch(
        ES_URL,
        basic_auth=(ES_USER, ES_PASS),
        request_timeout=BULK_TIMEOUT,
        verify_certs=True,
    )
    return es


def get_index_metadata(index_base: str) -> Tuple[str, str]:
    """
    Returns (index_prefix, key_field_name) for a given ES base index name.
    """
    for cfg in TABLE_CONFIG.values():
        if cfg["index_base"] == index_base:
            return cfg["prefix"], cfg["key_field"]
    raise ValueError(f"Unknown index base: {index_base}")


def generate_dataset_summary_mapping() -> Dict[str, Any]:
    """Dynamically generate mapping for dataset-summary from be/dataset_metadata.json"""
    mapping = {
        "settings": {
            "analysis": {
                "analyzer": {"en": {"type": "standard", "stopwords": "_english_"}},
                "normalizer": {
                    "lc_ascii": {
                        "type": "custom",
                        "char_filter": [],
                        "filter": ["lowercase", "asciifolding"],
                    }
                },
            }
        },
        "mappings": {"properties": {}},
    }

    script_dir = os.path.dirname(os.path.abspath(__file__))
    metadata_path = os.path.join(script_dir, "../../be/dataset_metadata.json")
    if not os.path.exists(metadata_path):
        metadata_path = "dataset_metadata.json"

    with open(metadata_path, "r") as f:
        meta = json.load(f)

    for field in meta["fields"]:
        es_field = field["es_field"]
        es_type = field.get("es_type", "keyword")

        prop = {"type": es_type}
        if es_type == "keyword":
            prop["normalizer"] = "lc_ascii"
        elif es_type == "text":
            prop["analyzer"] = "en"

        mapping["mappings"]["properties"][es_field] = prop

    # Explicitly define data_modalities with keyword and text sub-field
    mapping["mappings"]["properties"]["data_modalities"] = {
        "type": "keyword",
        "ignore_above": 256,
        "normalizer": "lc_ascii",
        "fields": {"text": {"type": "text", "analyzer": "en"}},
    }

    return mapping


def get_mapping(index_prefix: str) -> Dict[str, Any]:
    if index_prefix == "dataset":
        return generate_dataset_summary_mapping()
    mapping_file = f"{index_prefix}-summary_settings+mapping.json"
    if not os.path.exists(mapping_file):
        raise FileNotFoundError(f"Mapping file not found: {mapping_file}")
    with open(mapping_file, "r") as f:
        return json.load(f)


def ensure_index(es: Elasticsearch, index: str, prefix: str) -> None:
    try:
        if es.indices.exists(index=index):
            logging.info("Index %s exists.", index)
            return

        logging.info("Index %s does not exist. Creating...", index)
        mapping_body = get_mapping(prefix)

        es.indices.create(index=index, body=mapping_body)
        logging.info("Index %s created successfully.", index)

    except ApiError as e:
        logging.error(f"Elasticsearch ApiError: {e}")
        raise


# ------------- BQ helpers ---------------
def stream_rows_from_bq(
    project: str, dataset: str, table: str
) -> Iterable[Dict[str, Any]]:
    client = bigquery.Client(project=project)
    table_ref = f"{project}.{dataset}.{table}"
    logging.info("Reading BigQuery rows: %s", table_ref)
    for row in client.list_rows(table_ref):
        yield dict(row)


# --------- Transform / Actions ----------
def get_typed_fields(index_prefix: str) -> tuple[list[str], list[str], list[str]]:
    int_fields: list[str] = []
    float_fields: list[str] = []
    nested_fields: list[str] = []

    settings_mapping = get_mapping(index_prefix)
    for k, v in settings_mapping["mappings"]["properties"].items():
        if v["type"] == "integer":
            int_fields.append(k)
        elif v["type"] == "float":
            float_fields.append(k)
        elif v["type"] == "nested":
            nested_fields.append(k)
            for n_k, n_v in v["properties"].items():
                if n_v["type"] == "integer":
                    int_fields.append(n_k)
                elif n_v["type"] == "float":
                    float_fields.append(n_k)
    return int_fields, float_fields, nested_fields


def _coerce_num(v, to_float=False):
    if v is None:
        return None
    try:
        return float(v) if to_float else int(v)
    except Exception:
        return v


def transform_row(
    row: Dict[str, Any], es_index: str, prefix: str, key_field: str
) -> Tuple[str, Dict[str, Any]]:
    """
    Convert a BQ row to ES document.
    Ensures numerics are numeric and nested arrays are cleaned.
    """
    doc: Dict[str, Any] = {}

    if key_field != "summary":
        symbol = row.get(key_field)
        if not symbol:
            raise ValueError(f"Row missing '{key_field}'")
    else:
        symbol = "summary"

    numeric_int_fields, numeric_float_fields, nested_fields = get_typed_fields(prefix)
    for k, v in row.items():
        if k in numeric_int_fields:
            doc[k] = _coerce_num(v)
        elif k in numeric_float_fields:
            doc[k] = _coerce_num(v, True)
        elif k in nested_fields and isinstance(v, list):
            # Handle nested objects: clean and coerce numeric fields
            cleaned = []
            for item in v:
                if not isinstance(item, dict):
                    continue
                obj = dict(item)
                for fields_list, is_float in (
                    (numeric_int_fields, False),
                    (numeric_float_fields, True),
                ):
                    for nf in fields_list:
                        if nf in obj and obj[nf] is not None:
                            obj[nf] = _coerce_num(obj[nf], is_float)
                cleaned.append(obj)
            doc[k] = cleaned
        elif isinstance(v, list):
            # Handle regular lists: filter out None and empty strings
            doc[k] = [x for x in v if x not in (None, "")]
        else:
            doc[k] = v

    return symbol, doc


def actions_generator(
    rows_iter: Iterable[Dict[str, Any]], es_index: str, prefix: str, key_field: str
) -> Iterable[Dict[str, Any]]:
    for row in rows_iter:
        try:
            _id, doc = transform_row(row, es_index, prefix, key_field)
            yield {"_op_type": "index", "_index": es_index, "_id": _id, "_source": doc}
        except ValueError:
            pass


def prune_old_indexes(es: Elasticsearch) -> None:
    """
    Prune indexes for the summary families.
    Keep live version (pointed by alias) + up to 2 older versions.
    """
    logging.info("Starting index pruning...")
    index_bases = [cfg["index_base"] for cfg in TABLE_CONFIG.values()]

    for base in index_bases:
        pattern = f"*-{base}"
        try:
            indices = es.indices.get(index=pattern).body
        except ApiError:
            logging.info("No indexes found for %s", base)
            continue

        # Sort indices by name (YYYY-MM-DD prefix) descending
        index_names = sorted(indices.keys(), reverse=True)
        if not index_names:
            continue

        # Find the live index (pointed to by the alias)
        try:
            alias_info = es.indices.get_alias(name=base).body
            live_index = list(alias_info.keys())[0] if alias_info else None
        except ApiError:
            live_index = None

        logging.info("Index %s: live version is %s", base, live_index)

        if not live_index:
            continue

        # Since live is the latest, we keep it + 2 previous versions (first 3 in sorted list)
        to_keep = index_names[:3]
        to_delete = [idx for idx in index_names if idx not in to_keep]

        for idx in to_delete:
            logging.info("Deleting old index: %s", idx)
            es.indices.delete(index=idx)


def main() -> int:
    if not ES_URL:
        logging.error("ES_URL is not set")
        return 2

    date_str = datetime.datetime.now().strftime("%Y-%m-%d")
    es = make_es_client()

    sync_results = {}  # index_base -> new_index

    for table, cfg in TABLE_CONFIG.items():
        base = cfg["index_base"]
        prefix = cfg["prefix"]
        key_field = cfg["key_field"]

        es_index = f"{date_str}-{base}"

        logging.info("Starting sync for %s -> %s", table, es_index)

        try:
            # If current date index already exists, delete it first
            if es.indices.exists(index=es_index):
                logging.info(
                    "Index %s already exists for today. Deleting it first...", es_index
                )
                es.indices.delete(index=es_index)

            ensure_index(es, es_index, prefix)

            # Use BigQuery client to get total row count for progress bar
            bq_client = bigquery.Client(project=BQ_PROJECT)
            table_ref = f"{BQ_PROJECT}.{BQ_DATASET}.{table}"
            bq_table = bq_client.get_table(table_ref)
            total_rows = bq_table.num_rows

            rows_iter = stream_rows_from_bq(BQ_PROJECT, BQ_DATASET, table)

            logging.info(
                "Starting bulk indexing into %s (total: %s) …", es_index, total_rows
            )

            pbar = tqdm(total=total_rows, desc=table, unit="rows")

            def actions_with_progress():
                for action in actions_generator(rows_iter, es_index, prefix, key_field):
                    pbar.update(1)
                    yield action

            success, errors = helpers.bulk(
                es,
                actions_with_progress(),
                chunk_size=BULK_CHUNK_SIZE,
                max_retries=BULK_MAX_RETRIES,
                request_timeout=BULK_TIMEOUT,
                raise_on_error=False,
                stats_only=False,
            )
            pbar.close()

            if errors:
                sample = errors[:5] if isinstance(errors, list) else errors
                logging.error(
                    "Bulk completed with item errors. Sample: %s",
                    json.dumps(sample, indent=2)[:1200],
                )

            logging.info("Sync for %s done. Successful: %s", table, success)
            sync_results[base] = es_index

        except Exception as e:
            logging.error("Failed to sync table %s: %s", table, e)
            return 1

    # If all successful, move aliases
    if len(sync_results) == len(TABLE_CONFIG):
        logging.info("All tables synced successfully. Moving aliases...")
        actions = []
        for base, new_index in sync_results.items():
            # Remove existing alias from any old indexes
            try:
                old_indices = es.indices.get_alias(name=base).body
                for old_idx in old_indices:
                    actions.append({"remove": {"index": old_idx, "alias": base}})
            except ApiError:
                pass

            # Add new alias
            actions.append({"add": {"index": new_index, "alias": base}})

        if actions:
            es.indices.update_aliases(body={"actions": actions})
            logging.info("Aliases updated successfully.")

        # Prune old indexes
        prune_old_indexes(es)

    return 0


if __name__ == "__main__":
    sys.exit(main())
