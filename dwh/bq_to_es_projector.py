#!/usr/bin/env python3
"""
Projector: BigQuery -> Elasticsearch (Python client) for
<contrast|dataset|target|gene>_summary.

Reads rows from:
  <BQ_PROJECT>.<BQ_DATASET>.<contrast|dataset|target|gene>_summary

Writes to ES index:
  <contrast|dataset|target|gene>-summary

Env (to be set *in addition* to dev_secrets or prod_secrets):
  BQ_DATASET=perturb_seq
  BQ_TABLE=<contrast|dataset|target|gene>_summary

Optional:
  ES_INDEX=<contrast|dataset|target|gene>-summary
  BULK_CHUNK_SIZE=2000
  BULK_MAX_RETRIES=5
  BULK_TIMEOUT=120
"""

import os
import sys
import json
import logging
from typing import Any, Dict, Iterable, Tuple

from google.cloud import bigquery
from elasticsearch import Elasticsearch, helpers, ApiError

# ---------------- Config ----------------
BQ_PROJECT = os.getenv("BQ_PROJECT", "prj-ext-dev-pertcat-437314")
BQ_DATASET = os.getenv("BQ_DATASET", "perturb_seq")
BQ_TABLE = os.getenv("BQ_TABLE", "target_summary")

ES_URL = (os.getenv("ES_URL") or "").rstrip("/")
ES_INDEX = os.getenv("ES_INDEX", "target-summary-v1")

ES_API_KEY = os.getenv("ES_API_KEY")
ES_USER = os.getenv("ES_USERNAME")
ES_PASS = os.getenv("ES_PASSWORD")

BULK_CHUNK_SIZE = int(os.getenv("BULK_CHUNK_SIZE", "2000"))
BULK_MAX_RETRIES = int(os.getenv("BULK_MAX_RETRIES", "5"))
BULK_TIMEOUT = int(os.getenv("BULK_TIMEOUT", "120"))

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)


def _get_index_meta(index: str) -> Tuple[str, str]:
    """Returns (`index_name`, `key`) from full index name."""
    index_name = index.split("-")[0]
    if index_name == "contrast":
        key = "contrast_id"
    elif index_name == "dataset":
        key = "dataset_id"
    elif index_name == "target":
        key = "perturbed_target_symbol"
    elif index_name == "gene":
        key = "gene"
    else:
        raise ValueError(f"Unknown index prefix: {index_name} (from {index})")
    return index_name, key


# ------------- ES client ----------------
def make_es_client() -> Elasticsearch:
    if not ES_URL or not ES_USER or not ES_PASS:
        raise RuntimeError("ES credentials are not set")
    es = Elasticsearch(
        ES_URL, basic_auth=(ES_USER, ES_PASS), request_timeout=BULK_TIMEOUT
    )
    return es


def ensure_index(es: Elasticsearch, index: str, index_name: str) -> None:
    """Creates index if it doesn't exist."""
    try:
        if es.indices.exists(index=index):
            logging.info("Index %s exists.", index)
            return

        logging.info("Index '%s' not found, creating...", index)
        settings_file = f"{index_name}-summary_settings+mapping.json"
        try:
            with open(settings_file) as f:
                settings_mapping = json.load(f)
            es.indices.create(index=index, body=settings_mapping)
            logging.info("Index '%s' created with provided settings.", index)
        except FileNotFoundError:
            logging.warning(
                "Settings file '%s' not found. Creating index '%s' with default settings.",
                settings_file,
                index,
            )
            es.indices.create(index=index)
            logging.info("Index '%s' created with default settings.", index)

    except ApiError as e:
        logging.error("Elasticsearch ApiError ensuring index %s: %s", index, e)
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
def get_typed_fields(index_name: str) -> tuple[list[str], list[str], list[str]]:
    int_fields: list[str] = []
    float_fields: list[str] = []
    nested_fields: list[str] = []
    settings_file = f"{index_name}-summary_settings+mapping.json"
    try:
        with open(settings_file) as f:
            settings_mapping = json.load(f)

        mappings = settings_mapping.get("mappings", {})
        properties = mappings.get("properties", {})
        for k, v in properties.items():
            if v.get("type") == "integer":
                int_fields.append(k)
            elif v.get("type") == "float":
                float_fields.append(k)
            elif v.get("type") == "nested":
                nested_fields.append(k)
                nested_properties = v.get("properties", {})
                for n_k, n_v in nested_properties.items():
                    if n_v.get("type") == "integer":
                        int_fields.append(n_k)
                    elif n_v.get("type") == "float":
                        float_fields.append(n_k)
    except FileNotFoundError:
        logging.warning(
            "Settings file '%s' not found. Field types will not be coerced based on schema.",
            settings_file,
        )
    return int_fields, float_fields, nested_fields


def _coerce_num(v, to_float=False):
    if v is None:
        return None
    try:
        return float(v) if to_float else int(v)
    except Exception:
        return v


def transform_row(
    row: Dict[str, Any], key: str, typed_fields: tuple
) -> Tuple[str, Dict[str, Any]]:
    """
    Convert a BQ row to ES document.
    Ensures numerics are numeric and nested arrays are cleaned.
    """
    doc: Dict[str, Any] = {}
    symbol = row.get(key)
    if not symbol:
        raise ValueError(f"Row missing '{key}'")

    numeric_int_fields, numeric_float_fields, nested_fields = typed_fields
    for k, v in row.items():
        if k in numeric_int_fields:
            doc[k] = _coerce_num(v)
        elif k in numeric_float_fields:
            doc[k] = _coerce_num(v, True)
        elif isinstance(v, list):
            if v and isinstance(v[0], dict):
                cleaned = []
                for item in v:
                    if not isinstance(item, dict):
                        continue
                    obj = dict(item)
                    for fields_list, is_float in {
                        numeric_int_fields: False,
                        numeric_float_fields: True,
                    }.items():
                        for nf in fields_list:
                            if nf in obj and obj[nf] is not None:
                                obj[nf] = _coerce_num(obj[nf], is_float)
                    cleaned.append(obj)
                doc[k] = cleaned
            else:
                doc[k] = [x for x in v if x not in (None, "")]
        else:
            doc[k] = v

    return symbol, doc


def actions_generator(
    rows_iter: Iterable[Dict[str, Any]], key: str, typed_fields: tuple
) -> Iterable[Dict[str, Any]]:
    for row in rows_iter:
        _id, doc = transform_row(row, key, typed_fields)
        yield {"_op_type": "index", "_index": ES_INDEX, "_id": _id, "_source": doc}


# ----------------- Main ------------------
def main() -> int:
    if not ES_URL:
        logging.error("ES_URL is not set")
        return 2

    es = make_es_client()
    try:
        index_name, key = _get_index_meta(ES_INDEX)
    except ValueError as e:
        logging.error(str(e))
        return 1

    ensure_index(es, ES_INDEX, index_name)

    typed_fields = get_typed_fields(index_name)

    rows_iter = stream_rows_from_bq(BQ_PROJECT, BQ_DATASET, BQ_TABLE)

    logging.info("Starting bulk indexing into %s …", ES_INDEX)
    success, errors = helpers.bulk(
        es,
        actions_generator(rows_iter, key, typed_fields),
        chunk_size=BULK_CHUNK_SIZE,
        max_retries=BULK_MAX_RETRIES,
        request_timeout=BULK_TIMEOUT,
        raise_on_error=False,  # collect item errors
        stats_only=False,
    )

    if errors:
        # errors is a list of per-item failure dicts (can be large); show first few
        sample = errors[:5] if isinstance(errors, list) else errors
        logging.error(
            "Bulk completed with item errors. Sample: %s",
            json.dumps(sample, indent=2)[:1200],
        )

    logging.info("Bulk done. Successful actions: %s", success)
    return 0


if __name__ == "__main__":
    sys.exit(main())
