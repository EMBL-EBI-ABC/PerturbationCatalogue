"""Read-only preflight checks for the ENSG dev DWH pipeline."""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass
from typing import Any

import psycopg2
from elasticsearch import Elasticsearch
from google.api_core.exceptions import Forbidden, NotFound
from google.cloud import bigquery


SQL_IDENTIFIER_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

PRODUCTION_BQ_DATASETS = {"unified_data", "reference"}
PRODUCTION_ES_ALIASES = {
    "dataset-summary",
    "target-summary",
    "landing-page-summary",
}
BQ_TABLES_TO_REPORT = [
    "crispr_data",
    "mave_data",
    "dataset_summary",
    "target_summary",
    "landing_page_summary",
]
PG_TABLES_TO_REPORT = {
    "crispr_data": "crispr_data",
    "mave_data": "mave_data",
    "perturb_seq_dea": "perturb_seq_dea",
    "perturb_seq_gsea": "perturb_seq_gsea",
}
PG_SUMMARY_VIEWS_TO_REPORT = {
    "perturbation": "perturb_seq_summary_perturbation",
    "effect": "perturb_seq_summary_effect",
    "dataset": "perturb_seq_summary_dataset",
}
SYNC_STATE_TABLE = "sync_state"
OPENTARGETS_TARGETS_TABLE = "opentargets_targets"


@dataclass(frozen=True)
class PipelineConfig:
    """Environment-derived pipeline target names."""

    project: str
    bq_location: str
    bq_dataset: str
    pg_conn: str
    pg_connect_timeout: int
    es_url: str
    es_username: str
    es_password: str
    es_index_set: str
    es_indexes: dict[str, str]


def required_env(name: str) -> str:
    """Return a required environment variable."""
    value = os.getenv(name)
    if not value:
        raise RuntimeError(f"{name} is required")
    return value


def env_int(name: str, default: int) -> int:
    """Return an integer environment variable."""
    value = os.getenv(name)
    if not value:
        return default
    try:
        return int(value)
    except ValueError as exc:
        raise RuntimeError(f"{name} must be an integer") from exc


def validate_sql_identifier(value: str, label: str) -> None:
    """Reject unsafe SQL identifier names."""
    if not SQL_IDENTIFIER_RE.fullmatch(value):
        raise RuntimeError(f"{label} must be a plain SQL identifier: {value!r}")


def reject_legacy_migration_name(value: str, label: str) -> None:
    """Reject legacy migration assets that should not be reused."""
    if "gene_id_migration" in value:
        raise RuntimeError(f"{label} must not use legacy gene_id_migration assets")


def require_name_contains(value: str, expected: str, label: str) -> None:
    """Require a dev namespace marker in a configured object name."""
    if expected not in value:
        raise RuntimeError(f"{label} must contain {expected!r}: {value!r}")


def load_config() -> PipelineConfig:
    """Load target names from the environment."""
    es_index_set = os.getenv("ES_INDEX_SET", "")
    return PipelineConfig(
        project=required_env("GCLOUD_PROJECT"),
        bq_location=os.getenv("BQ_LOCATION", "EU"),
        bq_dataset=required_env("BQ_DATASET"),
        pg_conn=required_env("PG_CONN_INTERNAL"),
        pg_connect_timeout=env_int("PG_CONNECT_TIMEOUT", 10),
        es_url=required_env("ES_URL").rstrip("/"),
        es_username=required_env("ES_USERNAME"),
        es_password=required_env("ES_PASSWORD"),
        es_index_set=es_index_set,
        es_indexes={
            "dataset_summary": f"dataset-summary{es_index_set}",
            "target_summary": f"target-summary{es_index_set}",
            "landing_page_summary": f"landing-page-summary{es_index_set}",
        },
    )


def validate_dev_targets(config: PipelineConfig) -> list[str]:
    """Fail if configured destinations look like production or legacy targets."""
    checked: list[str] = []

    bq_targets = {"BQ_DATASET": config.bq_dataset}
    for label, value in bq_targets.items():
        validate_sql_identifier(value, label)
        reject_legacy_migration_name(value, label)

    for label, value in bq_targets.items():
        require_name_contains(value, "ensg_dev", label)
        if value in PRODUCTION_BQ_DATASETS:
            raise RuntimeError(f"{label} points at a production dataset: {value}")
        checked.append(f"{label}={value}")

    checked.append("PG_CONN_INTERNAL=<set>")

    for label, value in config.es_indexes.items():
        reject_legacy_migration_name(value, f"ES {label}")
        checked.append(f"ES {label}={value}")

    return checked


def bq_table_state(
    client: bigquery.Client,
    project: str,
    dataset: str,
    table: str,
    required_access: bool = True,
) -> dict[str, Any]:
    """Return metadata for a BigQuery table without querying row data."""
    table_id = f"{project}.{dataset}.{table}"
    try:
        bq_table = client.get_table(table_id)
    except NotFound:
        return {"exists": False, "table": table_id}
    except Forbidden as exc:
        if required_access:
            raise
        return {
            "accessible": False,
            "exists": None,
            "error": str(exc),
            "table": table_id,
        }
    return {
        "accessible": True,
        "exists": True,
        "table": table_id,
        "rows": bq_table.num_rows,
        "schema_columns": [field.name for field in bq_table.schema],
    }


def check_bq(config: PipelineConfig) -> dict[str, Any]:
    """Report BigQuery source, reference, and dev table state."""
    client = bigquery.Client(project=config.project, location=config.bq_location)
    result: dict[str, Any] = {
        "dev_dataset": config.bq_dataset,
        "reference_dataset": config.bq_dataset,
        "location": config.bq_location,
        "dev_tables": {},
        "reference_table": {},
        "production_baseline": {},
    }

    for table in BQ_TABLES_TO_REPORT:
        result["dev_tables"][table] = bq_table_state(
            client, config.project, config.bq_dataset, table
        )

    result["reference_table"] = bq_table_state(
        client,
        config.project,
        config.bq_dataset,
        OPENTARGETS_TARGETS_TABLE,
    )

    for table in [
        "crispr_data",
        "mave_data",
        "perturb_seq_dea",
        "perturb_seq_gsea",
        "dataset_summary",
        "target_summary",
        "landing_page_summary",
    ]:
        result["production_baseline"][table] = bq_table_state(
            client,
            config.project,
            "unified_data",
            table,
            required_access=False,
        )

    return result


def pg_object_state(cursor: Any, object_name: str) -> dict[str, Any]:
    """Return approximate PostgreSQL object state without scanning rows."""
    cursor.execute("SELECT to_regclass(%s)", (object_name,))
    regclass = cursor.fetchone()[0]
    if regclass is None:
        return {"exists": False, "object": object_name}

    cursor.execute(
        """
        SELECT c.relkind, c.reltuples::bigint
        FROM pg_class c
        JOIN pg_namespace n ON n.oid = c.relnamespace
        WHERE n.nspname = current_schema()
          AND c.relname = %s
        """,
        (object_name,),
    )
    row = cursor.fetchone()
    relkind, reltuples = row if row else (None, None)
    return {
        "exists": True,
        "object": object_name,
        "relkind": relkind,
        "estimated_rows": reltuples,
    }


def check_pg(config: PipelineConfig) -> dict[str, Any]:
    """Report PostgreSQL dev objects and production baselines."""
    result: dict[str, Any] = {
        "dev_tables": {},
        "dev_summary_views": {},
        "sync_state": {},
    }

    with psycopg2.connect(
        config.pg_conn,
        connect_timeout=config.pg_connect_timeout,
    ) as conn:
        with conn.cursor() as cursor:
            for logical_name, object_name in PG_TABLES_TO_REPORT.items():
                result["dev_tables"][logical_name] = pg_object_state(
                    cursor, object_name
                )
            for logical_name, object_name in PG_SUMMARY_VIEWS_TO_REPORT.items():
                result["dev_summary_views"][logical_name] = pg_object_state(
                    cursor, object_name
                )
            result["sync_state"] = pg_object_state(cursor, SYNC_STATE_TABLE)

    return result


def es_alias_state(es: Elasticsearch, alias: str) -> dict[str, Any]:
    """Return Elasticsearch alias/index state using read-only APIs."""
    exists = bool(es.indices.exists_alias(name=alias))
    result: dict[str, Any] = {"exists": exists, "alias": alias}
    if not exists:
        return result

    alias_response = es.indices.get_alias(name=alias)
    result["indices"] = sorted(alias_response.keys())
    try:
        result["docs"] = es.count(index=alias)["count"]
    except Exception as exc:
        result["count_error"] = str(exc)
    return result


def es_index_state(es: Elasticsearch, index: str) -> dict[str, Any]:
    """Return Elasticsearch standalone index state using read-only APIs."""
    exists = bool(es.indices.exists(index=index))
    result: dict[str, Any] = {"exists": exists, "index": index}
    if exists:
        result["docs"] = es.count(index=index)["count"]
    return result


def check_es(config: PipelineConfig) -> dict[str, Any]:
    """Report Elasticsearch dev aliases and production baselines."""
    es = Elasticsearch(
        config.es_url,
        basic_auth=(config.es_username, config.es_password),
        request_timeout=30,
        verify_certs=True,
    )
    target_kind = "indexes" if config.es_index_set else "aliases"
    result: dict[str, Any] = {target_kind: {}, "production_baseline": {}}
    state = es_index_state if config.es_index_set else es_alias_state
    for logical_name, index in config.es_indexes.items():
        result[target_kind][logical_name] = state(es, index)
    for alias in sorted(PRODUCTION_ES_ALIASES):
        result["production_baseline"][alias] = es_alias_state(es, alias)
    return result


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Read-only safety checks for the ENSG dev DWH pipeline."
    )
    parser.add_argument("--skip-bq", action="store_true")
    parser.add_argument("--skip-pg", action="store_true")
    parser.add_argument("--skip-es", action="store_true")
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Only validate configured target names; do not connect to services.",
    )
    return parser.parse_args()


def main() -> None:
    """Run preflight checks and print a JSON report."""
    args = parse_args()
    config = load_config()
    checked_targets = validate_dev_targets(config)

    report: dict[str, Any] = {
        "status": "ok",
        "checked_targets": checked_targets,
        "bq": None,
        "pg": None,
        "es": None,
    }

    if not args.validate_only:
        if not args.skip_bq:
            report["bq"] = check_bq(config)
        if not args.skip_pg:
            report["pg"] = check_pg(config)
        if not args.skip_es:
            report["es"] = check_es(config)

    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
