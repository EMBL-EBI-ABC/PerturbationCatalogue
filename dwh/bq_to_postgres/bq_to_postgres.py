import argparse
import io
import logging
import os
import uuid
from datetime import timezone
from typing import List, Dict, Tuple, Optional, Any, Set

import psycopg2
from psycopg2 import sql
from google.cloud import bigquery, storage
from tqdm import tqdm
import pyarrow as pa
import pyarrow.csv as pa_csv
import pyarrow.parquet as pq

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# ------------------------------------------------------------------------------
# Configuration & Constants
# ------------------------------------------------------------------------------

TABLES_TO_SYNC = [
    "crispr_data",
    "mave_data",
    "perturb_seq_dea",
    "perturb_seq_gsea",
]

# ponytail: only perturb-seq is partitioned to avoid table clutter from 1000s of small datasets.
PARTITIONED_TABLES = {
    "perturb_seq_dea",
    "perturb_seq_gsea",
}

SYNC_QUERIES = {
    "crispr_data": {
        "export_query": r"""
            SELECT
                dataset_id,
                sample_id,
                perturbed_target_ensg,
                score_name,
                score_value,
                significant,
                significance_criteria,
                ingested_at as max_ingested_at
            FROM `{project}.crispr.data`
            WHERE dataset_id = '{dataset_id}'
        """,
        "ts_query": r"""
            SELECT dataset_id, MAX(ingested_at) as latest_ts, COUNT(*) as row_count
            FROM `{project}.crispr.data`
            GROUP BY dataset_id
        """,
    },
    "mave_data": {
        "export_query": r"""
            SELECT
                dataset_id,
                sample_id,
                perturbed_target_ensg,
                score_name,
                score_value,
                perturbation_name,
                cast(
                    regexp_extract(perturbation_name, r'p\.[a-zA-Z]+(\d+)') as int64
                ) as perturbation_position,
                regexp_extract(perturbation_name, r'p\.([a-zA-Z]+)\d+') as perturbation_aa_wt,
                regexp_extract(
                    perturbation_name, r'p\.[a-zA-Z]+\d+([a-zA-Z=]+)'
                ) as perturbation_aa_change,
                ingested_at as max_ingested_at
            FROM `{project}.mavedb.data`
            WHERE dataset_id = '{dataset_id}'
        """,
        "ts_query": r"""
            SELECT dataset_id, MAX(ingested_at) as latest_ts, COUNT(*) as row_count
            FROM `{project}.mavedb.data`
            GROUP BY dataset_id
        """,
    },
    "perturb_seq_dea": {
        "export_query": r"""
            SELECT
                dataset_id,
                perturbed_target_ensg,
                effect_gene_ensg,
                padj,
                log2foldchange,
                score_name,
                score_value,
                cell_type,
                max_ingested_at
            FROM `{project}.perturb_seq.pertpy_dea`
            WHERE dataset_id = '{dataset_id}'
        """,
        "ts_query": r"""
            SELECT dataset_id, MAX(max_ingested_at) as latest_ts, COUNT(*) as row_count
            FROM `{project}.perturb_seq.pertpy_dea`
            GROUP BY dataset_id
        """,
    },
    "perturb_seq_gsea": {
        "export_query": r"""
            SELECT
                dataset_id,
                term,
                perturbed_target_ensg,
                es,
                nes,
                pval,
                sidak,
                fdr,
                geneset_size,
                leading_edge,
                cell_type,
                max_ingested_at
            FROM `{project}.perturb_seq.pertpy_gsea`
            WHERE dataset_id = '{dataset_id}'
        """,
        "ts_query": r"""
            SELECT dataset_id, MAX(max_ingested_at) as latest_ts, COUNT(*) as row_count
            FROM `{project}.perturb_seq.pertpy_gsea`
            GROUP BY dataset_id
        """,
    },
}

# Defines indexes for each table.
# Format: { table_name: [ (index_name, create_sql_template), ... ] }
INDEX_DEFINITIONS = {
    "perturb_seq_dea": [
        (
            "idx_perturbation_dea",
            "CREATE INDEX {idx} ON {table} (perturbed_target_ensg, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_phenotype_dea",
            "CREATE INDEX {idx} ON {table} (effect_gene_ensg, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_perturbation_phenotype_dea",
            "CREATE INDEX {idx} ON {table} (perturbed_target_ensg, effect_gene_ensg, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_perturb_seq_dea_dataset_id_padj",
            "CREATE INDEX {idx} ON {table} (dataset_id, padj) WHERE effect_gene_ensg IS NOT NULL",
        ),
    ],
    "perturb_seq_gsea": [
        (
            "idx_perturbation_gsea",
            "CREATE INDEX {idx} ON {table} (perturbed_target_ensg, dataset_id, fdr, nes)",
        ),
    ],
    "crispr_data": [
        (
            "idx_crispr_data_dataset",
            "CREATE INDEX {idx} ON {table} (dataset_id)",
        ),
        (
            "idx_crispr_data_target",
            "CREATE INDEX {idx} ON {table} (perturbed_target_ensg)",
        ),
    ],
    "mave_data": [
        (
            "idx_mave_data_dataset",
            "CREATE INDEX {idx} ON {table} (dataset_id)",
        ),
        (
            "idx_mave_data_target",
            "CREATE INDEX {idx} ON {table} (perturbed_target_ensg, dataset_id)",
        ),
    ],
}

# Names of materialized views to refresh for each table.
TABLE_MATERIALIZED_VIEWS = {
    "perturb_seq_dea": [
        "perturb_seq_summary_perturbation",
        "perturb_seq_summary_effect",
        "perturb_seq_summary_dataset",
    ],
}

MATERIALIZED_VIEW_DEFINITIONS = {
    "perturb_seq_summary_perturbation": {
        "create_sql": """
            CREATE MATERIALIZED VIEW perturb_seq_summary_perturbation AS
            SELECT
                dataset_id,
                perturbed_target_ensg,
                COUNT(*) AS n_total,
                COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
                COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up
            FROM perturb_seq_dea
            WHERE padj <= 0.05 AND perturbed_target_ensg IS NOT NULL
            GROUP BY dataset_id, perturbed_target_ensg;
        """,
        "index_sql": "CREATE UNIQUE INDEX idx_perturb_seq_summary_perturbation_pk ON perturb_seq_summary_perturbation (dataset_id, perturbed_target_ensg);",
    },
    "perturb_seq_summary_effect": {
        "create_sql": """
            CREATE MATERIALIZED VIEW perturb_seq_summary_effect AS
            SELECT
                dataset_id,
                effect_gene_ensg,
                COUNT(*) AS n_total,
                COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
                COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up,
                AVG(score_value) AS avg_score
            FROM perturb_seq_dea
            WHERE padj <= 0.05 AND effect_gene_ensg IS NOT NULL
            GROUP BY dataset_id, effect_gene_ensg;
        """,
        "index_sql": "CREATE UNIQUE INDEX idx_perturb_seq_summary_effect_pk ON perturb_seq_summary_effect (dataset_id, effect_gene_ensg);",
    },
    "perturb_seq_summary_dataset": {
        "create_sql": """
            CREATE MATERIALIZED VIEW perturb_seq_summary_dataset AS
            SELECT
                dataset_id,
                COUNT(*) AS n_total
            FROM perturb_seq_dea
            WHERE effect_gene_ensg IS NOT NULL
            GROUP BY dataset_id;
        """,
        "index_sql": "CREATE UNIQUE INDEX idx_perturb_seq_summary_dataset_pk ON perturb_seq_summary_dataset (dataset_id);",
    },
}


def parse_force_pg_tables(raw_tables: Optional[str]) -> Set[str]:
    """Parse comma-separated logical table names that should be fully reloaded."""
    if not raw_tables:
        return set()

    tables = {table.strip() for table in raw_tables.split(",") if table.strip()}
    unknown_tables = tables - set(TABLES_TO_SYNC)
    if unknown_tables:
        valid_tables = ", ".join(TABLES_TO_SYNC)
        unknown = ", ".join(sorted(unknown_tables))
        raise RuntimeError(
            f"Unknown forced PG table(s): {unknown}. Valid values: {valid_tables}"
        )
    return tables


# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------


def get_pg_type(field):
    """Maps BigQuery types to PostgreSQL types."""
    bq_type = field.field_type
    pg_type = {
        "STRING": "TEXT",
        "INTEGER": "BIGINT",
        "FLOAT": "DOUBLE PRECISION",
        "BOOLEAN": "BOOLEAN",
        "TIMESTAMP": "TIMESTAMP WITHOUT TIME ZONE",
        "DATE": "DATE",
    }.get(bq_type, "TEXT")

    if field.mode == "REPEATED":
        pg_type += "[]"

    return pg_type


def get_partition_table_name(table_name: str, dataset_id: str) -> str:
    """Returns a safe, unique, and deterministic Postgres partition table name under 63 bytes."""
    safe_ds = "".join(c if c.isalnum() else "_" for c in dataset_id)
    name = f"{table_name}_{safe_ds}"
    if len(name) > 63:
        import hashlib

        h = hashlib.md5(dataset_id.encode()).hexdigest()[:8]
        name = f"{table_name[:45]}_{safe_ds[:8]}_{h}"
    return name


def check_sync_state_schema_is_valid(cursor) -> bool:
    """Checks if the sync_state table exists and has the expected schema and primary key."""
    cursor.execute(
        "SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'sync_state')"
    )
    if not cursor.fetchone()[0]:
        return False

    cursor.execute(
        """
        SELECT column_name, data_type
        FROM information_schema.columns
        WHERE table_name = 'sync_state'
        ORDER BY column_name;
        """
    )
    cols = {row[0].lower(): row[1].upper() for row in cursor.fetchall()}
    expected = {
        "table_name": "TEXT",
        "dataset_id": "TEXT",
        "last_synced_at": "TIMESTAMP WITHOUT TIME ZONE",
    }
    for col, expected_type in expected.items():
        if col not in cols:
            logging.info(f"sync_state schema mismatch: column {col} missing.")
            return False
        pg_type = cols[col]
        if col in ("table_name", "dataset_id"):
            if "CHAR" not in pg_type and "TEXT" not in pg_type:
                logging.info(
                    f"sync_state schema mismatch: column {col} type {pg_type} is not TEXT."
                )
                return False
        elif col == "last_synced_at":
            if "TIMESTAMP" not in pg_type:
                logging.info(
                    f"sync_state schema mismatch: column {col} type {pg_type} is not TIMESTAMP."
                )
                return False

    cursor.execute(
        """
        SELECT kcu.column_name
        FROM information_schema.table_constraints tc
        JOIN information_schema.key_column_usage kcu
          ON tc.constraint_name = kcu.constraint_name
          AND tc.table_schema = kcu.table_schema
        WHERE tc.table_name = 'sync_state'
          AND tc.constraint_type = 'PRIMARY KEY';
        """
    )
    pk_cols = {row[0].lower() for row in cursor.fetchall()}
    if pk_cols != {"table_name", "dataset_id"}:
        logging.info(
            f"sync_state schema mismatch: primary key columns are {pk_cols} != (table_name, dataset_id)."
        )
        return False

    return True


def ensure_sync_state_exists(cursor):
    """Ensures sync_state table exists with correct schema. If it exists but schema is invalid, drops and recreates it and wipes everything."""
    if not check_sync_state_schema_is_valid(cursor):
        logging.warning(
            "sync_state table does not exist or has an invalid schema. Wiping all synced tables and recreating sync_state."
        )
        for table_name in TABLES_TO_SYNC:
            cursor.execute(
                sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(
                    sql.Identifier(table_name)
                )
            )
        cursor.execute("DROP TABLE IF EXISTS sync_state CASCADE;")
        cursor.execute(
            """
            CREATE TABLE sync_state (
                table_name TEXT NOT NULL,
                dataset_id TEXT NOT NULL,
                last_synced_at TIMESTAMP WITHOUT TIME ZONE,
                PRIMARY KEY (table_name, dataset_id)
            );
            """
        )


def get_all_sync_states(cursor) -> Dict[str, Dict[str, Any]]:
    """Gets the current sync state for all tables and datasets from Postgres."""
    ensure_sync_state_exists(cursor)
    states = {}
    cursor.execute("SELECT table_name, dataset_id, last_synced_at FROM sync_state")
    for table_name, dataset_id, last_synced_at in cursor.fetchall():
        if table_name not in states:
            states[table_name] = {}
        states[table_name][dataset_id] = last_synced_at
    return states


def get_bq_latest_timestamps_and_counts(
    bq_client, bq_dataset, table_name, bq_location
) -> Dict[str, Tuple[Any, int]]:
    """Gets the latest max_ingested_at and row count for every dataset_id in a BQ table."""
    query = SYNC_QUERIES[table_name]["ts_query"].format(
        project=bq_client.project, bq_dataset=bq_dataset
    )
    query_job = bq_client.query(query, location=bq_location)
    results = query_job.result()
    return {row.dataset_id: (row.latest_ts, row.row_count) for row in results}


def export_dataset_to_gcs(
    bq_client, bq_dataset, table_name, bq_location, gcs_bucket, gcs_prefix, dataset_id
) -> str:
    """Exports a specific dataset from BigQuery to GCS as Parquet via a temp table."""
    destination_uri = f"gs://{gcs_bucket}/{gcs_prefix}-*.parquet"
    dataset_ref = bq_client.dataset(bq_dataset)

    logging.info(f"        Exporting {dataset_id} to GCS...")

    # Query to a temporary table to filter by dataset_id
    temp_table_id = f"temp_sync_{uuid.uuid4().hex}"
    temp_table_ref = dataset_ref.table(temp_table_id)

    query = SYNC_QUERIES[table_name]["export_query"].format(
        project=bq_client.project, bq_dataset=bq_dataset, dataset_id=dataset_id
    )
    job_config = bigquery.QueryJobConfig(destination=temp_table_ref)

    query_job = None
    extract_job = None

    try:
        query_job = bq_client.query(query, job_config=job_config, location=bq_location)
        query_job.result()

        # Extract temp table to GCS
        extract_config = bigquery.ExtractJobConfig(destination_format="PARQUET")
        extract_job = bq_client.extract_table(
            temp_table_ref,
            destination_uri,
            location=bq_location,
            job_config=extract_config,
        )
        extract_job.result()
    except BaseException as e:
        if query_job and not query_job.done():
            logging.info("        Cancelling BigQuery query job...")
            query_job.cancel()
        if extract_job and not extract_job.done():
            logging.info("        Cancelling BigQuery extract job...")
            extract_job.cancel()
        raise e
    finally:
        # Clean up temp table
        bq_client.delete_table(temp_table_ref, not_found_ok=True)

    return destination_uri


def _convert_array_column(col):
    """Convert a PyArrow list-typed column to Postgres array literal strings."""
    result = []
    for val in col.to_pylist():
        if val is None:
            result.append(None)
        else:
            escaped = []
            for x in val:
                x_str = str(x).replace("\\", "\\\\").replace('"', '\\"')
                escaped.append(f'"{x_str}"')
            result.append("{" + ",".join(escaped) + "}")
    return pa.array(result, type=pa.string())


def _prepare_table_for_copy(table, bq_schema):
    """Prepare a PyArrow table for Postgres COPY."""
    array_cols = {f.name for f in bq_schema if f.mode == "REPEATED"}

    if array_cols:
        new_columns = []
        for i, name in enumerate(table.column_names):
            col = table.column(i)
            if name in array_cols:
                new_columns.append(_convert_array_column(col))
            else:
                new_columns.append(col)
        table = pa.table(
            {name: col for name, col in zip(table.column_names, new_columns)}
        )

    buf = io.BytesIO()
    write_options = pa_csv.WriteOptions(
        include_header=False,
        delimiter="\t",
    )
    pa_csv.write_csv(table, buf, write_options=write_options)
    buf.seek(0)
    return buf


def cleanup_gcs(gcs_bucket, gcs_prefix, gcs_client):
    """Removes temporary files from GCS."""
    logging.info(f"      - Cleaning up GCS files...")
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_prefix))
    for blob in blobs:
        blob.delete()


def delete_dataset_from_pg(cursor, pg_table, dataset_id):
    """Drops the partition table for a given dataset_id from Postgres (extremely fast delete)."""
    partition_name = get_partition_table_name(pg_table, dataset_id)
    logging.info(
        f"        Dropping partition {partition_name} of {pg_table} if exists..."
    )
    cursor.execute(
        sql.SQL("DROP TABLE IF EXISTS {}").format(sql.Identifier(partition_name))
    )


def update_sync_state(cursor, pg_table, dataset_id, timestamp):
    """Updates or inserts the sync state for a specific dataset."""
    cursor.execute(
        """
        INSERT INTO sync_state (table_name, dataset_id, last_synced_at)
        VALUES (%s, %s, %s)
        ON CONFLICT (table_name, dataset_id) DO UPDATE
        SET last_synced_at = EXCLUDED.last_synced_at;
    """,
        (pg_table, dataset_id, timestamp),
    )


def ensure_pg_table_exists(cursor, pg_table, bq_schema):
    """Ensures the target table exists in Postgres as a partitioned or normal table, creating it if necessary."""
    cursor.execute(
        "SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = %s)",
        (pg_table,),
    )
    if not cursor.fetchone()[0]:
        columns = [f"{field.name} {get_pg_type(field)}" for field in bq_schema]
        if pg_table in PARTITIONED_TABLES:
            logging.info(
                f"Table {pg_table} does not exist. Creating as partitioned by dataset_id."
            )
            cursor.execute(
                sql.SQL("CREATE TABLE {} ({}) PARTITION BY LIST (dataset_id)").format(
                    sql.Identifier(pg_table),
                    sql.SQL(", ").join(map(sql.SQL, columns)),
                )
            )
        else:
            logging.info(
                f"Table {pg_table} does not exist. Creating as standard table."
            )
            cursor.execute(
                sql.SQL("CREATE TABLE {} ({})").format(
                    sql.Identifier(pg_table),
                    sql.SQL(", ").join(map(sql.SQL, columns)),
                )
            )


def ensure_partition_exists(cursor, table_name, dataset_id):
    """Ensures a LIST partition exists on the parent partitioned table for dataset_id."""
    partition_name = get_partition_table_name(table_name, dataset_id)
    cursor.execute(
        "SELECT EXISTS (SELECT FROM pg_class WHERE relname = %s);",
        (partition_name,),
    )
    if not cursor.fetchone()[0]:
        logging.info(
            f"Creating partition {partition_name} of {table_name} for dataset_id='{dataset_id}'"
        )
        cursor.execute(
            sql.SQL(
                "CREATE TABLE IF NOT EXISTS {} PARTITION OF {} FOR VALUES IN (%s)"
            ).format(
                sql.Identifier(partition_name),
                sql.Identifier(table_name),
            ),
            (dataset_id,),
        )


def check_pg_schema_matches_bq(cursor, pg_table, bq_schema) -> bool:
    """Checks if the existing Postgres table schema matches BQ schema."""
    cursor.execute(
        """
        SELECT column_name, data_type
        FROM information_schema.columns
        WHERE table_name = %s
        ORDER BY column_name;
        """,
        (pg_table,),
    )
    pg_cols = {row[0].lower(): row[1].upper() for row in cursor.fetchall()}
    if not pg_cols:
        return False

    expected_cols = {}
    for field in bq_schema:
        col_name = field.name.lower()
        expected_type = get_pg_type(field).upper()
        if "TIMESTAMP" in expected_type:
            expected_type = "TIMESTAMP WITHOUT TIME ZONE"
        elif "DOUBLE" in expected_type:
            expected_type = "DOUBLE PRECISION"
        elif "TEXT[]" in expected_type:
            expected_type = "ARRAY"
        expected_cols[col_name] = expected_type

    if set(pg_cols.keys()) != set(expected_cols.keys()):
        logging.info(f"Schema mismatch for {pg_table}: column sets differ.")
        return False

    for col, expected_type in expected_cols.items():
        pg_type = pg_cols[col]
        if pg_type != expected_type:
            if pg_type == "ARRAY" and "[]" in get_pg_type(
                next(f for f in bq_schema if f.name.lower() == col)
            ):
                continue
            logging.info(
                f"Schema mismatch for {pg_table}.{col}: PG type {pg_type} != BQ type {expected_type}"
            )
            return False

    return True


def check_pg_table_is_partitioned(cursor, pg_table) -> bool:
    """Checks if the existing Postgres table is partitioned by LIST on dataset_id."""
    cursor.execute(
        """
        SELECT pt.partstrat, pg_get_partkeydef(c.oid)
        FROM pg_partitioned_table pt
        JOIN pg_class c ON pt.partrelid = c.oid
        WHERE c.relname = %s;
        """,
        (pg_table,),
    )
    row = cursor.fetchone()
    if not row:
        logging.info(f"Table {pg_table} is not partitioned in Postgres.")
        return False
    strat, key_def = row
    is_valid = strat == "l" and "dataset_id" in key_def.lower()
    if not is_valid:
        logging.info(
            f"Table {pg_table} is partitioned, but not by LIST on dataset_id (strat: {strat}, key: {key_def})."
        )
    return is_valid


def drop_dependent_views_and_indexes(cursor, table_name):
    """Finds and drops all materialized views and indexes depending on the table."""
    # 1. Find all dependent materialized views
    cursor.execute(
        """
        SELECT DISTINCT dependee.relname AS mv_name
        FROM pg_depend d
        JOIN pg_rewrite r ON d.objid = r.oid
        JOIN pg_class dependee ON r.ev_class = dependee.oid
        JOIN pg_class dependent ON d.refobjid = dependent.oid
        WHERE dependent.relname = %s AND dependee.relkind = 'm';
        """,
        (table_name,),
    )
    mvs = [row[0] for row in cursor.fetchall()]
    for mv in mvs:
        logging.info(f"    Dropping dependent materialized view {mv}...")
        cursor.execute(
            sql.SQL("DROP MATERIALIZED VIEW IF EXISTS {} CASCADE").format(
                sql.Identifier(mv)
            )
        )

    # 2. Find and drop all indexes
    cursor.execute(
        """
        SELECT c.relname AS index_name
        FROM pg_index i
        JOIN pg_class c ON c.oid = i.indexrelid
        JOIN pg_class t ON t.oid = i.indrelid
        WHERE t.relname = %s AND i.indisprimary = FALSE;
        """,
        (table_name,),
    )
    indexes = [row[0] for row in cursor.fetchall()]
    for idx in indexes:
        logging.info(f"    Dropping index {idx}...")
        cursor.execute(sql.SQL("DROP INDEX IF EXISTS {}").format(sql.Identifier(idx)))


def create_defined_indexes_and_views(cursor, table_name):
    """Recreates only the indexes and materialized views defined in DWH pipeline configs."""
    # 1. Recreate indexes
    if table_name in INDEX_DEFINITIONS:
        logging.info(f"    Creating defined indexes for {table_name}...")
        for index_name, index_sql_template in INDEX_DEFINITIONS[table_name]:
            logging.info(f"        Creating index {index_name}...")
            index_sql = index_sql_template.format(
                idx=sql.Identifier(index_name).as_string(cursor.connection),
                table=sql.Identifier(table_name).as_string(cursor.connection),
            )
            cursor.execute(index_sql)

    # 2. Recreate materialized views
    if table_name in TABLE_MATERIALIZED_VIEWS:
        logging.info(f"    Creating defined materialized views for {table_name}...")
        for mv_name in TABLE_MATERIALIZED_VIEWS[table_name]:
            if mv_name in MATERIALIZED_VIEW_DEFINITIONS:
                logging.info(f"        Creating materialized view {mv_name}...")
                mv_def = MATERIALIZED_VIEW_DEFINITIONS[mv_name]
                cursor.execute(mv_def["create_sql"])
                if mv_def.get("index_sql"):
                    cursor.execute(mv_def["index_sql"])


# ------------------------------------------------------------------------------
# Core Logic Class
# ------------------------------------------------------------------------------


class TableSynchronizer:
    def __init__(
        self,
        drop_and_recreate_indexes: bool,
        pg_conn_str: str,
        bq_dataset: str,
        bq_location: str,
        gcs_bucket: str,
    ):
        self.drop_and_recreate_indexes = drop_and_recreate_indexes
        self.pg_conn_str = pg_conn_str
        self.bq_dataset = bq_dataset
        self.bq_location = bq_location
        self.gcs_bucket = gcs_bucket
        self.bq_client = bigquery.Client()
        self.gcs_client = storage.Client()

    def sync_table(self, table_name: str, plan: Dict[str, Any]):
        logging.info(f"Syncing table {table_name}...")
        try:
            self._sync_unified(table_name, plan)
            logging.info(f"Successfully synced {table_name}.")

        except BaseException as e:
            logging.error(f"Error syncing {table_name}: {e}.")
            raise e

    def _get_bq_schema(self, table_name):
        query = (
            SYNC_QUERIES[table_name]["export_query"].format(
                project=self.bq_client.project,
                bq_dataset=self.bq_dataset,
                dataset_id="_dummy_limit_0_",
            )
            + " LIMIT 0"
        )
        job = self.bq_client.query(query, location=self.bq_location)
        result = job.result()
        return result.schema

    def _export_and_prepare_dataset(self, table_name, ds_id, bq_schema):
        """Exports a dataset to GCS, downloads and converts shards to a single TSV file on disk, and cleans up GCS."""
        gcs_prefix = f"tmp/{table_name}/{ds_id}/{uuid.uuid4().hex}"
        import tempfile

        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        try:
            export_dataset_to_gcs(
                self.bq_client,
                self.bq_dataset,
                table_name,
                self.bq_location,
                self.gcs_bucket,
                gcs_prefix,
                ds_id,
            )

            bucket = self.gcs_client.get_bucket(self.gcs_bucket)
            blobs = list(bucket.list_blobs(prefix=gcs_prefix))

            if not blobs:
                tmp_file.close()
                try:
                    os.unlink(tmp_file.name)
                except Exception:
                    pass
                return ds_id, None

            for blob in blobs:
                with tempfile.NamedTemporaryFile() as shard_file:
                    blob.download_to_file(shard_file)
                    shard_file.seek(0)
                    table = pq.read_table(shard_file.name)

                    chunk_buf = _prepare_table_for_copy(table, bq_schema)
                    tmp_file.write(chunk_buf.getvalue())
                    chunk_buf.close()

            tmp_file.seek(0)
            return ds_id, tmp_file
        except BaseException as e:
            tmp_file.close()
            if os.path.exists(tmp_file.name):
                try:
                    os.unlink(tmp_file.name)
                except Exception:
                    pass
            raise e
        finally:
            cleanup_gcs(self.gcs_bucket, gcs_prefix, self.gcs_client)

    def _load_tsv_to_pg(self, cursor, pg_table, ds_id, tsv_file):
        """Loads a prepared TSV file into the Postgres table (partition or standard table) for ds_id."""
        if tsv_file is None:
            return

        try:
            if pg_table in PARTITIONED_TABLES:
                delete_dataset_from_pg(cursor, pg_table, ds_id)
                ensure_partition_exists(cursor, pg_table, ds_id)

            copy_sql = sql.SQL(
                "COPY {} FROM STDIN WITH (FORMAT CSV, DELIMITER E'\\t', QUOTE '\"', NULL '')"
            ).format(sql.Identifier(pg_table))

            cursor.copy_expert(copy_sql, tsv_file)
        finally:
            tsv_file.close()
            if hasattr(tsv_file, "name") and os.path.exists(tsv_file.name):
                try:
                    os.unlink(tsv_file.name)
                except Exception:
                    pass

    def _sync_unified(self, table_name: str, plan: Dict[str, Any]):
        """
        Unified Sync Mode (strictly sequential):
        - Single transaction for:
            1. Drop dependent materialized views and indexes (if requested or mandatory)
            2. Data Updates (Delete old rows/partition + create empty partition if partitioned + COPY) for all datasets
            3. Recreate defined indexes and materialized views
            4. Update sync_state
        """
        bq_schema = self._get_bq_schema(table_name)
        datasets = sorted(plan["to_update"] + plan["to_insert"])

        # Main Transaction Block
        with psycopg2.connect(self.pg_conn_str) as conn:
            with conn.cursor() as cursor:
                ensure_pg_table_exists(cursor, table_name, bq_schema)

                logging.info(
                    f"    Starting transaction for {len(datasets)} datasets..."
                )

                # A. Drop dependent Views and Indexes
                drop_dependent_views_and_indexes(cursor, table_name)

                # B. Ingest Loop (strictly sequential)
                logging.info("    Running sequential dataset preparation and load...")

                pbar = tqdm(
                    datasets,
                    desc=f"    Syncing {table_name}",
                    unit="dataset",
                )
                for ds_id in pbar:
                    try:
                        _, tsv_buf = self._export_and_prepare_dataset(
                            table_name, ds_id, bq_schema
                        )
                        self._load_tsv_to_pg(cursor, table_name, ds_id, tsv_buf)
                        bq_ts = plan["bq_info"][ds_id][0]
                        update_sync_state(cursor, table_name, ds_id, bq_ts)
                    except Exception as e:
                        logging.error(f"Failed to process dataset {ds_id}: {e}")
                        raise e
                pbar.close()

                # C. Recreate only defined Indexes and Materialized Views
                create_defined_indexes_and_views(cursor, table_name)


# ------------------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="BigQuery to PostgreSQL Sync Script")
    parser.add_argument(
        "--bq-dataset",
        default=os.getenv("BQ_DATASET"),
        help="BigQuery dataset name (env: BQ_DATASET)",
    )
    parser.add_argument(
        "--bq-location",
        default=os.getenv("BQ_LOCATION"),
        help="BigQuery location (env: BQ_LOCATION)",
    )
    parser.add_argument(
        "--pg-conn",
        default=os.getenv("PG_CONN"),
        help="PostgreSQL connection string (env: PG_CONN)",
    )
    parser.add_argument(
        "--gcs-bucket",
        default=os.getenv("GCLOUD_TMP_BUCKET"),
        help="GCS bucket for temporary files (env: GCLOUD_TMP_BUCKET)",
    )
    parser.add_argument(
        "--drop-and-recreate-indexes",
        action="store_true",
        help="Drop indexes before ingestion and recreate them afterwards (in the same transaction).",
    )
    parser.add_argument(
        "--force-pg-tables",
        default="",
        help="Comma-separated logical table names to fully reload even when sync_state is current.",
    )

    args = parser.parse_args()
    force_pg_tables = parse_force_pg_tables(args.force_pg_tables)

    # Validate required arguments
    missing = []
    if not args.bq_dataset:
        missing.append("--bq-dataset / BQ_DATASET")
    if not args.bq_location:
        missing.append("--bq-location / BQ_LOCATION")
    if not args.pg_conn:
        missing.append("--pg-conn / PG_CONN")
    if not args.gcs_bucket:
        missing.append("--gcs-bucket / GCLOUD_TMP_BUCKET")
    if missing:
        parser.error("The following arguments are required: " + ", ".join(missing))

    # Initialize synchronizer
    synchronizer = TableSynchronizer(
        drop_and_recreate_indexes=args.drop_and_recreate_indexes,
        pg_conn_str=args.pg_conn,
        bq_dataset=args.bq_dataset,
        bq_location=args.bq_location,
        gcs_bucket=args.gcs_bucket,
    )

    # Database connection for planning
    conn = psycopg2.connect(args.pg_conn)
    conn.autocommit = True

    # ponytail: unconditionally terminate any other backend processes to prevent hanging locks.
    logging.info("Terminating other database sessions to prevent locks...")
    with conn.cursor() as cursor:
        cursor.execute(
            "SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE pid != pg_backend_pid() AND datname = current_database();"
        )

    try:
        with conn.cursor() as cursor:
            # fetch current state
            pg_states = get_all_sync_states(cursor)
            bq_client = synchronizer.bq_client

            for table_name in TABLES_TO_SYNC:
                logging.info(f"Checking {table_name}...")

                # 1. Get BQ state
                try:
                    bq_info = get_bq_latest_timestamps_and_counts(
                        bq_client, args.bq_dataset, table_name, args.bq_location
                    )
                except Exception as e:
                    logging.warning(
                        f"Could not fetch BQ info for {table_name}: {e}. Skipping."
                    )
                    continue

                # Schema and partition validations
                bq_schema = synchronizer._get_bq_schema(table_name)
                schema_matches = check_pg_schema_matches_bq(
                    cursor, table_name, bq_schema
                )
                is_partitioned = check_pg_table_is_partitioned(cursor, table_name)
                should_be_partitioned = table_name in PARTITIONED_TABLES

                # Determine if hard reset is required
                hard_reset_required = not schema_matches or (
                    is_partitioned != should_be_partitioned
                )

                # 2. Compare with PG state
                pg_info = pg_states.get(table_name, {})
                to_insert = []
                to_update = []

                if not hard_reset_required and table_name not in PARTITIONED_TABLES:
                    # Check if any dataset is new or updated compared to pg_info
                    any_changes = False
                    for ds_id, (bq_ts, bq_count) in bq_info.items():
                        if ds_id not in pg_info:
                            any_changes = True
                            break
                        pg_ts = pg_info[ds_id]
                        if pg_ts and pg_ts.tzinfo is None:
                            pg_ts = pg_ts.replace(tzinfo=timezone.utc)
                        if bq_ts > pg_ts:
                            any_changes = True
                            break
                    if any_changes:
                        logging.info(
                            f"  Changes detected in non-partitioned {table_name}. Triggering full table reload."
                        )
                        hard_reset_required = True

                if hard_reset_required:
                    logging.warning(
                        f"  Hard reset triggered for {table_name}: schema matches BQ={schema_matches}, partitioned={is_partitioned} (should be partitioned={should_be_partitioned})."
                    )
                    # Dropping table with CASCADE drops it completely along with partitions and dependent views
                    cursor.execute(
                        sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(
                            sql.Identifier(table_name)
                        )
                    )
                    # Clear sync state records for this table
                    cursor.execute(
                        "DELETE FROM sync_state WHERE table_name = %s", (table_name,)
                    )
                    # Reimport all datasets
                    to_insert = list(bq_info.keys())
                elif table_name in force_pg_tables:
                    logging.info(
                        f"  Force reload requested for {table_name}; all BigQuery datasets will be reloaded."
                    )
                    to_update = list(bq_info)
                else:
                    for ds_id, (bq_ts, bq_count) in bq_info.items():
                        if ds_id not in pg_info:
                            to_insert.append(ds_id)
                        else:
                            pg_ts = pg_info[ds_id]
                            if pg_ts and pg_ts.tzinfo is None:
                                pg_ts = pg_ts.replace(tzinfo=timezone.utc)
                            if bq_ts > pg_ts:
                                to_update.append(ds_id)

                if not to_insert and not to_update:
                    logging.info(f"  {table_name} is up to date.")
                    continue

                # 3. Create plan
                plan = {
                    "to_insert": to_insert,
                    "to_update": to_update,
                    "bq_info": bq_info,
                }

                logging.info(f"  Plan for {table_name}:")
                logging.info(f"    To Insert: {len(to_insert)}")
                logging.info(f"    To Update: {len(to_update)}")

                # 4. Execute Sync
                synchronizer.sync_table(table_name, plan)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
