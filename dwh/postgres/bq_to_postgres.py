import argparse
import concurrent.futures
import contextlib
import enum
import io
import logging
import os
import sys
import uuid
import time
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict, Tuple, Optional, Any

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

# Defines indexes for each table.
# Format: { table_name: [ (index_name, create_sql_template), ... ] }
INDEX_DEFINITIONS = {
    "perturb_seq_dea": [
        (
            "idx_perturbation_dea",
            "CREATE INDEX {idx} ON {table} (perturbed_target_symbol, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_phenotype_dea",
            "CREATE INDEX {idx} ON {table} (gene, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_perturbation_phenotype_dea",
            "CREATE INDEX {idx} ON {table} (perturbed_target_symbol, gene, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_perturb_seq_dea_dataset_id_padj",
            "CREATE INDEX {idx} ON {table} (dataset_id, padj) WHERE gene IS NOT NULL",
        ),
    ],
    "perturb_seq_gsea": [
        (
            "idx_perturbation_gsea",
            "CREATE INDEX {idx} ON {table} (perturbed_target_symbol, dataset_id, fdr, nes)",
        ),
    ],
    "crispr_data": [
        (
            "idx_crispr_data_dataset",
            "CREATE INDEX {idx} ON {table} (dataset_id)",
        ),
        (
            "idx_crispr_data_target",
            "CREATE INDEX {idx} ON {table} (perturbed_target_symbol)",
        ),
    ],
    "mave_data": [
        (
            "idx_mave_data_dataset",
            "CREATE INDEX {idx} ON {table} (dataset_id)",
        ),
        (
            "idx_mave_data_target",
            "CREATE INDEX {idx} ON {table} (perturbed_target_symbol, dataset_id)",
        ),
    ],
}

# Materialized view definitions.
# Format: { base_table_name: [ (view_name, create_sql_template, [(index_name, index_sql_template), ...]), ... ] }
MATERIALIZED_VIEW_DEFINITIONS = {
    "perturb_seq_dea": [
        (
            "perturb_seq_summary_perturbation",
            """CREATE MATERIALIZED VIEW {view} AS
SELECT
    dataset_id,
    perturbed_target_symbol,
    COUNT(*) AS n_total,
    COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
    COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up
FROM {source_table}
WHERE padj <= 0.05
GROUP BY dataset_id, perturbed_target_symbol""",
            [
                (
                    "idx_perturb_seq_summary_perturbation_pk",
                    "CREATE UNIQUE INDEX {idx} ON {view} (dataset_id, perturbed_target_symbol)",
                ),
            ],
        ),
        (
            "perturb_seq_summary_effect",
            """CREATE MATERIALIZED VIEW {view} AS
SELECT
    dataset_id,
    gene,
    COUNT(*) AS n_total,
    COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
    COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up,
    AVG(score_value) AS avg_score
FROM {source_table}
WHERE padj <= 0.05
GROUP BY dataset_id, gene""",
            [
                (
                    "idx_perturb_seq_summary_effect_pk",
                    "CREATE UNIQUE INDEX {idx} ON {view} (dataset_id, gene)",
                ),
            ],
        ),
        (
            "perturb_seq_summary_dataset",
            """CREATE MATERIALIZED VIEW {view} AS
SELECT
    dataset_id,
    COUNT(*) AS n_total
FROM {source_table}
WHERE gene IS NOT NULL
GROUP BY dataset_id""",
            [
                (
                    "idx_perturb_seq_summary_dataset_pk",
                    "CREATE UNIQUE INDEX {idx} ON {view} (dataset_id)",
                ),
            ],
        ),
    ],
}

# Number of threads for concurrent GCS blob download + Parquet-to-CSV conversion.
GCS_DOWNLOAD_WORKERS = os.cpu_count() or 4


class IngestionMode(enum.Enum):
    LIVE_NONBLOCKING = "live_nonblocking"
    COPY_NONBLOCKING = "copy_nonblocking"
    DIRECT_BLOCKING = "direct_blocking"

    def __str__(self):
        return self.value


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


def get_all_sync_states(cursor) -> Dict[str, Dict[str, Any]]:
    """Gets the current sync state for all tables and datasets from Postgres."""
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS sync_state (
            table_name TEXT NOT NULL,
            dataset_id TEXT NOT NULL,
            last_synced_at TIMESTAMP WITHOUT TIME ZONE,
            PRIMARY KEY (table_name, dataset_id)
        );
    """
    )
    states = {}
    cursor.execute("SELECT table_name, dataset_id, last_synced_at FROM sync_state")
    for table_name, dataset_id, last_synced_at in cursor.fetchall():
        if table_name not in states:
            states[table_name] = {}
        states[table_name][dataset_id] = last_synced_at
    return states


def get_bq_latest_timestamps_and_counts(
    bq_client, bq_dataset, bq_table, bq_location
) -> Dict[str, Tuple[Any, int]]:
    """Gets the latest max_ingested_at and row count for every dataset_id in a BQ table."""
    query = f"""
        SELECT dataset_id, MAX(max_ingested_at) as latest_ts, COUNT(*) as row_count
        FROM `{bq_dataset}.{bq_table}`
        GROUP BY dataset_id
    """
    query_job = bq_client.query(query, location=bq_location)
    results = query_job.result()
    return {row.dataset_id: (row.latest_ts, row.row_count) for row in results}


def export_dataset_to_gcs(
    bq_client, bq_dataset, bq_table, bq_location, gcs_bucket, gcs_prefix, dataset_id
) -> str:
    """Exports a specific dataset from BigQuery to GCS as Parquet via a temp table."""
    destination_uri = f"gs://{gcs_bucket}/{gcs_prefix}-*.parquet"
    dataset_ref = bq_client.dataset(bq_dataset)

    # Query to a temporary table to filter by dataset_id
    temp_table_id = f"temp_sync_{uuid.uuid4().hex}"
    temp_table_ref = dataset_ref.table(temp_table_id)

    query = f"SELECT * FROM `{bq_dataset}.{bq_table}` WHERE dataset_id = '{dataset_id}'"
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
    except (Exception, KeyboardInterrupt) as e:
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


def _download_and_convert_blob(blob, bq_schema):
    """Download a single Parquet shard and convert to a TSV buffer."""
    data = blob.download_as_bytes()
    table = pq.read_table(io.BytesIO(data))
    return _prepare_table_for_copy(table, bq_schema)


def load_parquet_from_gcs_to_pg(cursor, pg_table, gcs_bucket, gcs_prefix, bq_schema):
    """Loads Parquet files from GCS into Postgres using COPY."""
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_prefix))

    if not blobs:
        return

    copy_sql = sql.SQL(
        "COPY {} FROM STDIN WITH (FORMAT CSV, DELIMITER E'\\t', QUOTE '\"', NULL '')"
    ).format(sql.Identifier(pg_table))

    max_workers = GCS_DOWNLOAD_WORKERS
    blob_iter = iter(blobs)
    pbar = tqdm(total=len(blobs), desc="        Loading shards", leave=False)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        pending = {}
        for blob in iter(lambda: next(blob_iter, None), None):
            fut = pool.submit(_download_and_convert_blob, blob, bq_schema)
            pending[fut] = blob
            if len(pending) >= max_workers:
                break

        while pending:
            done, _ = concurrent.futures.wait(
                pending, return_when=concurrent.futures.FIRST_COMPLETED
            )
            for fut in done:
                tsv_buf = fut.result()
                cursor.copy_expert(copy_sql, tsv_buf)
                tsv_buf.close()
                del pending[fut]
                pbar.update(1)

                next_blob = next(blob_iter, None)
                if next_blob is not None:
                    new_fut = pool.submit(
                        _download_and_convert_blob, next_blob, bq_schema
                    )
                    pending[new_fut] = next_blob

    pbar.close()


def cleanup_gcs(gcs_bucket, gcs_prefix):
    """Removes temporary files from GCS."""
    logging.info(f"      - Cleaning up GCS files...")
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_prefix))
    for blob in blobs:
        blob.delete()


def delete_dataset_from_pg(cursor, pg_table, dataset_id):
    """Deletes all rows for a given dataset_id from a Postgres table."""
    cursor.execute(
        sql.SQL("DELETE FROM {} WHERE dataset_id = %s").format(
            sql.Identifier(pg_table)
        ),
        (dataset_id,),
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
    """Ensures the target table exists in Postgres, creating it if necessary."""
    cursor.execute(
        "SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = %s)",
        (pg_table,),
    )
    if not cursor.fetchone()[0]:
        logging.info(f"Table {pg_table} does not exist. Creating.")
        columns = [f"{field.name} {get_pg_type(field)}" for field in bq_schema]
        cursor.execute(
            sql.SQL("CREATE TABLE {} ({})").format(
                sql.Identifier(pg_table),
                sql.SQL(", ").join(map(sql.SQL, columns)),
            )
        )


def copy_table_schema(cursor, src_table, dst_table):
    """Creates a dst_table with the same schema as src_table (no data, no indexes)."""
    logging.info(f"      - Creating {dst_table} (copy of {src_table} schema)...")
    cursor.execute(
        sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(sql.Identifier(dst_table))
    )
    cursor.execute(
        sql.SQL("CREATE TABLE {} (LIKE {})").format(
            sql.Identifier(dst_table), sql.Identifier(src_table)
        )
    )


def drop_indexes(cursor, table_name, suffix=""):
    """Drops all indexes for a given table based on INDEX_DEFINITIONS."""
    if table_name not in INDEX_DEFINITIONS:
        return
    logging.info(f"      - Dropping indexes for {table_name}{suffix}...")
    for index_name, _ in INDEX_DEFINITIONS[table_name]:
        idx = f"{index_name}{suffix}"
        cursor.execute(sql.SQL("DROP INDEX IF EXISTS {}").format(sql.Identifier(idx)))


def create_indexes(cursor, table_name, suffix=""):
    """Creates all indexes for a given table based on INDEX_DEFINITIONS."""
    if table_name not in INDEX_DEFINITIONS:
        return
    logging.info(
        f"      - Creating indexes for {table_name}{suffix} (this may take a while)..."
    )
    for index_name, index_sql_template in INDEX_DEFINITIONS[table_name]:
        idx = f"{index_name}{suffix}"
        logging.info(f"        Creating {idx}...")
        index_sql = index_sql_template.format(
            idx=sql.Identifier(idx).as_string(cursor.connection),
            table=sql.Identifier(f"{table_name}{suffix}").as_string(cursor.connection),
        )
        cursor.execute(index_sql)


def create_materialized_views(cursor, table_name, suffix=""):
    """Creates materialized views (with optional suffix) from definitions."""
    if table_name not in MATERIALIZED_VIEW_DEFINITIONS:
        return
    for view_name, create_sql_template, view_indexes in MATERIALIZED_VIEW_DEFINITIONS[
        table_name
    ]:
        view = f"{view_name}{suffix}"
        source = f"{table_name}{suffix}"
        logging.info(f"      - Creating materialized view {view}...")
        cursor.execute(
            sql.SQL("DROP MATERIALIZED VIEW IF EXISTS {}").format(sql.Identifier(view))
        )
        create_sql = create_sql_template.format(
            view=sql.Identifier(view).as_string(cursor.connection),
            source_table=sql.Identifier(source).as_string(cursor.connection),
        )
        cursor.execute(create_sql)
        for mv_index_name, mv_index_sql_template in view_indexes:
            idx = f"{mv_index_name}{suffix}"
            logging.info(f"        Creating index {idx}...")
            mv_index_sql = mv_index_sql_template.format(
                idx=sql.Identifier(idx).as_string(cursor.connection),
                view=sql.Identifier(view).as_string(cursor.connection),
            )
            cursor.execute(mv_index_sql)


def refresh_materialized_views(cursor, table_name, concurrently=False):
    """Refreshes materialized views associated with the table."""
    if table_name not in MATERIALIZED_VIEW_DEFINITIONS:
        return
    for view_name, _, _ in MATERIALIZED_VIEW_DEFINITIONS[table_name]:
        logging.info(f"      - Refreshing materialized view {view_name}...")
        conc_clause = "CONCURRENTLY " if concurrently else ""
        try:
            cursor.execute(
                sql.SQL("REFRESH MATERIALIZED VIEW {}{}").format(
                    sql.SQL(conc_clause), sql.Identifier(view_name)
                )
            )
        except psycopg2.Error as e:
            logging.warning(
                f"        Failed to refresh {view_name} {conc_clause.strip()}: {e}. Trying without CONCURRENTLY."
            )
            if concurrently:
                cursor.connection.rollback()  # Required to recover from error in transaction if any (though we are usually autocommit here if using concurrent)
                # If we were in a transaction block, we can't retry easily without rollback.
                # Assuming this runs in a state where we can retry.
                cursor.execute(
                    sql.SQL("REFRESH MATERIALIZED VIEW {}").format(
                        sql.Identifier(view_name)
                    )
                )


def swap_table(conn, table_name):
    """Atomically swaps {table_name}_upd into {table_name}."""
    upd = f"{table_name}_upd"
    logging.info(f"      - Swapping {upd} -> {table_name} (atomic transaction)...")

    with conn:
        with conn.cursor() as cursor:
            # 1. Drop original table (CASCADE drops dependent MVs)
            cursor.execute(
                sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(
                    sql.Identifier(table_name)
                )
            )
            # 2. Rename _upd table
            cursor.execute(
                sql.SQL("ALTER TABLE {} RENAME TO {}").format(
                    sql.Identifier(upd), sql.Identifier(table_name)
                )
            )
            # 3. Rename indexes
            if table_name in INDEX_DEFINITIONS:
                for index_name, _ in INDEX_DEFINITIONS[table_name]:
                    cursor.execute(
                        sql.SQL("ALTER INDEX {} RENAME TO {}").format(
                            sql.Identifier(f"{index_name}_upd"),
                            sql.Identifier(index_name),
                        )
                    )
            # 4. Rename MVs and their indexes
            if table_name in MATERIALIZED_VIEW_DEFINITIONS:
                for view_name, _, view_indexes in MATERIALIZED_VIEW_DEFINITIONS[
                    table_name
                ]:
                    cursor.execute(
                        sql.SQL("ALTER MATERIALIZED VIEW {} RENAME TO {}").format(
                            sql.Identifier(f"{view_name}_upd"),
                            sql.Identifier(view_name),
                        )
                    )
                    for mv_index_name, _ in view_indexes:
                        cursor.execute(
                            sql.SQL("ALTER INDEX {} RENAME TO {}").format(
                                sql.Identifier(f"{mv_index_name}_upd"),
                                sql.Identifier(mv_index_name),
                            )
                        )
    logging.info(f"      - Swap complete for {table_name}.")


# ------------------------------------------------------------------------------
# Core Logic Class
# ------------------------------------------------------------------------------


class TableSynchronizer:
    def __init__(
        self,
        mode: IngestionMode,
        pg_conn_str: str,
        bq_dataset: str,
        bq_location: str,
        gcs_bucket: str,
    ):
        self.mode = mode
        self.pg_conn_str = pg_conn_str
        self.bq_dataset = bq_dataset
        self.bq_location = bq_location
        self.gcs_bucket = gcs_bucket
        self.bq_client = bigquery.Client()

    def sync_table(self, table_name: str, plan: Dict[str, Any]):
        logging.info(f"Syncing table {table_name} in mode {self.mode}...")
        try:
            if self.mode == IngestionMode.LIVE_NONBLOCKING:
                self._sync_live_nonblocking(table_name, plan)
            elif self.mode == IngestionMode.COPY_NONBLOCKING:
                self._sync_copy_nonblocking(table_name, plan)
            elif self.mode == IngestionMode.DIRECT_BLOCKING:
                self._sync_direct_blocking(table_name, plan)

            logging.info(f"Successfully synced {table_name}.")

        except Exception as e:
            if isinstance(e, KeyboardInterrupt):
                raise e
            logging.error(f"Error syncing {table_name}: {e}.")
            raise e

    def _get_bq_schema(self, table_name):
        return self.bq_client.get_table(f"{self.bq_dataset}.{table_name}").schema

    def _ingest_dataset_logic(self, cursor, pg_table, table_name, ds_id, bq_schema):
        """Standard ingestion: export, delete, copy."""
        gcs_prefix = f"tmp/{table_name}/{ds_id}/{uuid.uuid4().hex}"
        logging.info(f"    Processing {ds_id}...")
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
            delete_dataset_from_pg(cursor, pg_table, ds_id)
            load_parquet_from_gcs_to_pg(
                cursor,
                pg_table,
                self.gcs_bucket,
                gcs_prefix,
                bq_schema,
            )
        finally:
            cleanup_gcs(self.gcs_bucket, gcs_prefix)

    def _sync_live_nonblocking(self, table_name: str, plan: Dict[str, Any]):
        """
        Mode 1: LIVE_NONBLOCKING
        - Single transaction for all data updates (Delete + Insert).
        - No table copying.
        - No index dropping (updates are slower but non-blocking).
        - Updates sync_state in same transaction.
        - Refreshes MVs afterwards (concurrently if possible).
        """
        bq_schema = self._get_bq_schema(table_name)
        datasets = sorted(plan["to_update"] + plan["to_insert"])

        # 1. Data Update Transaction
        with psycopg2.connect(self.pg_conn_str) as conn:
            with conn.cursor() as cursor:
                ensure_pg_table_exists(cursor, table_name, bq_schema)

                # Perform all updates in one transaction
                logging.info(
                    f"    Starting data update transaction for {len(datasets)} datasets..."
                )
                for ds_id in tqdm(
                    datasets, desc=f"    Syncing {table_name}", unit="dataset"
                ):
                    # Ingest data
                    self._ingest_dataset_logic(
                        cursor, table_name, table_name, ds_id, bq_schema
                    )

                    # Update sync state
                    bq_ts = plan["bq_info"][ds_id][0]
                    update_sync_state(cursor, table_name, ds_id, bq_ts)

        # 2. Materialized View Refresh (outside transaction)
        # We try CONCURRENTLY to be non-blocking.
        with psycopg2.connect(self.pg_conn_str) as conn:
            conn.autocommit = (
                True  # Required for REFRESH MATERIALIZED VIEW CONCURRENTLY
            )
            with conn.cursor() as cursor:
                refresh_materialized_views(cursor, table_name, concurrently=True)

    def _sync_copy_nonblocking(self, table_name: str, plan: Dict[str, Any]):
        """
        Mode 2: COPY_NONBLOCKING (Shadow Table Strategy)
        - Create _upd table (copy).
        - Drop indexes on _upd.
        - Ingest data into _upd.
        - Recreate indexes/MVs on _upd.
        - Update sync_state updates are collected and applied during swap.
        - Atomic SWAP.
        """
        upd_table = f"{table_name}_upd"
        bq_schema = self._get_bq_schema(table_name)
        datasets = sorted(plan["to_update"] + plan["to_insert"])

        sync_state_updates = []  # List of (ds_id, timestamp) to apply at the end

        # Connection for preparation (autocommit to avoid long open tx)
        conn = psycopg2.connect(self.pg_conn_str)
        conn.autocommit = True
        try:
            with conn.cursor() as cursor:
                ensure_pg_table_exists(cursor, table_name, bq_schema)

                # A. Prepare Shadow Table
                copy_table_schema(cursor, table_name, upd_table)
                # Copy data *except* the datasets we are about to update?
                # The prompt says "delete old rows and add new rows".
                # Efficient approach: Copy ALL data, then delete specific datasets?
                # Or Copy filtering out datasets? Filtering in SQL is better.
                # However, original logic was: Copy everything, then DELETE rows for datasets being updated.
                # Let's stick to that for simplicity and correctness.
                logging.info(
                    f"      - Copying data from {table_name} to {upd_table}..."
                )
                cursor.execute(
                    sql.SQL("INSERT INTO {} SELECT * FROM {}").format(
                        sql.Identifier(upd_table), sql.Identifier(table_name)
                    )
                )

                # B. Drop indexes on _upd (should be empty, but just in case)
                drop_indexes(cursor, table_name, suffix="_upd")

                # C. Ingest Loop (Non-transactional on shadow table)
                for ds_id in tqdm(
                    datasets, desc=f"    Syncing {upd_table}", unit="dataset"
                ):
                    self._ingest_dataset_logic(
                        cursor, upd_table, table_name, ds_id, bq_schema
                    )
                    # Collect sync state update
                    bq_ts = plan["bq_info"][ds_id][0]
                    sync_state_updates.append((ds_id, bq_ts))

                # D. Recreate Indexes on _upd
                create_indexes(cursor, table_name, suffix="_upd")

                # E. Create MVs on _upd
                create_materialized_views(cursor, table_name, suffix="_upd")

        except Exception:
            # Cleanup on failure
            with psycopg2.connect(self.pg_conn_str) as clean_conn:
                clean_conn.autocommit = True
                with clean_conn.cursor() as clean_cur:
                    clean_cur.execute(
                        sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(
                            sql.Identifier(upd_table)
                        )
                    )
            raise
        finally:
            conn.close()

        # F. Atomic Swap & Sync State Update
        with psycopg2.connect(self.pg_conn_str) as swap_conn:
            try:
                with swap_conn.cursor() as cursor:
                    # We need to update sync_state inside this transaction
                    logging.info("      - Updating sync_state entries...")
                    for ds_id, bq_ts in sync_state_updates:
                        update_sync_state(cursor, table_name, ds_id, bq_ts)

                # Perform the swap (logic handles internal transaction management,
                # but swap_table expects to manage its own transaction block context usually
                # or be part of one. The `swap_table` function does `with conn:`, which starts a transaction.
                # So we should pass the connection.
                swap_table(swap_conn, table_name)
            except Exception:
                # If swap fails, clean up
                with psycopg2.connect(self.pg_conn_str) as clean_conn:
                    clean_conn.autocommit = True
                    with clean_conn.cursor() as clean_cur:
                        clean_cur.execute(
                            sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(
                                sql.Identifier(upd_table)
                            )
                        )
                raise

    def _sync_direct_blocking(self, table_name: str, plan: Dict[str, Any]):
        """
        Mode 3: DIRECT_BLOCKING
        - Huge single transaction.
        - Drop indexes.
        - Delete old -> Ingest new.
        - Recreate indexes.
        - Refresh MVs.
        - Update sync_state.
        """
        bq_schema = self._get_bq_schema(table_name)
        datasets = sorted(plan["to_update"] + plan["to_insert"])

        with psycopg2.connect(self.pg_conn_str) as conn:
            with conn.cursor() as cursor:
                ensure_pg_table_exists(cursor, table_name, bq_schema)

                logging.info(f"    Starting blocking transaction for {table_name}...")

                # 1. Drop Indexes
                drop_indexes(cursor, table_name)

                # 2. Ingest Data
                for ds_id in tqdm(
                    datasets, desc=f"    Syncing {table_name}", unit="dataset"
                ):
                    self._ingest_dataset_logic(
                        cursor, table_name, table_name, ds_id, bq_schema
                    )

                    # Update sync state
                    bq_ts = plan["bq_info"][ds_id][0]
                    update_sync_state(cursor, table_name, ds_id, bq_ts)

                # 3. Recreate Indexes
                create_indexes(cursor, table_name)

                # 4. Refresh Materialized Views (Standard sync refresh inside transaction)
                if table_name in MATERIALIZED_VIEW_DEFINITIONS:
                    for view_name, _, _ in MATERIALIZED_VIEW_DEFINITIONS[table_name]:
                        logging.info(
                            f"      - Refreshing materialized view {view_name} (blocking)..."
                        )
                        cursor.execute(
                            sql.SQL("REFRESH MATERIALIZED VIEW {}").format(
                                sql.Identifier(view_name)
                            )
                        )


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bq-dataset", required=True)
    parser.add_argument("--bq-location", required=True)
    parser.add_argument("--pg-conn", required=True)
    parser.add_argument("--gcs-bucket", required=True)
    parser.add_argument(
        "--ingestion-mode",
        type=IngestionMode,
        choices=list(IngestionMode),
        default=IngestionMode.COPY_NONBLOCKING,
        help="Choose operation mode: live_nonblocking (slow, no copy, no index drop), copy_nonblocking (default, copy table, atomic swap), direct_blocking (fastest, blocks table, drops indexes).",
    )
    parser.add_argument(
        "--yes", action="store_true", help="Proceed without confirmation"
    )
    args = parser.parse_args()

    mode = args.ingestion_mode
    # Since argparse with type=Enum returns the Enum member, we can use it directly.
    # However, depending on python version it might be different. Let's ensure compatibility.
    # If type=IngestionMode is used, args.ingestion_mode will be an IngestionMode member.

    logging.info(f"Starting sync in mode: {mode}")

    bq_client = bigquery.Client()
    synchronizer = TableSynchronizer(
        mode, args.pg_conn, args.bq_dataset, args.bq_location, args.gcs_bucket
    )

    try:
        # Phase 1: Planning
        with psycopg2.connect(args.pg_conn) as plan_conn:
            with plan_conn.cursor() as cursor:
                logging.info("Fetching current sync states from Postgres...")
                all_states = get_all_sync_states(cursor)

                sync_plan = {}
                total_datasets = 0

                for table_name in TABLES_TO_SYNC:
                    logging.info(f"Checking BQ table {table_name}...")
                    bq_info = get_bq_latest_timestamps_and_counts(
                        bq_client, args.bq_dataset, table_name, args.bq_location
                    )
                    current_table_states = all_states.get(table_name, {})

                    to_update = []
                    to_insert = []

                    for ds_id, (bq_ts, row_count) in bq_info.items():
                        if ds_id in current_table_states:
                            pg_ts = current_table_states[ds_id]
                            # Simple timezone-aware comparison
                            if pg_ts is None or bq_ts > pg_ts.replace(
                                tzinfo=bq_ts.tzinfo
                            ):
                                to_update.append(ds_id)
                        else:
                            to_insert.append(ds_id)

                    if to_update or to_insert:
                        sync_plan[table_name] = {
                            "to_update": to_update,
                            "to_insert": to_insert,
                            "bq_info": bq_info,
                        }
                        total_datasets += len(to_update) + len(to_insert)

        # Summary and Confirmation
        if not sync_plan:
            logging.info("Everything is up to date. Nothing to sync.")
            return

        print("\n--- Sync Summary ---")
        for table_name, plan in sync_plan.items():
            print(f"Table: {table_name}")
            if plan["to_update"]:
                print(f"  Datasets to update: {len(plan['to_update'])}")
                print(
                    f"    {', '.join(plan['to_update'][:5])}{'...' if len(plan['to_update']) > 5 else ''}"
                )
            if plan["to_insert"]:
                print(f"  Datasets to insert: {len(plan['to_insert'])}")
                print(
                    f"    {', '.join(plan['to_insert'][:5])}{'...' if len(plan['to_insert']) > 5 else ''}"
                )

        print(f"\nTotal datasets to process: {total_datasets}")
        print(f"Ingestion Mode: {mode}")
        print("--------------------\n")

        if not args.yes:
            try:
                confirm = input("Proceed with sync? (y/N): ")
            except EOFError:
                confirm = "n"
            if confirm.lower() != "y":
                logging.info("Sync cancelled by user.")
                sys.exit(0)

        # Phase 2: Execution
        for table_name, plan in sync_plan.items():
            synchronizer.sync_table(table_name, plan)

        logging.info("All sync operations completed.")

    except psycopg2.Error as e:
        logging.error(f"Postgres connection error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("\nExiting...")
        sys.exit(0)


if __name__ == "__main__":
    main()
