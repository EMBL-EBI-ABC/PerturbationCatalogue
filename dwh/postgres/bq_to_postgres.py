import argparse
import concurrent.futures
import logging
import io
import os
import sys
import uuid
from concurrent.futures import ThreadPoolExecutor
from google.cloud import bigquery, storage
import psycopg2
from psycopg2 import sql
from tqdm import tqdm
import pyarrow as pa
import pyarrow.csv as pa_csv
import pyarrow.parquet as pq


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

TABLES_TO_SYNC = [
    "crispr_data",
    "mave_data",
    "perturb_seq_dea",
    "perturb_seq_gsea",
]

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

# Materialized view definitions, keyed by base table name.
# Each entry is (view_name, create_sql_template, [(index_name, index_sql_template), ...]).
# Templates use {view} for the view name and {source_table} for the base table.
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


def get_all_sync_states(cursor):
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


def get_bq_latest_timestamps_and_counts(bq_client, bq_dataset, bq_table, bq_location):
    """Gets the latest max_ingested_at and row count for every dataset_id in a BQ table."""
    query = f"""
        SELECT dataset_id, MAX(max_ingested_at) as latest_ts, COUNT(*) as row_count
        FROM `{bq_dataset}.{bq_table}`
        GROUP BY dataset_id
    """
    query_job = bq_client.query(query, location=bq_location)
    results = query_job.result()
    return {row.dataset_id: (row.latest_ts, row.row_count) for row in results}


# Number of threads for concurrent GCS blob download + Parquet-to-CSV conversion.
# Auto-detected from available CPUs; each worker holds one shard in memory.
GCS_DOWNLOAD_WORKERS = os.cpu_count() or 4


def export_dataset_to_gcs(
    bq_client, bq_dataset, bq_table, bq_location, gcs_bucket, gcs_prefix, dataset_id
):
    """Exports a specific dataset from BigQuery to GCS as Parquet."""
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
    """Convert a PyArrow list-typed column to Postgres array literal strings.

    E.g. ["a", "b"] -> '{"a","b"}'
    """
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
    """Prepare a PyArrow table for Postgres COPY: handle array columns and
    convert to a CSV (tab-delimited) bytes buffer using PyArrow's C++ writer.
    """
    # Build a set of array column names for targeted conversion
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

    # Write to TSV bytes buffer using PyArrow's C++ CSV writer (very fast)
    buf = io.BytesIO()
    write_options = pa_csv.WriteOptions(
        include_header=False,
        delimiter="\t",
    )
    pa_csv.write_csv(table, buf, write_options=write_options)
    buf.seek(0)
    return buf


def _download_and_convert_blob(blob, bq_schema):
    """Download a single Parquet shard from GCS and convert to a TSV buffer.

    This runs in a worker thread to overlap GCS I/O with CPU work.
    Returns a BytesIO buffer ready for COPY FROM STDIN.
    """
    data = blob.download_as_bytes()
    table = pq.read_table(io.BytesIO(data))
    return _prepare_table_for_copy(table, bq_schema)


def load_parquet_from_gcs_to_pg(cursor, pg_table, gcs_bucket, gcs_prefix, bq_schema):
    """Loads Parquet files from GCS into Postgres using COPY.

    Uses concurrent GCS downloads and PyArrow's native C++ CSV writer for
    maximum throughput. Blobs are downloaded and converted to TSV in parallel
    threads; the main thread feeds buffers to Postgres sequentially.

    Memory is bounded: at most ``GCS_DOWNLOAD_WORKERS`` shard buffers exist in
    memory at any time. New work is only submitted as previous results are
    consumed by Postgres COPY.
    """
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_prefix))

    copy_sql = sql.SQL(
        "COPY {} FROM STDIN WITH (FORMAT CSV, DELIMITER E'\\t', QUOTE '\"', NULL '')"
    ).format(sql.Identifier(pg_table))

    max_workers = GCS_DOWNLOAD_WORKERS
    blob_iter = iter(blobs)
    pbar = tqdm(total=len(blobs), desc="        Loading shards", leave=False)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        # Seed the pool with an initial batch (bounded to max_workers)
        pending = {}
        for blob in iter(lambda: next(blob_iter, None), None):
            fut = pool.submit(_download_and_convert_blob, blob, bq_schema)
            pending[fut] = blob
            if len(pending) >= max_workers:
                break

        # Sliding window: consume one result, submit one new blob
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

                # Submit next blob if available
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


def copy_table(cursor, src_table, dst_table):
    """Creates a copy of a table with all its data.

    The new table has the same column definitions but no indexes or constraints.
    """
    logging.info(f"      - Copying {src_table} -> {dst_table}...")
    cursor.execute(
        sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(sql.Identifier(dst_table))
    )
    cursor.execute(
        sql.SQL("CREATE TABLE {} (LIKE {})").format(
            sql.Identifier(dst_table), sql.Identifier(src_table)
        )
    )
    cursor.execute(
        sql.SQL("INSERT INTO {} SELECT * FROM {}").format(
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
        logging.info(f"        Dropping {idx}...")
        cursor.execute(sql.SQL("DROP INDEX IF EXISTS {}").format(sql.Identifier(idx)))


def create_indexes(cursor, table_name, suffix=""):
    """Creates all indexes for a given table based on INDEX_DEFINITIONS.

    When suffix is provided, creates indexes with the suffix in both the index
    name and the target table name.
    """
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
    """Creates materialized views (with optional suffix) from definitions.

    Views are created against {table_name}{suffix} and named {view_name}{suffix}.
    """
    if table_name not in MATERIALIZED_VIEW_DEFINITIONS:
        return
    for view_name, create_sql_template, view_indexes in MATERIALIZED_VIEW_DEFINITIONS[
        table_name
    ]:
        view = f"{view_name}{suffix}"
        source = f"{table_name}{suffix}"
        logging.info(f"      - Creating materialized view {view}...")
        # Drop if exists (e.g. from a previous failed run)
        cursor.execute(
            sql.SQL("DROP MATERIALIZED VIEW IF EXISTS {}").format(sql.Identifier(view))
        )
        create_sql = create_sql_template.format(
            view=sql.Identifier(view).as_string(cursor.connection),
            source_table=sql.Identifier(source).as_string(cursor.connection),
        )
        cursor.execute(create_sql)
        # Create indexes on the materialized view
        for mv_index_name, mv_index_sql_template in view_indexes:
            idx = f"{mv_index_name}{suffix}"
            logging.info(f"        Creating index {idx}...")
            mv_index_sql = mv_index_sql_template.format(
                idx=sql.Identifier(idx).as_string(cursor.connection),
                view=sql.Identifier(view).as_string(cursor.connection),
            )
            cursor.execute(mv_index_sql)


def swap_table(conn, table_name):
    """Atomically swaps {table_name}_upd into {table_name} in a single short transaction.

    Steps (all inside one transaction):
    1. DROP the original table CASCADE (also drops its materialized views)
    2. RENAME the _upd table to the original name
    3. RENAME each _upd materialized view to the original name
    4. RENAME each _upd index to the original name
    """
    upd = f"{table_name}_upd"
    logging.info(f"      - Swapping {upd} -> {table_name} (short transaction)...")

    with conn:
        with conn.cursor() as cursor:
            # 1. Drop original table (CASCADE drops dependent materialized views)
            cursor.execute(
                sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(
                    sql.Identifier(table_name)
                )
            )

            # 2. Rename _upd table to original
            cursor.execute(
                sql.SQL("ALTER TABLE {} RENAME TO {}").format(
                    sql.Identifier(upd), sql.Identifier(table_name)
                )
            )

            # 3. Rename indexes on the table
            if table_name in INDEX_DEFINITIONS:
                for index_name, _ in INDEX_DEFINITIONS[table_name]:
                    cursor.execute(
                        sql.SQL("ALTER INDEX {} RENAME TO {}").format(
                            sql.Identifier(f"{index_name}_upd"),
                            sql.Identifier(index_name),
                        )
                    )

            # 4. Rename materialized views and their indexes
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bq-dataset", required=True)
    parser.add_argument("--bq-location", required=True)
    parser.add_argument("--pg-conn", required=True)
    parser.add_argument("--gcs-bucket", required=True)
    parser.add_argument(
        "--yes", action="store_true", help="Proceed without confirmation"
    )
    args = parser.parse_args()

    bq_client = bigquery.Client()

    try:
        # Phase 1: Planning (Read-only / Meta-data fetch)
        with psycopg2.connect(args.pg_conn) as plan_conn:
            with plan_conn.cursor() as cursor:
                # 1. Fetch current sync states for all tables
                logging.info("Fetching current sync states from Postgres...")
                all_states = get_all_sync_states(cursor)

                # 2. Plan sync for each table
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

        # 3. Summary and confirmation
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
        print("--------------------\n")

        if not args.yes:
            try:
                confirm = input("Proceed with sync? (y/N): ")
            except EOFError:
                confirm = "n"
            if confirm.lower() != "y":
                logging.info("Sync cancelled by user.")
                sys.exit(0)

        # 4. Execution (Per-table: non-transactional work on _upd copy, then atomic swap)
        overall_pbar = tqdm(
            total=total_datasets, desc="Overall Progress", unit="dataset"
        )

        for table_name, plan in sync_plan.items():
            logging.info(f"Syncing table {table_name}...")
            upd_table = f"{table_name}_upd"

            try:
                # Phase A: Non-transactional work on the _upd copy.
                # We use autocommit so each statement commits immediately.
                # This avoids holding a long transaction lock on the original table.
                conn = psycopg2.connect(args.pg_conn)
                conn.autocommit = True
                try:
                    with conn.cursor() as cursor:
                        bq_table_obj = bq_client.get_table(
                            f"{args.bq_dataset}.{table_name}"
                        )
                        ensure_pg_table_exists(cursor, table_name, bq_table_obj.schema)

                        # A1. Copy original table -> _upd table
                        copy_table(cursor, table_name, upd_table)

                        # A2. Drop indexes on _upd table (they were not copied, but
                        #     drop just in case from a previous partial run)
                        drop_indexes(cursor, table_name, suffix="_upd")

                        # A3. For each dataset: delete old rows, load new data
                        all_to_process = sorted(plan["to_update"] + plan["to_insert"])

                        for ds_id in all_to_process:
                            bq_timestamp = plan["bq_info"][ds_id][0]

                            # Delete old data if this is an update
                            if ds_id in plan["to_update"]:
                                delete_dataset_from_pg(cursor, upd_table, ds_id)

                            # Export from BQ and load into _upd table
                            gcs_prefix = f"tmp/{table_name}/{ds_id}/{uuid.uuid4().hex}"
                            logging.info(f"    Processing {ds_id}:")

                            try:
                                logging.info(
                                    f"      - Exporting from BigQuery to GCS..."
                                )
                                export_dataset_to_gcs(
                                    bq_client,
                                    args.bq_dataset,
                                    table_name,
                                    args.bq_location,
                                    args.gcs_bucket,
                                    gcs_prefix,
                                    ds_id,
                                )

                                logging.info(f"      - Loading from GCS to Postgres...")
                                load_parquet_from_gcs_to_pg(
                                    cursor,
                                    upd_table,
                                    args.gcs_bucket,
                                    gcs_prefix,
                                    bq_table_obj.schema,
                                )
                            finally:
                                cleanup_gcs(args.gcs_bucket, gcs_prefix)

                            # Update sync state (uses the original table name as key)
                            update_sync_state(cursor, table_name, ds_id, bq_timestamp)

                            overall_pbar.update(1)

                        # A4. Recreate indexes on _upd table
                        create_indexes(cursor, table_name, suffix="_upd")

                        # A5. Create materialized views from _upd table
                        create_materialized_views(cursor, table_name, suffix="_upd")

                finally:
                    conn.close()

                # Phase B: Atomic swap (short transaction).
                conn = psycopg2.connect(args.pg_conn)
                try:
                    swap_table(conn, table_name)
                finally:
                    conn.close()

                logging.info(f"Successfully synced {table_name}.")

            except Exception as e:
                if isinstance(e, KeyboardInterrupt):
                    raise e  # Allow Ctrl+C to stop everything
                logging.error(f"Error syncing {table_name}: {e}. Skipping this table.")
                # Clean up _upd table if it exists
                try:
                    cleanup_conn = psycopg2.connect(args.pg_conn)
                    cleanup_conn.autocommit = True
                    with cleanup_conn.cursor() as cleanup_cursor:
                        cleanup_cursor.execute(
                            sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(
                                sql.Identifier(upd_table)
                            )
                        )
                    cleanup_conn.close()
                except Exception:
                    pass

        overall_pbar.close()
        logging.info("Multi-table sync operations completed.")

    except psycopg2.Error as e:
        logging.error(f"Postgres connection error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("\nExiting...")
        sys.exit(0)


if __name__ == "__main__":
    main()
