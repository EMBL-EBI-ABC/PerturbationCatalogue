import argparse
import concurrent.futures
import contextlib
import enum
import io
import logging
import os
import uuid
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
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

# Number of threads for concurrent GCS blob download + Parquet-to-CSV conversion.
GCS_DOWNLOAD_WORKERS = os.cpu_count() or 4


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


def _download_and_convert_blob(blob, bq_schema):
    """Download a single Parquet shard and convert to a TSV buffer."""
    data = blob.download_as_bytes()
    table = pq.read_table(io.BytesIO(data))
    return _prepare_table_for_copy(table, bq_schema)


def load_parquet_from_gcs_to_pg(
    cursor, pg_table, gcs_bucket, gcs_prefix, bq_schema, gcs_client
):
    """Loads Parquet files from GCS into Postgres using COPY."""
    bucket = gcs_client.get_bucket(gcs_bucket)
    logging.info(f"        Listing blobs in {gcs_prefix}...")
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


def cleanup_gcs(gcs_bucket, gcs_prefix, gcs_client):
    """Removes temporary files from GCS."""
    logging.info(f"      - Cleaning up GCS files...")
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_prefix))
    for blob in blobs:
        blob.delete()


def delete_dataset_from_pg(cursor, pg_table, dataset_id):
    """Deletes all rows for a given dataset_id from a Postgres table."""
    logging.info(f"        Deleting {dataset_id} from {pg_table}...")
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


def drop_indexes(cursor, table_name):
    """Drops all indexes for a given table based on INDEX_DEFINITIONS."""
    if table_name not in INDEX_DEFINITIONS:
        return
    logging.info(f"    Dropping indexes for {table_name}...")
    for index_name, _ in INDEX_DEFINITIONS[table_name]:
        idx = index_name
        logging.info(f"        Dropping {idx}...")
        cursor.execute(sql.SQL("DROP INDEX IF EXISTS {}").format(sql.Identifier(idx)))


def create_indexes(cursor, table_name):
    """Creates all indexes for a given table based on INDEX_DEFINITIONS."""
    if table_name not in INDEX_DEFINITIONS:
        return
    logging.info(
        f"      - Creating indexes for {table_name} (this may take a while)..."
    )
    for index_name, index_sql_template in INDEX_DEFINITIONS[table_name]:
        idx = index_name
        logging.info(f"        Creating {idx}...")
        index_sql = index_sql_template.format(
            idx=sql.Identifier(idx).as_string(cursor.connection),
            table=sql.Identifier(table_name).as_string(cursor.connection),
        )
        cursor.execute(index_sql)


def ensure_perturb_seq_summary_views(cursor):
    """Create ENSG summary materialized views and unique indexes when missing."""
    cursor.execute(
        sql.SQL(
            """
            CREATE MATERIALIZED VIEW IF NOT EXISTS perturb_seq_summary_perturbation AS
            SELECT
                dataset_id,
                perturbed_target_ensg,
                COUNT(*) AS n_total,
                COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
                COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up
            FROM perturb_seq_dea
            WHERE padj <= 0.05 AND perturbed_target_ensg IS NOT NULL
            GROUP BY dataset_id, perturbed_target_ensg
            """
        )
    )
    cursor.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS idx_perturb_seq_summary_perturbation_pk
        ON perturb_seq_summary_perturbation (dataset_id, perturbed_target_ensg)
        """
    )

    cursor.execute(
        sql.SQL(
            """
            CREATE MATERIALIZED VIEW IF NOT EXISTS perturb_seq_summary_effect AS
            SELECT
                dataset_id,
                effect_gene_ensg,
                COUNT(*) AS n_total,
                COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
                COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up,
                AVG(score_value) AS avg_score
            FROM perturb_seq_dea
            WHERE padj <= 0.05 AND effect_gene_ensg IS NOT NULL
            GROUP BY dataset_id, effect_gene_ensg
            """
        )
    )
    cursor.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS idx_perturb_seq_summary_effect_pk
        ON perturb_seq_summary_effect (dataset_id, effect_gene_ensg)
        """
    )

    cursor.execute(
        sql.SQL(
            """
            CREATE MATERIALIZED VIEW IF NOT EXISTS perturb_seq_summary_dataset AS
            SELECT dataset_id, COUNT(*) AS n_total
            FROM perturb_seq_dea
            WHERE effect_gene_ensg IS NOT NULL
            GROUP BY dataset_id
            """
        )
    )
    cursor.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS idx_perturb_seq_summary_dataset_pk
        ON perturb_seq_summary_dataset (dataset_id)
        """
    )


def refresh_materialized_views(cursor, table_name, concurrently=False):
    """Refreshes materialized views associated with the table."""
    if table_name not in TABLE_MATERIALIZED_VIEWS:
        return
    if table_name == "perturb_seq_dea":
        ensure_perturb_seq_summary_views(cursor)
    for view_name in TABLE_MATERIALIZED_VIEWS[table_name]:
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
                # If we were in a transaction block, we can't retry easily without rollback.
                # Assuming this runs in a state where we can retry (autocommit=True for MV refresh).
                cursor.execute(
                    sql.SQL("REFRESH MATERIALIZED VIEW {}").format(
                        sql.Identifier(view_name)
                    )
                )


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
                self.gcs_client,
            )
        finally:
            cleanup_gcs(self.gcs_bucket, gcs_prefix, self.gcs_client)

    def _sync_unified(self, table_name: str, plan: Dict[str, Any]):
        """
        Unified Sync Mode:
        - Single transaction for:
            1. (Optional) Drop Indexes
            2. Data Updates (Delete + Insert) for all datasets
            3. (Optional) Recreate Indexes
            4. Update sync_state
        - (Separate) Refresh MVs concurrently.
        """
        bq_schema = self._get_bq_schema(table_name)
        datasets = sorted(plan["to_update"] + plan["to_insert"])

        # 1. Main Transaction Block
        with psycopg2.connect(self.pg_conn_str) as conn:
            with conn.cursor() as cursor:
                ensure_pg_table_exists(cursor, table_name, bq_schema)

                logging.info(
                    f"    Starting transaction for {len(datasets)} datasets..."
                )

                # A. Drop Indexes (if requested)
                if self.drop_and_recreate_indexes:
                    drop_indexes(cursor, table_name)

                # B. Ingest Loop
                for ds_id in tqdm(
                    datasets, desc=f"    Syncing {table_name}", unit="dataset"
                ):
                    # Data Ingestion
                    self._ingest_dataset_logic(
                        cursor, table_name, table_name, ds_id, bq_schema
                    )

                    # Update sync state
                    bq_ts = plan["bq_info"][ds_id][0]
                    update_sync_state(cursor, table_name, ds_id, bq_ts)

                # C. Recreate Indexes (if requested)
                if self.drop_and_recreate_indexes:
                    create_indexes(cursor, table_name)

        # 2. Materialized View Refresh (concurrently, separate connection)
        # Only if we successfully committed the transaction above.
        with psycopg2.connect(self.pg_conn_str) as conn:
            conn.autocommit = (
                True  # Required for REFRESH MATERIALIZED VIEW CONCURRENTLY
            )
            with conn.cursor() as cursor:
                refresh_materialized_views(cursor, table_name, concurrently=True)


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

                # 2. Compare with PG state
                pg_info = pg_states.get(table_name, {})
                to_insert = []
                to_update = []

                if table_name in force_pg_tables:
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
