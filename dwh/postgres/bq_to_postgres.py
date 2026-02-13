import argparse
import logging
import io
import sys
import uuid
from google.cloud import bigquery, storage
import psycopg2
from psycopg2 import sql
from tqdm import tqdm
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
            "CREATE INDEX idx_perturbation_dea ON public.perturb_seq_dea (perturbed_target_symbol, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_phenotype_dea",
            "CREATE INDEX idx_phenotype_dea ON public.perturb_seq_dea (gene, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_perturbation_phenotype_dea",
            "CREATE INDEX idx_perturbation_phenotype_dea ON public.perturb_seq_dea (perturbed_target_symbol, gene, dataset_id, padj, score_value, log2foldchange)",
        ),
        (
            "idx_perturb_seq_dea_dataset_id_padj",
            "CREATE INDEX idx_perturb_seq_dea_dataset_id_padj ON public.perturb_seq_dea (dataset_id, padj) WHERE gene IS NOT NULL",
        ),
    ],
    "perturb_seq_gsea": [
        (
            "idx_perturbation_gsea",
            "CREATE INDEX idx_perturbation_gsea ON public.perturb_seq_gsea (perturbed_target_symbol, dataset_id, fdr, nes)",
        ),
    ],
    "crispr_data": [
        (
            "idx_crispr_data_dataset",
            "CREATE INDEX idx_crispr_data_dataset ON public.crispr_data (dataset_id)",
        ),
        (
            "idx_crispr_data_target",
            "CREATE INDEX idx_crispr_data_target ON public.crispr_data (perturbed_target_symbol)",
        ),
    ],
    "mave_data": [
        (
            "idx_mave_data_dataset",
            "CREATE INDEX idx_mave_data_dataset ON public.mave_data (dataset_id)",
        ),
        (
            "idx_mave_data_target",
            "CREATE INDEX idx_mave_data_target ON public.mave_data (perturbed_target_symbol, dataset_id)",
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


def format_row(row, schema):
    """Formats a row for Postgres COPY command (TSV)."""
    formatted = []
    for i, field in enumerate(schema):
        val = row[i]
        if val is None:
            formatted.append("")
        elif field.mode == "REPEATED":
            # Convert list to Postgres array literal: {"val1", "val2"}
            inner = [str(x) if not isinstance(x, str) else x for x in val]
            escaped = []
            for x in inner:
                x_str = str(x).replace('"', '\\"')
                escaped.append(f'"{x_str}"')
            formatted.append(f'{{{",".join(escaped)}}}')
        elif field.field_type == "BOOLEAN":
            formatted.append("true" if val else "false")
        else:
            formatted.append(str(val))
    return "\t".join(formatted)


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


def load_parquet_from_gcs_to_pg(cursor, pg_table, gcs_bucket, gcs_prefix, bq_schema):
    """Loads Parquet files from GCS into Postgres using COPY."""
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_prefix))

    for blob in tqdm(blobs, desc="        Loading shards", leave=False):
        with blob.open("rb") as f:
            parquet_file = pq.ParquetFile(f)
            for i in range(parquet_file.num_row_groups):
                table = parquet_file.read_row_group(i)
                rows = table.to_pylist()

                tsv_data = io.StringIO()
                for row_dict in rows:
                    row_val = [row_dict.get(field.name) for field in bq_schema]
                    tsv_data.write(format_row(row_val, bq_schema) + "\n")

                tsv_data.seek(0)
                cursor.copy_expert(
                    sql.SQL("COPY {} FROM STDIN WITH (FORMAT TEXT, NULL '')").format(
                        sql.Identifier(pg_table)
                    ),
                    tsv_data,
                )


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


def drop_indexes(cursor, table_name):
    """Drops all indexes for a given table based on INDEX_DEFINITIONS."""
    if table_name not in INDEX_DEFINITIONS:
        return
    logging.info(f"      - Dropping indexes for {table_name}...")
    for index_name, _ in INDEX_DEFINITIONS[table_name]:
        cursor.execute(
            sql.SQL("DROP INDEX IF EXISTS {}").format(sql.Identifier(index_name))
        )


def create_indexes(cursor, table_name):
    """Creates all indexes for a given table based on INDEX_DEFINITIONS."""
    if table_name not in INDEX_DEFINITIONS:
        return
    logging.info(
        f"      - Reinstating indexes for {table_name} (this may take a while)..."
    )
    for _, index_sql in INDEX_DEFINITIONS[table_name]:
        cursor.execute(index_sql)


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
        with psycopg2.connect(args.pg_conn) as conn:
            with conn.cursor() as cursor:
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

                # 4. Sequential Execution (one dataset at a time)
                overall_pbar = tqdm(
                    total=total_datasets, desc="Overall Progress", unit="dataset"
                )

                for table_name, plan in sync_plan.items():
                    logging.info(f"Syncing table {table_name}...")

                    bq_table_obj = bq_client.get_table(
                        f"{args.bq_dataset}.{table_name}"
                    )
                    ensure_pg_table_exists(cursor, table_name, bq_table_obj.schema)

                    # Drop indexes before bulk update
                    drop_indexes(cursor, table_name)

                    all_to_process = sorted(plan["to_update"] + plan["to_insert"])

                    for ds_id in all_to_process:
                        row_count = plan["bq_info"][ds_id][1]
                        bq_timestamp = plan["bq_info"][ds_id][0]

                        try:
                            # Delete if update
                            if ds_id in plan["to_update"]:
                                delete_dataset_from_pg(cursor, table_name, ds_id)

                            # Parquet-based sync via GCS
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
                                    table_name,
                                    args.gcs_bucket,
                                    gcs_prefix,
                                    bq_table_obj.schema,
                                )
                            finally:
                                cleanup_gcs(args.gcs_bucket, gcs_prefix)

                            # Update sync states
                            update_sync_state(cursor, table_name, ds_id, bq_timestamp)

                            overall_pbar.update(1)

                        except (Exception, KeyboardInterrupt) as e:
                            conn.rollback()
                            overall_pbar.close()
                            if isinstance(e, KeyboardInterrupt):
                                logging.error(
                                    f"\nSync of {table_name}/{ds_id} interrupted. Rolling back..."
                                )
                            else:
                                logging.error(
                                    f"\nError syncing {table_name}/{ds_id}: {e}. Rolling back..."
                                )
                            sys.exit(1)

                    # Rebuild indexes after all data for this table is loaded
                    create_indexes(cursor, table_name)

                    # COMMIT ONCE PER TABLE
                    conn.commit()

                overall_pbar.close()
                logging.info("Multi-table sync completed successfully.")

    except psycopg2.Error as e:
        logging.error(f"Postgres connection error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("\nExiting...")
        sys.exit(0)


if __name__ == "__main__":
    main()
