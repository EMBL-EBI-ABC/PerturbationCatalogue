import argparse
import logging
import uuid
import io
import sys

from google.cloud import bigquery, storage
import psycopg2
from psycopg2 import sql
from tqdm import tqdm
import pyarrow.parquet as pq

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


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


def get_all_sync_states(pg_conn):
    """Gets the current sync state for all tables and datasets."""
    states = {}
    with psycopg2.connect(pg_conn) as conn:
        with conn.cursor() as cursor:
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
            cursor.execute(
                "SELECT table_name, dataset_id, last_synced_at FROM sync_state"
            )
            for table_name, dataset_id, last_synced_at in cursor.fetchall():
                if table_name not in states:
                    states[table_name] = {}
                states[table_name][dataset_id] = last_synced_at
    return states


def get_bq_latest_timestamps(bq_client, bq_dataset, bq_table, bq_location):
    """Gets the latest max_ingested_at for every dataset_id in a BQ table."""
    query = f"""
        SELECT dataset_id, MAX(max_ingested_at) as latest_ts
        FROM `{bq_dataset}.{bq_table}`
        GROUP BY dataset_id
    """
    query_job = bq_client.query(query, location=bq_location)
    results = query_job.result()
    return {row.dataset_id: row.latest_ts for row in results}


def export_bq_to_gcs(
    bq_client,
    bq_dataset,
    bq_table,
    bq_location,
    gcs_bucket,
    gcs_file_path_prefix,
    dataset_ids,
):
    """Exports specific datasets from BigQuery to GCS in Parquet format."""
    dataset_ref = bq_client.dataset(bq_dataset)
    destination_uri = f"gs://{gcs_bucket}/{gcs_file_path_prefix}-*.parquet"

    job_config = bigquery.ExtractJobConfig(destination_format="PARQUET")

    # Query to a temporary table, then export.
    temp_table_id = f"temp_export_{uuid.uuid4().hex}"
    temp_table_ref = dataset_ref.table(temp_table_id)

    # Convert dataset_ids to a SQL-compatible list
    dataset_ids_str = ", ".join([f"'{ds}'" for ds in dataset_ids])
    query = f"""
    SELECT *
    FROM `{bq_dataset}.{bq_table}`
    WHERE dataset_id IN ({dataset_ids_str})
    """
    query_job_config = bigquery.QueryJobConfig(destination=temp_table_ref)

    query_job = bq_client.query(
        query, job_config=query_job_config, location=bq_location
    )
    query_job.result()  # Wait for the query to finish

    extract_job = bq_client.extract_table(
        temp_table_ref,
        destination_uri,
        location=bq_location,
        job_config=job_config,
    )
    extract_job.result()  # Wait for the extract to finish

    bq_client.delete_table(temp_table_ref)  # Clean up temp table
    return destination_uri


def format_row(row, schema):
    formatted = []
    for i, field in enumerate(schema):
        val = row[i]
        if val is None:
            formatted.append("")
        elif field.mode == "REPEATED":
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


def load_parquet_to_pg(cursor, pg_table, gcs_bucket, gcs_file_path_prefix, schema):
    """Loads parquet shards from GCS into Postgres using COPY."""
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_file_path_prefix))

    for blob in tqdm(blobs, desc=f"Loading shards to {pg_table}", leave=False):
        with blob.open("rb") as f:
            parquet_file = pq.ParquetFile(f)
            for i in range(parquet_file.num_row_groups):
                table = parquet_file.read_row_group(i)
                rows = table.to_pylist()

                tsv_data = io.StringIO()
                for row_dict in rows:
                    row_val = [row_dict.get(field.name) for field in schema]
                    tsv_data.write(format_row(row_val, schema) + "\n")

                tsv_data.seek(0)
                cursor.copy_expert(
                    sql.SQL("COPY {} FROM STDIN WITH (FORMAT TEXT, NULL '')").format(
                        sql.Identifier(pg_table)
                    ),
                    tsv_data,
                )


def update_sync_state(pg_conn, pg_table, dataset_id, timestamp):
    """Updates or inserts the sync state for a specific dataset."""
    with psycopg2.connect(pg_conn) as conn:
        with conn.cursor() as cursor:
            cursor.execute(
                """
                INSERT INTO sync_state (table_name, dataset_id, last_synced_at)
                VALUES (%s, %s, %s)
                ON CONFLICT (table_name, dataset_id) DO UPDATE
                SET last_synced_at = EXCLUDED.last_synced_at;
            """,
                (pg_table, dataset_id, timestamp),
            )


def delete_dataset_from_pg(pg_conn, pg_table, dataset_id):
    """Deletes all rows for a given dataset_id from a Postgres table."""
    with psycopg2.connect(pg_conn) as conn:
        with conn.cursor() as cursor:
            cursor.execute(
                sql.SQL("DELETE FROM {} WHERE dataset_id = %s").format(
                    sql.Identifier(pg_table)
                ),
                (dataset_id,),
            )


def ensure_pg_table_exists(pg_conn, pg_table, bq_schema):
    """Ensures the target table exists in Postgres, creating it if necessary."""
    with psycopg2.connect(pg_conn) as conn:
        with conn.cursor() as cursor:
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


def cleanup_gcs(gcs_bucket, gcs_file_path_prefix):
    """Removes temporary files from GCS."""
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_file_path_prefix))
    for blob in blobs:
        blob.delete()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bq-dataset", required=True)
    parser.add_argument("--bq-table", required=True)
    parser.add_argument("--bq-location", required=True)
    parser.add_argument("--pg-conn", required=True)
    parser.add_argument("--pg-table", required=True)
    parser.add_argument("--gcs-bucket", required=True)
    parser.add_argument(
        "--yes", action="store_true", help="Proceed without confirmation"
    )
    args = parser.parse_args()

    bq_client = bigquery.Client()

    # 1. Fetch current sync states
    logging.info("Fetching current sync states from Postgres...")
    all_states = get_all_sync_states(args.pg_conn)
    current_table_states = all_states.get(args.pg_table, {})

    # 2. Fetch latest timestamps from BQ
    logging.info(f"Fetching latest dataset timestamps from BQ table {args.bq_table}...")
    bq_latest = get_bq_latest_timestamps(
        bq_client, args.bq_dataset, args.bq_table, args.bq_location
    )

    # 3. Determine datasets to process
    to_update = []  # Exists in PG but timestamp in BQ is newer
    to_insert = []  # Does not exist in PG

    for ds_id, bq_ts in bq_latest.items():
        if ds_id in current_table_states:
            pg_ts = current_table_states[ds_id]
            # BQ timestamp is offset-aware usually, but let's be careful.
            # Comparison assumes both are comparable.
            if pg_ts is None or bq_ts > pg_ts.replace(tzinfo=bq_ts.tzinfo):
                to_update.append(ds_id)
        else:
            to_insert.append(ds_id)

    # 4. Summary and confirmation
    total_to_process = len(to_update) + len(to_insert)
    if total_to_process == 0:
        logging.info("Everything is up to date. Nothing to sync.")
        return

    print("\n--- Sync Summary ---")
    print(f"Table: {args.pg_table}")
    print(f"Datasets to update (delete + re-ingest): {len(to_update)}")
    if to_update:
        print(f"  {', '.join(to_update[:10])}{'...' if len(to_update) > 10 else ''}")
    print(f"Datasets to insert (new): {len(to_insert)}")
    if to_insert:
        print(f"  {', '.join(to_insert[:10])}{'...' if len(to_insert) > 10 else ''}")
    print(f"Total datasets to process: {total_to_process}")
    print("--------------------\n")

    if not args.yes:
        confirm = input("Proceed with sync? (y/N): ")
        if confirm.lower() != "y":
            logging.info("Sync cancelled by user.")
            sys.exit(0)

    # 5. Execution
    bq_table_obj = bq_client.get_table(f"{args.bq_dataset}.{args.bq_table}")
    ensure_pg_table_exists(args.pg_conn, args.pg_table, bq_table_obj.schema)

    all_to_process = sorted(to_update + to_insert)
    # We process in chunks to avoid massive GCS exports if many datasets
    chunk_size = 10  # Tuneable

    pbar = tqdm(total=len(all_to_process), desc="Synchronizing datasets")

    for i in range(0, len(all_to_process), chunk_size):
        chunk = all_to_process[i : i + chunk_size]
        gcs_prefix = f"tmp/{args.bq_dataset}_{args.bq_table}_{uuid.uuid4()}"

        try:
            # Delete if update
            for ds_id in chunk:
                if ds_id in to_update:
                    delete_dataset_from_pg(args.pg_conn, args.pg_table, ds_id)

            # Export chunk to GCS
            export_bq_to_gcs(
                bq_client,
                args.bq_dataset,
                args.bq_table,
                args.bq_location,
                args.gcs_bucket,
                gcs_prefix,
                chunk,
            )

            # Load to Postgres
            with psycopg2.connect(args.pg_conn) as conn:
                with conn.cursor() as cursor:
                    load_parquet_to_pg(
                        cursor,
                        args.pg_table,
                        args.gcs_bucket,
                        gcs_prefix,
                        bq_table_obj.schema,
                    )

            # Update sync states
            for ds_id in chunk:
                update_sync_state(args.pg_conn, args.pg_table, ds_id, bq_latest[ds_id])

            pbar.update(len(chunk))

        finally:
            cleanup_gcs(args.gcs_bucket, gcs_prefix)

    pbar.close()
    logging.info("Sync completed successfully.")


if __name__ == "__main__":
    main()
