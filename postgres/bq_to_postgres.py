import argparse
import logging
import uuid
from datetime import datetime

from google.cloud import bigquery, storage
import psycopg2
from psycopg2 import sql
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)


def get_last_synced_at(pg_conn, pg_table):
    """Gets the last synced timestamp from the sync_state table."""
    with psycopg2.connect(pg_conn) as conn:
        with conn.cursor() as cursor:
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS sync_state (
                    table_name TEXT NOT NULL PRIMARY KEY,
                    last_synced_at TIMESTAMP WITHOUT TIME ZONE
                );
            """
            )
            cursor.execute(
                "SELECT last_synced_at FROM sync_state WHERE table_name = %s",
                (pg_table,),
            )
            result = cursor.fetchone()
            return result[0] if result else None


def export_bq_to_gcs(
    bq_client,
    bq_dataset,
    bq_table,
    bq_location,
    gcs_bucket,
    gcs_file_path_prefix,
    last_synced_at,
):
    """Exports data from BigQuery to a GCS bucket."""
    dataset_ref = bq_client.dataset(bq_dataset)
    table_ref = dataset_ref.table(bq_table)
    destination_uri = f"gs://{gcs_bucket}/{gcs_file_path_prefix}-*.csv"

    if last_synced_at:
        # Incremental load: query to a temporary table, then export.
        temp_table_id = f"temp_export_{uuid.uuid4().hex}"
        temp_table_ref = dataset_ref.table(temp_table_id)

        query = f"""
        SELECT *
        FROM `{bq_dataset}.{bq_table}`
        WHERE max_ingested_at > DATETIME('{last_synced_at.isoformat()}')
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
        )
        extract_job.result()  # Wait for the extract to finish

        bq_client.delete_table(temp_table_ref)  # Clean up temp table

    else:
        # Full load: direct export
        extract_job = bq_client.extract_table(
            table_ref,
            destination_uri,
            location=bq_location,
        )
        extract_job.result()

    logging.info(f"Exported data to gs://{gcs_bucket}/{gcs_file_path_prefix}-*.csv")


def load_to_postgres(
    pg_conn,
    pg_table,
    gcs_bucket,
    gcs_file_path_prefix,
    last_synced_at,
    bq_client,
    bq_dataset,
    bq_table_name,
):
    """Loads data from a GCS file into a PostgreSQL table."""
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_file_path_prefix))

    with psycopg2.connect(pg_conn) as conn:
        with conn.cursor() as cursor:
            if last_synced_at:
                for blob in tqdm(blobs, desc="Loading shards"):
                    with blob.open("r") as f:
                        cursor.copy_expert(
                            sql.SQL(
                                "COPY {} FROM STDIN WITH (FORMAT CSV, HEADER)"
                            ).format(sql.Identifier(pg_table)),
                            f,
                        )
            else:
                # Full load
                staging_table = f"staging_{uuid.uuid4().hex}"
                table = bq_client.get_table(f"{bq_dataset}.{bq_table_name}")
                schema = table.schema
                columns = [
                    f"{field.name} {get_pg_type(field.field_type)}" for field in schema
                ]

                cursor.execute(
                    sql.SQL("CREATE TABLE {} ({})").format(
                        sql.Identifier(staging_table),
                        sql.SQL(", ").join(map(sql.SQL, columns)),
                    )
                )

                for blob in tqdm(blobs, desc="Loading shards"):
                    with blob.open("r") as f:
                        cursor.copy_expert(
                            sql.SQL(
                                "COPY {} FROM STDIN WITH (FORMAT CSV, HEADER)"
                            ).format(sql.Identifier(staging_table)),
                            f,
                        )

                cursor.execute(
                    sql.SQL("DROP TABLE IF EXISTS {}").format(sql.Identifier(pg_table))
                )
                cursor.execute(
                    sql.SQL("ALTER TABLE {} RENAME TO {}").format(
                        sql.Identifier(staging_table), sql.Identifier(pg_table)
                    )
                )


def get_pg_type(bq_type):
    """Maps BigQuery types to PostgreSQL types."""
    return {
        "STRING": "TEXT",
        "INTEGER": "BIGINT",
        "FLOAT": "DOUBLE PRECISION",
        "BOOLEAN": "BOOLEAN",
        "TIMESTAMP": "TIMESTAMP WITHOUT TIME ZONE",
        "DATE": "DATE",
    }.get(bq_type, "TEXT")


def update_sync_state(pg_conn, pg_table, bq_client, bq_dataset, bq_table_name):
    """Updates the sync_state table with the latest timestamp."""
    query = f"SELECT MAX(max_ingested_at) FROM `{bq_dataset}.{bq_table_name}`"
    query_job = bq_client.query(query)
    max_ingested_at = list(query_job.result())[0][0]

    if max_ingested_at:
        with psycopg2.connect(pg_conn) as conn:
            with conn.cursor() as cursor:
                cursor.execute(
                    """
                    INSERT INTO sync_state (table_name, last_synced_at)
                    VALUES (%s, %s)
                    ON CONFLICT (table_name) DO UPDATE
                    SET last_synced_at = EXCLUDED.last_synced_at;
                """,
                    (pg_table, max_ingested_at),
                )


def cleanup_gcs(gcs_bucket, gcs_file_path_prefix):
    """Removes the temporary files from GCS."""
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_file_path_prefix))
    for blob in blobs:
        blob.delete()
    logging.info(f"Deleted files with prefix gs://{gcs_bucket}/{gcs_file_path_prefix}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bq-dataset", required=True)
    parser.add_argument("--bq-table", required=True)
    parser.add_argument("--bq-location", required=True)
    parser.add_argument("--pg-conn", required=True)
    parser.add_argument("--pg-table", required=True)
    parser.add_argument("--gcs-bucket", required=True)
    args = parser.parse_args()

    bq_client = bigquery.Client()
    gcs_file_path_prefix = f"tmp/{args.bq_dataset}_{args.bq_table}_{uuid.uuid4()}"

    last_synced_at = get_last_synced_at(args.pg_conn, args.pg_table)

    export_bq_to_gcs(
        bq_client,
        args.bq_dataset,
        args.bq_table,
        args.bq_location,
        args.gcs_bucket,
        gcs_file_path_prefix,
        last_synced_at,
    )
    load_to_postgres(
        args.pg_conn,
        args.pg_table,
        args.gcs_bucket,
        gcs_file_path_prefix,
        last_synced_at,
        bq_client,
        args.bq_dataset,
        args.bq_table,
    )
    update_sync_state(
        args.pg_conn, args.pg_table, bq_client, args.bq_dataset, args.bq_table
    )
    cleanup_gcs(args.gcs_bucket, gcs_file_path_prefix)


if __name__ == "__main__":
    main()
