import argparse
import logging
import uuid
import io
import json

from google.cloud import bigquery, storage
import psycopg2
from psycopg2 import sql
from tqdm import tqdm
import pyarrow.parquet as pq

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
    """Exports data from BigQuery to a GCS bucket in Parquet format."""
    dataset_ref = bq_client.dataset(bq_dataset)
    table_ref = dataset_ref.table(bq_table)
    destination_uri = f"gs://{gcs_bucket}/{gcs_file_path_prefix}-*.parquet"

    job_config = bigquery.ExtractJobConfig(destination_format="PARQUET")

    if last_synced_at:
        # Incremental load: query to a temporary table, then export.
        temp_table_id = f"temp_export_{uuid.uuid4().hex}"
        temp_table_ref = dataset_ref.table(temp_table_id)

        query = f"""
        SELECT *
        FROM `{bq_dataset}.{bq_table}`
        WHERE max_ingested_at > TIMESTAMP('{last_synced_at.isoformat()}')
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

    else:
        # Full load: direct export
        extract_job = bq_client.extract_table(
            table_ref,
            destination_uri,
            location=bq_location,
            job_config=job_config,
        )
        extract_job.result()

    logging.info(f"Exported data to gs://{gcs_bucket}/{gcs_file_path_prefix}-*.parquet")


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
    """Loads data from a GCS Parquet file into a PostgreSQL table."""
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(gcs_bucket)
    blobs = list(bucket.list_blobs(prefix=gcs_file_path_prefix))

    table = bq_client.get_table(f"{bq_dataset}.{bq_table_name}")
    schema = table.schema

    def format_row(row, schema):
        formatted = []
        for i, field in enumerate(schema):
            val = row[i]
            if val is None:
                formatted.append("")
            elif field.mode == "REPEATED":
                # Convert list to Postgres array literal: {"val1", "val2"}
                # We use json.dumps to handle basic quoting for strings.
                inner = [str(x) if not isinstance(x, str) else x for x in val]
                # Simple quoting for postgres array format
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

    def stream_parquet_to_pg(cursor, target_table, blobs, schema):
        for blob in tqdm(blobs, desc=f"Loading shards to {target_table}"):
            with blob.open("rb") as f:
                parquet_file = pq.ParquetFile(f)
                for i in range(parquet_file.num_row_groups):
                    table = parquet_file.read_row_group(i)
                    rows = table.to_pylist()

                    # Convert rows to TSV format for COPY
                    tsv_data = io.StringIO()
                    for row_dict in rows:
                        # Convert dict to ordered list based on schema
                        row_val = [row_dict.get(field.name) for field in schema]
                        tsv_data.write(format_row(row_val, schema) + "\n")

                    tsv_data.seek(0)
                    cursor.copy_expert(
                        sql.SQL(
                            "COPY {} FROM STDIN WITH (FORMAT TEXT, NULL '')"
                        ).format(sql.Identifier(target_table)),
                        tsv_data,
                    )

    with psycopg2.connect(pg_conn) as conn:
        with conn.cursor() as cursor:
            if last_synced_at:
                stream_parquet_to_pg(cursor, pg_table, blobs, schema)
            else:
                # Full load
                staging_table = f"staging_{uuid.uuid4().hex}"
                columns = [f"{field.name} {get_pg_type(field)}" for field in schema]

                cursor.execute(
                    sql.SQL("CREATE TABLE {} ({})").format(
                        sql.Identifier(staging_table),
                        sql.SQL(", ").join(map(sql.SQL, columns)),
                    )
                )

                stream_parquet_to_pg(cursor, staging_table, blobs, schema)

                cursor.execute(
                    sql.SQL("DROP TABLE IF EXISTS {}").format(sql.Identifier(pg_table))
                )
                cursor.execute(
                    sql.SQL("ALTER TABLE {} RENAME TO {}").format(
                        sql.Identifier(staging_table), sql.Identifier(pg_table)
                    )
                )


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


def update_sync_state(
    pg_conn, pg_table, bq_client, bq_dataset, bq_table_name, bq_location
):
    """Updates the sync_state table with the latest timestamp."""
    query = f"SELECT MAX(max_ingested_at) FROM `{bq_dataset}.{bq_table_name}`"
    query_job = bq_client.query(query, location=bq_location)
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
        args.pg_conn,
        args.pg_table,
        bq_client,
        args.bq_dataset,
        args.bq_table,
        args.bq_location,
    )
    cleanup_gcs(args.gcs_bucket, gcs_file_path_prefix)


if __name__ == "__main__":
    main()
