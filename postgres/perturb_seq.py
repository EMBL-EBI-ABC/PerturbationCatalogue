import os
import sys
from google.cloud import bigquery, storage
from google.cloud.sql.connector import Connector
import psycopg2


def get_env_var(name):
    """Gets an environment variable or exits if it's not set."""
    value = os.environ.get(name)
    if not value:
        print(f"Error: Environment variable {name} is not set.", file=sys.stderr)
        sys.exit(1)
    return value


def main():
    """Main function to migrate data from BigQuery to PostgreSQL."""
    print("--- Starting BigQuery to PostgreSQL Migration ---")

    # 1. Read environment variables
    try:
        project_id = get_env_var("GCLOUD_PROJECT")
        region = get_env_var("GCLOUD_REGION")
        ps_instance_name = get_env_var(
            "PS_INSTANCE_ID"
        )  # e.g. my-project:us-central1:my-instance
        ps_password = get_env_var("PS_PASSWORD")
        ps_table = get_env_var("PS_TABLE")
        bq_dataset = get_env_var("BQ_DATASET")
        bq_table = get_env_var("BQ_TABLE")
        warehouse_bucket = get_env_var("WAREHOUSE_BUCKET")
        ps_user = os.environ.get("PS_USER", "postgres")
        ps_db = os.environ.get("PS_DB", "postgres")
    except SystemExit as e:
        return  # Error message is already printed

    print("Step 1/7: Environment variables loaded.")

    # 2. BigQuery Export
    print("Step 2/7: Starting BigQuery export to GCS...")
    try:
        bq_client = bigquery.Client(project=project_id)
        destination_uri = f"gs://{warehouse_bucket}/bq_export.csv"
        dataset_ref = bigquery.DatasetReference(project_id, bq_dataset)
        table_ref = dataset_ref.table(bq_table)

        extract_job = bq_client.extract_table(
            table_ref,
            destination_uri,
            location=region,
        )
        extract_job.result()  # Waits for job to complete.
        print("Step 2/7: BigQuery export to GCS completed.")
    except Exception as e:
        print(f"Error during BigQuery export: {e}", file=sys.stderr)
        return

    # 3. PostgreSQL Connection
    print("Step 3/7: Connecting to PostgreSQL instance...")
    try:
        connector = Connector()
        conn = connector.connect(
            ps_instance_name,
            "psycopg2",
            user=ps_user,
            password=ps_password,
            db=ps_db,
        )
        cursor = conn.cursor()
        print("Step 3/7: PostgreSQL connection successful.")
    except Exception as e:
        print(f"Error connecting to PostgreSQL: {e}", file=sys.stderr)
        return

    # 4. PostgreSQL Table Creation
    print(f"Step 4/7: Creating table '{ps_table}' if it doesn't exist...")
    try:
        cursor.execute(f"DROP TABLE IF EXISTS {ps_table};")
        create_table_sql = f"""
        CREATE TABLE {ps_table} (
            dataset_id TEXT,
            perturbed_target_symbol TEXT,
            gene TEXT,
            log2foldchange DOUBLE PRECISION,
            padj DOUBLE PRECISION,
            basemean DOUBLE PRECISION,
            max_ingested_at TIMESTAMP WITH TIME ZONE
        );
        """
        cursor.execute(create_table_sql)
        conn.commit()
        print(f"Step 4/7: Table '{ps_table}' created successfully.")
    except Exception as e:
        print(f"Error creating PostgreSQL table: {e}", file=sys.stderr)
        conn.rollback()
        return

    # 5. Data Streaming from GCS to PostgreSQL
    print("Step 5/7: Starting data streaming from GCS to PostgreSQL...")
    try:
        storage_client = storage.Client(project=project_id)
        bucket = storage_client.bucket(warehouse_bucket)
        blob = bucket.blob("bq_export.csv")

        with blob.open("r") as gcs_file:
            next(gcs_file)  # Skip header row
            cursor.copy_from(gcs_file, ps_table, sep=",")

        conn.commit()
        print("Step 5/7: Data streaming completed successfully.")
    except Exception as e:
        print(f"Error during data streaming: {e}", file=sys.stderr)
        conn.rollback()
        return

    # 6. PostgreSQL Index Creation
    print(f"Step 6/7: Creating index on table '{ps_table}'...")
    try:
        index_sql = f"CREATE INDEX idx_{ps_table}_multi ON {ps_table} (dataset_id, perturbed_target_symbol, gene);"
        cursor.execute(index_sql)
        conn.commit()
        print("Step 6/7: Index created successfully.")
    except Exception as e:
        print(f"Error creating index: {e}", file=sys.stderr)
        conn.rollback()
        return

    # 7. Cleanup GCS
    print("Step 7/7: Deleting temporary CSV file from GCS...")
    try:
        blob.delete()
        print("Step 7/7: GCS cleanup successful.")
    except Exception as e:
        print(f"Warning: Failed to delete GCS file: {e}", file=sys.stderr)

    finally:
        if "conn" in locals() and conn:
            conn.close()
        if "connector" in locals() and connector:
            connector.close()
        print("--- Migration script finished. ---")


if __name__ == "__main__":
    main()
