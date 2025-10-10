import os
import sys
import subprocess
import argparse
from google.cloud import bigquery, storage
from tqdm import tqdm


def get_env_var(name):
    """Gets an environment variable or exits if it's not set."""
    value = os.environ.get(name)
    if not value:
        print(f"Error: Environment variable {name} is not set.", file=sys.stderr)
        sys.exit(1)
    return value


def run_command(command, silent=False):
    """Runs a command and handles errors."""
    if not silent:
        print(f"Running command: {' '.join(command)}")
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        print(f"Error running command: {' '.join(command)}", file=sys.stderr)
        print(process.stderr, file=sys.stderr)
        sys.exit(1)
    if not silent:
        print(process.stdout)
    return process.stdout


def execute_sql_from_gcs(
    project_id, instance_id, database_id, bucket_name, sql_content, sql_file_name
):
    """Uploads a SQL script to GCS and executes it using gcloud sql import sql."""
    print(f"Uploading and executing {sql_file_name}...")
    storage_client = storage.Client(project=project_id)
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(sql_file_name)
    blob.upload_from_string(sql_content)

    sql_uri = f"gs://{bucket_name}/{sql_file_name}"
    run_command(
        [
            "gcloud",
            "sql",
            "import",
            "sql",
            instance_id,
            sql_uri,
            f"--database={database_id}",
            "--quiet",
        ]
    )

    blob.delete()
    print(f"{sql_file_name} executed and cleaned up successfully.")


def main():
    """Main function to migrate data from BigQuery to PostgreSQL."""
    parser = argparse.ArgumentParser(
        description="Migrate data from BigQuery to PostgreSQL."
    )
    parser.add_argument(
        "--limit", type=int, help="Limit the number of CSV files to import."
    )
    args = parser.parse_args()

    print("--- Starting BigQuery to PostgreSQL Migration ---")

    # 1. Read environment variables
    project_id = get_env_var("GCLOUD_PROJECT")
    region = get_env_var("GCLOUD_REGION")
    ps_instance_id = get_env_var("PS_INSTANCE_ID")
    ps_table = get_env_var("PS_TABLE")
    bq_dataset = get_env_var("BQ_DATASET")
    bq_table = get_env_var("BQ_TABLE")
    warehouse_bucket = get_env_var("WAREHOUSE_BUCKET")
    ps_db = os.environ.get("PS_DB", "postgres")

    print("Step 1/5: Environment variables loaded.")

    # 2. BigQuery Export
    print("Step 2/5: Starting BigQuery export to GCS...")
    destination_uri = f"gs://{warehouse_bucket}/bq_export/bq_export_*.csv"
    job_config = bigquery.ExtractJobConfig()
    job_config.print_header = False

    bq_client = bigquery.Client(project=project_id)
    dataset_ref = bigquery.DatasetReference(project_id, bq_dataset)
    table_ref = dataset_ref.table(bq_table)

    extract_job = bq_client.extract_table(
        table_ref,
        destination_uri,
        job_config=job_config,
        location=region,
    )
    extract_job.result()  # Waits for job to complete.
    print("Step 2/5: BigQuery export to GCS completed.")

    # 3. Create Table in Cloud SQL
    print(f"Step 3/5: Creating table '{ps_table}' in Cloud SQL...")
    create_table_sql = f"""DROP TABLE IF EXISTS {ps_table}; CREATE TABLE {ps_table} (
        dataset_id TEXT,
        perturbed_target_symbol TEXT,
        gene TEXT,
        log2foldchange DOUBLE PRECISION,
        padj DOUBLE PRECISION,
        basemean DOUBLE PRECISION,
        max_ingested_at TIMESTAMP WITH TIME ZONE
    );"""
    execute_sql_from_gcs(
        project_id,
        ps_instance_id,
        ps_db,
        warehouse_bucket,
        create_table_sql,
        "create_table.sql",
    )

    # 4. Create Index
    print(f"Step 4/5: Creating index on table '{ps_table}'...")
    create_index_sql = f"CREATE INDEX idx_{ps_table}_multi ON {ps_table} (dataset_id, perturbed_target_symbol, gene);"
    execute_sql_from_gcs(
        project_id,
        ps_instance_id,
        ps_db,
        warehouse_bucket,
        create_index_sql,
        "create_index.sql",
    )

    # 5. Import data from GCS to Cloud SQL and Cleanup
    print("Step 5/5: Importing data from GCS to Cloud SQL and cleaning up...")
    storage_client = storage.Client(project=project_id)
    bucket = storage_client.bucket(warehouse_bucket)
    blobs = list(bucket.list_blobs(prefix="bq_export/"))
    if not blobs:
        raise Exception("BigQuery export created no files.")

    if args.limit:
        print(f"Limiting import to the first {args.limit} files.")
        blobs = blobs[: args.limit]

    print(f"Found {len(blobs)} exported file(s). Importing one by one...")

    for blob in tqdm(blobs, desc="Importing files"):
        file_uri = f"gs://{warehouse_bucket}/{blob.name}"
        run_command(
            [
                "gcloud",
                "sql",
                "import",
                "csv",
                ps_instance_id,
                file_uri,
                f"--database={ps_db}",
                f"--table={ps_table}",
                "--quiet",
            ],
            silent=True,
        )

    print("Data import completed.")

    storage_client = storage.Client(project=project_id)
    bucket = storage_client.bucket(warehouse_bucket)
    # Re-list blobs to be sure we get them all for deletion
    blobs_to_delete = bucket.list_blobs(prefix="bq_export/")
    for blob in blobs_to_delete:
        blob.delete()
    print("GCS export files cleaned up successfully.")

    print("--- Migration script finished. ---")


if __name__ == "__main__":
    main()
