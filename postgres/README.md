# PostgreSQL database

Create the database using the link: https://console.cloud.google.com/sql/instances/create;engine=PostgreSQL;template=POSTGRES_ENTERPRISE_PLUS_DATA_CACHE_ENABLED_DEV_TEMPLATE

Change the settings:
* Cloud SQL edition: Enterprise
* Edition preset: Development
* Instance ID: `$PS_INSTANCE_ID`
* Password: `$PS_PASSWORD`
* Zonal availability: Single zone
* Customise your instance
  + Machine configuration: Dedicated core, 4 vCPU, 16 GB
  + Data protection: disable Automated daily backups and Enable point-in-time recovery

## Data Migration from BigQuery

These instructions explain how to migrate the data from the BigQuery table specified in the environment variables to the PostgreSQL instance. This will use the default `postgres` database.

1.  **Connect to the PostgreSQL instance:**

    ```bash
    gcloud sql connect $PS_INSTANCE_ID --user=postgres
    ```

2.  **Create the table to hold the BigQuery data.** Make sure the `PS_TABLE` environment variable is set.

    ```sql
    CREATE TABLE "$PS_TABLE" (
        dataset_id TEXT,
        perturbed_target_symbol TEXT,
        gene TEXT,
        log2foldchange DOUBLE PRECISION,
        padj DOUBLE PRECISION,
        basemean DOUBLE PRECISION,
        max_ingested_at TIMESTAMP WITH TIME ZONE
    );
    ```

3.  **Export the BigQuery table to a CSV file in a Google Cloud Storage bucket.**

    First, if you don't have a GCS bucket, create one. Replace `your-bucket-name` with a unique bucket name.

    ```bash
    gsutil mb gs://your-bucket-name
    ```

    Now, export the data. This command should be run in a separate terminal.

    ```bash
    bq --location=$GCLOUD_REGION extract \
      --destination_format CSV \
      --field_delimiter ',' \
      --print_header=true \
      "$BQ_DATASET.$BQ_TABLE" \
      "gs://your-bucket-name/bq_export.csv"
    ```

4.  **Download the CSV file from GCS:**

    ```bash
    gsutil cp "gs://your-bucket-name/bq_export.csv" .
    ```

5.  **Import the data into the PostgreSQL table.**

    Run this command in the `psql` shell you connected to in step 1.

    ```sql
    \copy "$PS_TABLE" FROM 'bq_export.csv' WITH (FORMAT CSV, HEADER true);
    ```

6.  **Create indexes on the table for faster queries:**

    ```sql
    CREATE INDEX idx_ps_table_multi ON "$PS_TABLE" (dataset_id, perturbed_target_symbol, gene);
    ```
