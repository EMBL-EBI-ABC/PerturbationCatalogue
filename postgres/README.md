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

This section explains how to migrate data from a BigQuery table to the PostgreSQL instance using an automated script.

### Prerequisites

1.  **Set Environment Variables:**

    Export the following environment variables in your shell:

    ```bash
    export GCP_PROJECT="your-gcp-project-id"
    export GCLOUD_REGION="your-gcp-region" # e.g., us-central1
    export PS_INSTANCE_NAME="your-gcp-project-id:your-gcp-region:your-ps-instance-id"
    export PS_PASSWORD="your-postgres-password"
    export PS_TABLE="your-postgres-table-name"
    export BQ_DATASET="your-bigquery-dataset"
    export BQ_TABLE="your-bigquery-table"
    export WAREHOUSE_BUCKET="your-gcs-bucket-for-staging"
    ```

2.  **Install Dependencies:**

    The script requires several Python libraries. Install them using pip:

    ```bash
    pip install -r requirements.txt
    ```

### Running the Migration

Execute the Python script to start the migration:

```bash
python3 perturb_seq.py
```

The script will perform the following steps:
1.  Export the specified BigQuery table to a CSV file in your GCS bucket.
2.  Connect to your Cloud SQL for PostgreSQL instance.
3.  Create a new table (or drop and recreate if it exists).
4.  Stream the data from the GCS CSV file directly into the PostgreSQL table without saving it locally.
5.  Create an index on the new table to improve query performance.
6.  Delete the temporary CSV file from your GCS bucket.

The script will display its progress as it executes each step.

