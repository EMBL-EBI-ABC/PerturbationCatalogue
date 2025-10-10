# Create PostgreSQL instance

Create the database using the link: https://console.cloud.google.com/sql/instances/create;engine=PostgreSQL;template=POSTGRES_ENTERPRISE_PLUS_DATA_CACHE_ENABLED_DEV_TEMPLATE

Change the settings:
* Cloud SQL edition: Enterprise
* Edition preset: Development
* Instance ID: `$PS_INSTANCE_ID`
* Password: `$PS_PASSWORD`
* Zonal availability: Single zone
* Customise your instance
  + Machine configuration: Dedicated core, 4 vCPU, 16 GB
  + Data protection: disable "Automated daily backups" and "Enable point-in-time recovery"

# Migrate data from BigQuery

This section explains how to migrate data from a BigQuery table to the PostgreSQL instance using a Python script that orchestrates cloud-native services.

### Prerequisites

1.  **Google Cloud SDK (`gcloud`):**
    Ensure you have the Google Cloud SDK installed and authenticated:
    ```bash
    gcloud auth application-default login
    ```

2.  **Cloud SQL Permissions:**
    Your Cloud SQL instance needs permission to read from your GCS bucket. Run the following command once to grant this permission. This command finds your instance's service account and gives it the necessary role on your warehouse bucket.

    ```bash
    SERVICE_ACCOUNT_EMAIL=$(gcloud sql instances describe $PS_INSTANCE_ID --format='value(serviceAccountEmailAddress)')
    gsutil iam ch serviceAccount:$SERVICE_ACCOUNT_EMAIL:objectAdmin gs://$WAREHOUSE_BUCKET
    ```

3.  **Set Environment Variables:**
    Export the following environment variables:
    ```bash
    export GCLOUD_PROJECT="your-gcp-project-id"
    export GCLOUD_REGION="your-gcp-region"
    export PS_INSTANCE_ID="your-ps-instance-id"
    export PS_PASSWORD="your-postgres-password"
    export PS_TABLE="your-postgres-table-name"
    export BQ_DATASET="your-bigquery-dataset"
    export BQ_TABLE="your-bigquery-table"
    export WAREHOUSE_BUCKET="your-gcs-bucket-for-staging"
    # Optional:
    # export PS_USER="your-postgres-user" # defaults to 'postgres'
    # export PS_DB="your-postgres-db"   # defaults to 'postgres'
    ```

4.  **Install Dependencies:**
    Install the required Python libraries:
    ```bash
    pip install -r requirements.txt
    ```

### Running the Migration

Execute the Python script to start the migration:

```bash
python3 perturb_seq.py
```

To test with the first N chunks only:

```bash
python3 perturb_seq.py --limit 10
```
