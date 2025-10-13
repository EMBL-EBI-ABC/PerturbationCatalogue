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

## 1. Authenticate Google Cloud SDK
```bash
gcloud auth application-default login
```

## 2. Set up Cloud SQL permissions
This is required so that the Cloud SQL instance can read from the bucket. These commands need to be run once.
```bash
SERVICE_ACCOUNT_EMAIL=$(gcloud sql instances describe $PS_INSTANCE_ID --format='value(serviceAccountEmailAddress)')
gsutil iam ch serviceAccount:$SERVICE_ACCOUNT_EMAIL:objectAdmin gs://$WAREHOUSE_BUCKET
```

## 3. Export environment variables
These should come from project wide secrets:
```bash
export GCLOUD_PROJECT="your-gcp-project-id"
export GCLOUD_REGION="your-gcp-region"
export PS_INSTANCE_ID="your-ps-instance-id"
export PS_PASSWORD="your-postgres-password"
export PS_TABLE="your-postgres-table-name"
export BQ_DATASET="your-bigquery-dataset"
export BQ_TABLE="your-bigquery-table"
export WAREHOUSE_BUCKET="your-gcs-bucket-for-staging"
```

## 4. Install dependencies
```bash
pip install -r requirements.txt
```

## 5. Run the migration
```bash
python3 perturb_seq.py
```

To test with the first N chunks only:

```bash
python3 perturb_seq.py --limit 10
```
