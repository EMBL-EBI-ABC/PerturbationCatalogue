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
  + Data protection: disable Automated daily backups and Enable point-in-time recovery

# Migrate data from BigQuery

Set environment variables (from project wide secrets):
* GCLOUD_PROJECT
* GCLOUD_REGION
* PS_INSTANCE_ID
* PS_PASSWORD
* PS_TABLE
* BQ_DATASET
* BQ_TABLE
* WAREHOUSE_BUCKET

Install dependencies:
```bash
pip install -r requirements.txt
```

Migrate the data:

```bash
python3 perturb_seq.py
```
