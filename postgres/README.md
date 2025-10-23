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

## 1. Create a Google Cloud VM
```bash
gcloud compute instances create bq-to-pg-projector \
    --project=${GCLOUD_PROJECT} \
    --zone=${GCLOUD_ZONE} \
    --machine-type=e2-medium \
    --network=default \
    --scopes=https://www.googleapis.com/auth/cloud-platform
```

## 2. Copy files and SSH into the VM
```bash
gcloud compute scp --project=${GCLOUD_PROJECT} --zone=${GCLOUD_ZONE} postgres/requirements.txt postgres/bq_to_postgres.py   bq-to-pg-projector:~
gcloud compute ssh bq-to-pg-projector --project=${GCLOUD_PROJECT} --zone=${GCLOUD_ZONE}
```

## 3. Install dependencies
```bash
sudo apt install -y python3-pip python3-venv
python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
```

## 4. Run the script
```bash
python3 bq_to_postgres.py \
    --bq-dataset <your-bq-dataset> \
    --bq-table <your-bq-table> \
    --bq-location <your-bq-location> \
    --pg-conn "postgresql://<user>:<password>@<private-ip>:<port>/<dbname>" \
    --pg-table <your-pg-table> \
    --gcs-bucket <your-gcs-bucket>
```
