# Migrating the SQL instance from dev to prod

## 1. Create a manual backup in dev
```bash
dev_secrets
gcloud sql backups create \
  --instance=${PG_INSTANCE_ID} \
  --project=${GCLOUD_PROJECT} \
  --description="cross-project migration backup"
export BACKUP_ID=$(gcloud sql backups list \
  --instance="${PG_INSTANCE_ID}" \
  --project="${GCLOUD_PROJECT}" \
  --sort-by="~windowStartTime" \
  --limit=1 \
  --format="value(id)")
export DATABASE_VERSION=$(gcloud sql instances describe ${PG_INSTANCE_ID} --project=${GCLOUD_PROJECT} --format="value(databaseVersion)")
export BACKUP_INSTANCE=${PG_INSTANCE_ID}
export BACKUP_PROJECT=${GCLOUD_PROJECT}
```

## 2. Create a new clean instance in prod
```bash
prod_secrets
export NEW_INSTANCE_ID=perturbation-catalogue-data-$(date +%Y-%m-%d)
gcloud sql instances create ${NEW_INSTANCE_ID} \
  --project=${GCLOUD_PROJECT} \
  --database-version=${DATABASE_VERSION} \
  --region=${GCLOUD_REGION} \
  --tier=db-custom-4-16384 \
  --edition=ENTERPRISE \
  --availability-type=ZONAL \
  --storage-type=SSD \
  --storage-size=500GB \
  --storage-auto-increase \
  --retained-transaction-log-days=7 \
  --enable-point-in-time-recovery \
  --maintenance-window-day=SUN \
  --maintenance-window-hour=0 \
  --maintenance-release-channel=production \
  --database-flags=temp_file_limit=104857600,maintenance_work_mem=4194304,max_parallel_maintenance_workers=4 \
  --network=projects/${GCLOUD_PROJECT}/global/networks/default \
  --deletion-protection
gcloud sql users set-password postgres \
  --instance=${NEW_INSTANCE_ID} \
  --project=${GCLOUD_PROJECT} \
  --password="${PG_PASSWORD}"
```

## 3. Restore the backup to the new instance
```bash
prod_secrets
gcloud sql backups restore ${BACKUP_ID} \
  --backup-instance=${BACKUP_INSTANCE} \
  --backup-project=${BACKUP_PROJECT} \
  --restore-instance=${NEW_INSTANCE_ID} \
  --project=${GCLOUD_PROJECT}
# Might need to run a manual wait command suggested by the command above after that.
```
