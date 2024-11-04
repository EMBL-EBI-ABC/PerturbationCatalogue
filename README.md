# Pert_Cat_FE

## Local build
```bash
ng serve
```

## BigQuery API configuration
Generate API key:
* Go to https://console.cloud.google.com/apis/credentials
* Create Credentials → API key
* Restrict key to "BigQuery API"

Publish dataset:
* Go to dataset in BigQuery
* Share dataset
* Add allUsers with the BigQuery Data Viewer permission
