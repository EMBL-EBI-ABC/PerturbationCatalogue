# Perturbation Catalogue back-end

## Environment variables
Before running either of the deployment options below, run `pc_secrets dev`.
`ES_INDEX_SET` is optional and selects a suffixed Elasticsearch index set; it
defaults to the standard indexes when unset.
`RELEASE_BUCKET` is the manually published GCS bucket containing release
artifacts. The runtime service account must be allowed to sign URLs and read
objects in this bucket.

For the default Cloud Run identity, grant the following permissions. The
runtime service account needs to read the bucket and sign itself; the Cloud
Run service agent also needs to sign the runtime service account because Cloud
Run delegates the signing request through it:

```bash
PROJECT_NUMBER=$(gcloud projects describe "$GCLOUD_PROJECT" --format='value(projectNumber)')
RUNTIME_SA="${PROJECT_NUMBER}-compute@developer.gserviceaccount.com"
CLOUD_RUN_AGENT="service-${PROJECT_NUMBER}@serverless-robot-prod.iam.gserviceaccount.com"

gcloud storage buckets add-iam-policy-binding "gs://$RELEASE_BUCKET" \
  --member="serviceAccount:$RUNTIME_SA" \
  --role=roles/storage.objectViewer
gcloud iam service-accounts add-iam-policy-binding "$RUNTIME_SA" \
  --member="serviceAccount:$RUNTIME_SA" \
  --role=roles/iam.serviceAccountTokenCreator
gcloud iam service-accounts add-iam-policy-binding "$RUNTIME_SA" \
  --member="serviceAccount:$CLOUD_RUN_AGENT" \
  --role=roles/iam.serviceAccountTokenCreator
```

For local ADC, grant the local user `roles/storage.objectViewer` on the
release bucket and `roles/iam.serviceAccountTokenCreator` on the runtime
service account.

Unfiltered dataset downloads redirect to seven-day V4 signed URLs in that
bucket; filtered CSV downloads continue to run through the API.

If you are running locally and as such connecting to Postges externally, allow connections from your IP:
* https://console.cloud.google.com/sql/instances
* Go to instance
* Connections
* Networking
* Authorised networks → Use my IP → Save

## Local deployment

```bash
python3 -m venv fastapi-env
source fastapi-env/bin/activate
pip install -r requirements.txt
uvicorn main:app --reload
```

## Local tests

See [tests/README.md](tests/README.md). Run with `pc_secrets -v dev && be/fastapi-env/bin/pytest -q be/tests`.

## Docker deployment

```bash
docker build -t perturbation-catalogue-be .
docker run \
  -p 8000:8080 \
  -e ES_URL=${ES_URL} \
  -e ES_USERNAME=${ES_USERNAME} \
  -e ES_PASSWORD=${ES_PASSWORD} \
  -e ES_INDEX_SET=${ES_INDEX_SET:-} \
  -e PS_HOST=${PS_HOST} \
  -e PS_PORT=${PS_PORT} \
  -e PS_USER=${PS_USER} \
  -e PS_PASSWORD=${PS_PASSWORD} \
  -e PS_DB=${PS_DB} \
  perturbation-catalogue-be
```

## Google Cloud Run deployment

1. Go to https://console.cloud.google.com/run.
1. Deploy container → Continuously deploy from a repository (source or function).
1. Set up cloud build.
1. Choose this repository → Next.
1. Branch: `^main$`; Build type: Dockerfile; Source location: `/be/Dockerfile` → Save.
1. Service name: `perturbation-catalogue-be`.
1. Choose region.
1. Pick: Allow unauthenticated invokations.
1. Billing: Request-based.
1. Container(s), volumes, networking, security → Containers → Variables & Secrets → fill in environment variables: (see the environment variables section above)
1. Networking → Enable checkbox for "Connect to VPC for outbound traffic". Leave default settings, namely; Send traffic directly to a VPC; default network and subnet.
1. Click: Create.

The deployment can then be accessed at the URL shown on the build page.

Set up path trigger:

1. Go to https://console.cloud.google.com/cloud-build/triggers.
1. Edit the perturbation-catalogue-be trigger.
1. Click on “Show included and ignored files filters”.
1. Set “Included files filters (glob)” to `be/**`.
1. Click on “Save”.

Then, repeat the steps above with the following changes:
* Branch: `^dev$`
* Service name: `perturbation-catalogue-be-dev`

Then, repeat the steps above with the following changes:
* Branch: `^main$|^dev$`
* Reverse regex: checked
* Service name: `perturbation-catalogue-be-live`

This will create a deployment which will automatically deploy the latest commit.
