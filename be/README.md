# Perturbation Catalogue back-end

## Environment variables

For deploying any of the options below, set the following environmental variables:
- Elastic
  - ES_URL
  - ES_USERNAME
  - ES_PASSWORD
- Postgres
  - PS_HOST
  - PS_PORT
  - PS_USER
  - PS_PASSWORD
  - PS_DB

## Local deployment

```bash
python3 -m venv fastapi-env
source fastapi-env/bin/activate
pip install -r requirements.txt
uvicorn main:app --reload
```

## Docker deployment

```bash
docker build -t perturbation-catalogue-be .
docker run \
  -p 8000:8080 \
  -e ES_URL=${ES_URL} \
  -e ES_USERNAME=${ES_USERNAME} \
  -e ES_PASSWORD=${ES_PASSWORD} \
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
