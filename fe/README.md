# Dash front-end for the Perturbation Catalogue

## Environment variables

For deploying any of the options below, set the following environmental variable:
- `PERTURBATION_CATALOGUE_BE` to point to back-end deployment base URL, starting with https:// and without the trailing slash

## Local deployment
```bash
python3 -m venv env
source env/bin/activate
pip install --quiet -r requirements.txt
python3 app.py
```

## Docker deployment
```bash
docker build -t dash-app .
docker run \
  -p 80:80 \
  -e PERTURBATION_CATALOGUE_BE=${PERTURBATION_CATALOGUE_BE} \
  dash-app
```

## Google Cloud Run deployment

### Continuous deployment

1. Go to https://console.cloud.google.com/run.
1. Deploy container → Service → Continuously deploy from a repository (source or function).
1. Set up cloud build.
1. Choose this repository → Next.
1. Branch: `^main$`; Build type: Dockerfile; Source location: `/fe/Dockerfile` → Save.
1. Service name: `perturbation-catalogue-dash`.
1. Choose region.
1. Pick: Allow unauthenticated invokations.
1. Billing: Request-based.
1. Container(s), volumes, networking, security → Container(s) → Variables and Secrets → fill in environment variables: (see the environment variables section above)
1. Container(s), volumes, networking, security → Container(s) → Container port: 80
1. Click: Create.

The deployment can then be accessed at the URL shown on the build page.

### Set up path trigger

This ensures that the deployment is only updated when something in `fe` is modified.

1. Go to https://console.cloud.google.com/cloud-build/triggers.
1. Edit the perturbation-catalogue-dash trigger.
1. Click on “Show included and ignored files filters”.
1. Set “Included files filters (glob)” to `fe/**`.
1. Click on “Save”.

Then, repeat the steps above with the branch set to `.*` and service name to `perturbation-catalogue-dash-live`. This will create a deployment which will automatically deploy the latest commit, not just from the main branch.
