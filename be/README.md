# Pert_Cat_BE

## Local set up

```bash
python3 -m venv fastapi-env
source fastapi-env/bin/activate
pip install fastapi uvicorn elasticsearch[async]
uvicorn main:app --reload
```

## Local Docker build

```bash
docker build -t perturbation-catalogue-be .
docker run -p 8000:8000 perturbation-catalogue-be
```

## Configuring Google Cloud Run deployment

### Set up deployment

1. Go to https://console.cloud.google.com/run.
1. Deploy container → Service → Continuously deploy from a repository (source or function).
1. Set up cloud build.
1. Choose this repository → Next.
1. Build type: Dockerfile; Source location: `/be/Dockerfile` → Save.
1. Service name: `perturbation-catalogue-be`.
1. Choose region.
1. Pick: Allow unauthenticated invokations.
1. Pick: CPU is only allocated during request processing.
1. Container(s), volumes, networking, security → Container(s) → Variables and Secrets → fill in these three secret variables:
   - ES_URL
   - ES_USERNAME
   - ES_PASSWORD
1. Click: Create.

The deployment can then be accessed at the URL shown on the build page.