# Dash front-end for the Perturbation Catalogue

## Local run — direct
```bash
pip install -r requirements.txt && python3 app.py
```

## Local run — Docker
```bash
docker build -t dash-app . && docker run -p 80:80 dash-app
```

## Google Cloud Run deployment

### Continuous deployment

1. Go to https://console.cloud.google.com/run.
1. Deploy container → Service → Continuously deploy from a repository (source or function).
1. Set up cloud build.
1. Choose this repository → Next.
1. Build type: Dockerfile; Source location: `/fe/Dockerfile` → Save.
1. Service name: `perturbation-catalogue-dash`.
1. Choose region.
1. Pick: Allow unauthenticated invokations.
1. Billing: Request-based.
1. Container(s), volumes, networking, security → Container(s) → Container port: 80
1. Click: Create.

The deployment can then be accessed at the URL shown on the build page.

### Set up path trigger

This ensures that the deployment is only updated when something in `app` is modified.

1. Go to https://console.cloud.google.com/cloud-build/triggers.
1. Edit the perturbation-catalogue-app trigger.
1. Click on “Show included and ignored files filters”.
1. Set “Included files filters (glob)” to `app/**`.
1. Clik on “Save”.