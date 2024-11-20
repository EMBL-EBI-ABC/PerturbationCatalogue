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
docker build -t fastapi-hello-world .
docker run -p 8000:8000 fastapi-hello-world
```

## Configuring Google Cloud Run deployment

### Set up deployment

1. Go to: https://console.cloud.google.com/run
1. Deploy container
1. Service
1. Continuously deploy from a repository (source or function)
1. Set up cloud build
1. (If needed) Authenticate to GitHub
1. Choose this repository
1. I understand that GitHub content...
1. Next
1. Build type: Dockerfile
1. Save
1. Choose region
1. Allow unauthenticated invokations
1. CPU is only allocated during request processing
1. Create

The deployment can then be accessed at the URL shown on the build page.

### Add or modify secret variables

1. Go to: https://console.cloud.google.com/run
1. Open the back-end deployment.
1. Edit & Deploy New Revision.
1. Variables and Secrets
1. Add variable(s)
1. Deploy

Secret variables required for the service to work:

- ES_URL
- ES_USERNAME
- ES_PASSWORD
