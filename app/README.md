# Dash front-end for the Perturbation Catalogue

## Local run — direct
```bash
pip install -r requirements.txt && python3 app.py
```

## Local run — Docker
```bash
docker build -t dash-app . && docker run -p 80:80 dash-app
```
