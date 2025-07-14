# Perturbation Catalogue

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

This monorepo contains data processing, back-end, and front-end facilities for the Perturbation Catalogue project. Try it out:

## Production deployment (`main` branch)
* **API** → [https://www.ebi.ac.uk/perturbation-catalogue/mavedb/search](https://www.ebi.ac.uk/perturbation-catalogue/mavedb/search)
* **Website** → [https://www.ebi.ac.uk/perturbation-catalogue](https://www.ebi.ac.uk/perturbation-catalogue)

## Development deployment (`dev` branch)
* **API** → https://perturbation-catalogue-be-dev-959149465821.europe-west2.run.app/mavedb/search
* **Website** → https://perturbation-catalogue-dash-dev-959149465821.europe-west2.run.app


## Live deployment (latest commit)
* **API** → https://perturbation-catalogue-be-live-959149465821.europe-west2.run.app/mavedb/search
* **Website** → https://perturbation-catalogue-dash-live-959149465821.europe-west2.run.app

## Set up
To automatically enforce the code style before committing, run:

```python
pip install pre-commit black
pre-commit install
```

You only need to run these commands once.
