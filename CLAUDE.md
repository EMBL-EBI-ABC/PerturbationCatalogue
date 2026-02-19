# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Perturbation Catalogue is a data platform for searching and analyzing genetic perturbation experiments (Perturb-seq, CRISPR screens, MAVE). It's a monorepo with three main components: a FastAPI backend (`be/`), a Dash frontend (`fe/`), and a data warehouse layer (`dwh/`).

- **Live site:** https://www.ebi.ac.uk/perturbation-catalogue
- **API docs:** https://www.ebi.ac.uk/perturbation-catalogue/v1/docs

## Architecture

```
Frontend (Dash, port 8050)  →  Backend (FastAPI, port 8000)  →  Elasticsearch (search/facets)
                                                              →  PostgreSQL (analytics/details)
                                                              →  Google Gemini (AI chat)

Data pipeline: BigQuery (raw) → dbt transforms → Elasticsearch + PostgreSQL
```

**Backend (`be/`):** Async FastAPI app. `main.py` handles app setup, config, and Elasticsearch queries. `data_query.py` contains all search/query logic with modality-specific routers. `ai_chat.py` implements Google Gemini chat with function calling. Connection pools (asyncpg, AsyncElasticsearch) are initialized in the lifespan handler and stored in `db_pools` dict.

**Frontend (`fe/`):** Dash multi-page app. `app.py` is the entry point. Pages are in `pages/` and auto-register via `dash.register_page()`. `utils.py` has API call helpers and caching. Styling uses dash-bootstrap-components + custom CSS in `assets/`.

**Data warehouse (`dwh/`):** dbt models in `bq_dbt/` transform BigQuery data. `postgres/bq_to_postgres.py` syncs to PostgreSQL. `bq_to_es_projector/` indexes into Elasticsearch.

## Commands

### Setup
```bash
# Pre-commit hooks (one time)
pip install pre-commit black
pre-commit install

# Backend
cd be && python3 -m venv fastapi-env && source fastapi-env/bin/activate && pip install -r requirements.txt

# Frontend
cd fe && python3 -m venv env && source env/bin/activate && pip install -r requirements.txt
```

### Running locally
```bash
# Backend (requires PG_HOST, PG_PORT, PG_USER, PG_PASSWORD, PG_DB, ES_URL, ES_USERNAME, ES_PASSWORD env vars)
cd be && uvicorn main:app --reload

# Frontend (requires PERTURBATION_CATALOGUE_BE env var pointing to backend URL)
cd fe && python3 app.py
```

### Code formatting
```bash
black be/ fe/              # Manual formatting
pre-commit run --all-files # Run all pre-commit hooks
```

### dbt (data warehouse)
```bash
cd dwh/bq_dbt
dbt run                           # Run all models
dbt run --select +dataset_summary # Run with dependencies
dbt test                          # Data quality tests
```

## LLM Agent Rules (from prompts/README.md)

- Do NOT run `git`, `ruff`, `pre-commit`, or try to tidy up the repository.
- If tasked with writing code, simply do that and nothing else.

## Key Patterns

- **Modality mapping:** API URL modalities are lowercase (`perturb-seq`, `crispr-screen`, `mave`) but Elasticsearch values are mixed case (`Perturb-seq`, `CRISPR screen`, `MAVE`). See `be-specification.md`.
- **Async everywhere:** Backend uses `async/await` for all DB and network calls. Connection pools are shared via `db_pools`.
- **Elasticsearch indexes:** `dataset-summary` (search), `target-summary` (targets), `landing-page-summary` (homepage stats).
- **Code style:** Black formatter enforced via pre-commit. No other linters.
- **Git branching:** `dev` is the default PR target branch. Feature branches use `name/feature-name` convention.

## Specifications

Detailed specs for backend and frontend live in `prompts/`:
- `prompts/be-specification.md` - Backend API implementation details
- `prompts/fe-specification.md` - Frontend implementation details
- `prompts/data-model.md` - Data model specification
- `prompts/postgres-summary.md` - PostgreSQL schema overview
