# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Perturbation Catalogue is a data platform for searching and analyzing genetic perturbation experiments (Perturb-seq, CRISPR screens, MAVE). It's a monorepo with three main components: a FastAPI backend (`be/`), a Dash frontend (`fe/`), and a data warehouse layer (`dwh/`).

- **Live site:** https://www.ebi.ac.uk/perturbation-catalogue
- **API docs:** https://www.ebi.ac.uk/perturbation-catalogue/v1/docs

## Architecture

```
Frontend (Dash, port 8050 local / 80 Docker)
  → Backend (FastAPI, port 8000 local / 8080 Docker)
      → Elasticsearch (search/facets)
      → PostgreSQL (analytics/details)
      → Google Gemini (AI chat via SSE)

Data pipeline: BigQuery (raw) → dbt transforms → Elasticsearch + PostgreSQL
Deployment: Google Cloud Run (main → prod, dev → dev)
```

**Backend (`be/`):** Async FastAPI app (Python 3.12). `main.py` handles app setup, Elasticsearch query logic (`build_elasticsearch_query`, `perform_search`, `parse_filters_from_params`), and the lifespan handler that initializes connection pools. `data_query.py` contains modality-specific routers and all PostgreSQL query logic. `ai_chat.py` implements Gemini chat with SSE streaming and optional Open Targets MCP integration. `models.py` has Pydantic models. Connection pools (asyncpg, AsyncElasticsearch) are initialized in the lifespan handler and shared via the `db_pools` dict (defined in `data_query.py`, populated in `main.py`, also used by `ai_chat.py`).

**Frontend (`fe/`):** Dash multi-page app mounted at `/perturbation-catalogue/`. `app.py` is the entry point (Gunicorn in production via `server = app.server`). Pages in `pages/` auto-register via `dash.register_page()`. `utils.py` has synchronous `requests`-based API helpers and in-memory caching. `assets/chat.js` handles the entire AI chat UI client-side (SSE consumption, Plotly/PDBe-Molstar/Cytoscape rendering). Styling uses dash-bootstrap-components + `assets/custom.css`.

**Data warehouse (`dwh/`):** dbt models in `bq_dbt/` transform BigQuery data. `postgres/bq_to_postgres.py` does incremental sync to PostgreSQL (BQ → GCS Parquet → PostgreSQL, with `sync_state` tracking). `bq_to_es_projector/` indexes into Elasticsearch. PostgreSQL requires three materialized views (`perturb_seq_summary_perturbation`, `perturb_seq_summary_effect`, `perturb_seq_summary_dataset`) — see `prompts/postgres-summary.md` for creation SQL.

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
# Optional AI env vars: GEMINI_API_KEY or GOOGLE_CLOUD_PROJECT (Vertex AI), GEMINI_MODEL (default: gemini-2.5-flash)
cd be && uvicorn main:app --reload

# Frontend (requires PERTURBATION_CATALOGUE_BE env var pointing to backend URL, no trailing slash)
cd fe && python3 app.py
```

### Code formatting
```bash
black be/ fe/              # Manual formatting
pre-commit run --all-files # Run all pre-commit hooks (Black only, targeting python3.10)
```

### Testing
No automated test suite exists. `deploy/be_load_test.py` is a standalone async load-testing script (not pytest).

### dbt (data warehouse)
```bash
cd dwh/bq_dbt
dbt run                           # Run all models
dbt run --select +dataset_summary # Run with dependencies
```

## LLM Agent Rules (from prompts/README.md)

- Do NOT run `git`, `ruff`, `pre-commit`, or try to tidy up the repository.
- If tasked with writing code, simply do that and nothing else.

## Key Patterns

- **Modality mapping:** API URL modalities are lowercase (`perturb-seq`, `crispr-screen`, `mave`) but Elasticsearch values are mixed case (`Perturb-seq`, `CRISPR screen`, `MAVE`). The ES `lc_ascii` normalizer lowercases all keyword aggregation keys, so `main.py:perform_search` rebuilds original casing from `_source` documents before returning facets.
- **Two-phase modality search:** `/v1/{modality}/search` first queries PostgreSQL for dataset IDs (sorted by statistical significance), then queries Elasticsearch for metadata/facets restricted to those IDs. Python code re-sorts ES results to match the Postgres significance order.
- **Dynamic field mapping:** `be/dataset_metadata.json` drives the `DatasetMetadata` and `CommonModalitySearchParams` Pydantic models — both are constructed at import time from this JSON. Adding new ES dataset fields only requires updating this file.
- **Numeric range filter syntax:** API uses underscore syntax: `param=min_max` (range), `param=val_` (gte), `param=_val` (lte).
- **Search modes:** `/search` supports `search_mode=targets` (queries `target-summary` index) or `search_mode=datasets` (queries `dataset-summary` index). Uses cursor-based `search_after` pagination for deep results.
- **Async backend, sync frontend:** Backend uses `async/await` everywhere. Frontend uses synchronous `requests.get()` — correct since Dash is not async.
- **AI chat sessions are in-memory:** `_sessions` dict on the backend, not persistent. Container restarts or Cloud Run scale-down clears all sessions. MCP integration gracefully degrades if the `mcp` package is not installed.
- **Elasticsearch indexes:** `dataset-summary` (search), `target-summary` (targets), `landing-page-summary` (homepage stats).
- **Code style:** Black formatter enforced via pre-commit. No other linters.
- **Git branching:** `dev` is the default PR target branch. Feature branches use `name/feature-name` convention.

## Specifications

Detailed specs for backend and frontend live in `prompts/`:
- `prompts/be-specification.md` - Backend API implementation details
- `prompts/fe-specification.md` - Frontend implementation details
- `prompts/data-model.md` - Data model specification (API contract, filter/sort syntax, endpoint signatures)
- `prompts/postgres-summary.md` - PostgreSQL schema and materialized view creation SQL
