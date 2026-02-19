FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    postgresql-client \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy and install Python dependencies
COPY dwh/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy dbt project
COPY dwh/bq_dbt/ bq_dbt/

# Copy bq_to_postgres
COPY dwh/bq_to_postgres/bq_to_postgres.py bq_to_postgres/bq_to_postgres.py

# Copy bq_to_elastic
COPY dwh/bq_to_elastic/ bq_to_elastic/

# Copy dataset_metadata.json (used by ES projector to generate dataset-summary mapping)
COPY be/dataset_metadata.json be/dataset_metadata.json

# Install dbt packages
RUN cd bq_dbt && dbt deps --profiles-dir .
