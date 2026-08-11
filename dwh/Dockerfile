FROM python:3.11-slim

# Install system dependencies
RUN apt-get update -qq && \
    apt-get install -y -qq --no-install-recommends \
    git \
    postgresql-client \
    > /dev/null && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -q -r requirements.txt

# Copy dbt project
COPY bq_dbt/ bq_dbt/

# Copy bq_to_postgres
COPY bq_to_postgres/bq_to_postgres.py bq_to_postgres/bq_to_postgres.py

# Copy release artifact generation
COPY release/release.py release/release.py
COPY release/release_stage.py release/release_stage.py
COPY release/release_worker.py release/release_worker.py
COPY release/metadata.py release/metadata.py

# Copy bq_to_elastic
COPY bq_to_elastic/ bq_to_elastic/

# Copy dataset_metadata.json (used by ES projector to generate dataset-summary mapping)
# This file is copied into dwh/ by trigger_pipeline.sh before build submission
COPY be/dataset_metadata.json be/dataset_metadata.json

# Install dbt packages
RUN cd bq_dbt && dbt deps --profiles-dir .
