# BigQuery to Postgres Projector

## 1. Component Overview

Implement a Python script that loads data from a BigQuery table into a Google Cloud PostgreSQL instance:
- Save the script to `@postgres/bq_to_postgres.py`. The script must run on a Google Cloud VM. All script arguments must use the `--argument` form and be mandatory.
- Save its dependencies to `@postgres/requirements.txt`.
- Update the execution instructions in `@postgres/instructions.md` accordingly. Hypothetical or missing instructions are not acceptable. Provide exact commands to create the VM to run the script on. The VM must connect to Cloud VPC in order to be able to connect to the Cloud SQL instance by its private IP.

## 2. Script Arguments

### BigQuery Source
* `--bq-dataset`: The source dataset in BigQuery.
* `--bq-table`: The source table in BigQuery. The script must assume this table contains a timestamp column named `max_ingested_at`.
* `--bq-location`: The location of the BigQuery dataset.

### PostgreSQL Target
* `--pg-conn`: The connection URI for the target PostgreSQL database (e.g., `postgresql://user:password@host:port/dbname`).
* `--pg-table`: The name of the target table in PostgreSQL.

### Temporary Storage
* `--gcs-bucket`: The name of the Google Cloud Storage bucket to use for temporary data staging.

## 3. Data Transfer Workflow

To ensure high performance with large datasets, the script must follow this workflow:
1. **Export from BigQuery:** Export the required data from BigQuery to a new file (e.g., in CSV format) in the specified GCS bucket. The export URI must include a wildcard (`*`) to shard the output into multiple files.
2. **Load into PostgreSQL:** Use the PostgreSQL `COPY` command to efficiently bulk-load the data from the GCS file(s) into the target PostgreSQL table.
3. **Cleanup:** Remove the temporary file(s) from the GCS bucket after the operation completes.

## 4. Incremental Load and State Management

The script must perform incremental loads to avoid re-ingesting data. The synchronization state is maintained in a PostgreSQL table named `sync_state`.

### `sync_state` Table
If the `sync_state` table does not exist, the script must create it with the following schema:
- `table_name`: `TEXT`, `NOT NULL`, `PRIMARY KEY`. The name of the target table (`--pg-table`).
- `last_synced_at`: `TIMESTAMP WITHOUT TIME ZONE`. The latest `max_ingested_at` timestamp from the data successfully synced.

### Loading Logic
1. Query the `sync_state` table for the `last_synced_at` timestamp corresponding to `--pg-table`.
2. **Full Load:** If no entry exists for the table or `last_synced_at` is `NULL`, perform a full data load.
3. **Incremental Load:** Otherwise, export only those rows from BigQuery where `max_ingested_at` is greater than the `last_synced_at` timestamp.

## 5. Transactional Integrity

All data loading operations in PostgreSQL must be transactional to ensure data consistency.

### Full Load
1. Create a temporary staging table.
2. Load data from the GCS file(s) into this staging table.
3. In a single, atomic transaction: drop the original target table and rename the staging table to the target table name.

### Incremental Load
1. In a single transaction, insert the new data from the GCS file(s) into the main target table.

### State Update
After the data transaction is successfully committed, update the `sync_state` table. Set `last_synced_at` to the maximum `max_ingested_at` value from the records that were just loaded. This update must be an atomic `INSERT ... ON CONFLICT` operation to prevent race conditions.
