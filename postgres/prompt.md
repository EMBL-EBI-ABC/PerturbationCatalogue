# BigQuery to Postgres projector

## 1. Component overview

Implement a highly efficient Python script which will load data from BigQuery into a Google Cloud PostgreSQL instance and save it to @postgres/bq_to_postgres.py. Save the necessary requirements to @postgres/requirements.txt and instructions to run it to @postgres/instructions.md.

All operations must be performed in Google Cloud, never on a local machine. Take that into account when writing the instructions. The script must also be as efficient as possible to be able to move large amounts of data.

All arguments of the script must use the two dash form (--argument) and be mandatory, unless specified otherwise.

## 2. Input data in BigQuery
Arguments:
* bq_dataset: source dataset in BigQuery
* bq_table: source table in BigQuery

## 3. Output data in PostgreSQL
Arguments:
* pg_conn: connection information for the target database
* pg_table: target table in PostgreSQL

## 4. Incremental load

To avoid moving data which was already ingested, the script must always consider the synchronisation state of BigQuery and PostgreSQL. It is stored in the "sync_state" table of the PostgreSQL database. Instructions to work with the table are as follows.

First, if the `sync_state` table is missing, create it. It only has two fields:
- `table_name`: text, not null. It refers to the *target* table name in PostgreSQL, corresponding to pg_table above.
- `last_synced`: timestamp without timezone. It refers to the highest timestamp of the "max_ingested_at" field in the corresponding BigQuery data.

Then:
- If the target table is not mentioned in the `sync_state` table or the timestamp is empty, then the target table must be dropped and fully ingested from BigQuery.
- Otherwise, only those entries should be ingested from BigQuery such that last_synced (of the sync_state) < max_ingested_at (of the particular entry from BigQuery).

After the import is completed, the `sync_state` table must be updated.
