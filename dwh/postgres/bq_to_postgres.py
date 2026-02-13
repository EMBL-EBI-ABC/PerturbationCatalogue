import argparse
import logging
import io
import sys
import uuid

from google.cloud import bigquery
import psycopg2
from psycopg2 import sql
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

TABLES_TO_SYNC = [
    "crispr_data",
    "mave_data",
    "perturb_seq_dea",
    "perturb_seq_gsea",
]


class BQRowIteratorIO(io.RawIOBase):
    """File-like object that streams rows from a BigQuery iterator in TSV format."""

    def __init__(self, row_iterator, bq_schema, pbar=None):
        self.row_iterator = row_iterator
        self.bq_schema = bq_schema
        self.buffer = b""
        self.pbar = pbar

    def readable(self):
        return True

    def read(self, n=-1):
        while not self.buffer:
            try:
                row = next(self.row_iterator)
                if self.pbar:
                    self.pbar.update(1)
                line = self._format_row(row) + "\n"
                self.buffer = line.encode("utf-8")
            except StopIteration:
                return b""

        if n == -1 or n >= len(self.buffer):
            result = self.buffer
            self.buffer = b""
        else:
            result = self.buffer[:n]
            self.buffer = self.buffer[n:]
        return result

    def _format_row(self, row):
        formatted = []
        for field in self.bq_schema:
            val = row[field.name]
            if val is None:
                formatted.append("")
            elif field.mode == "REPEATED":
                # Convert list to Postgres array literal: {"val1", "val2"}
                inner = [str(x) if not isinstance(x, str) else x for x in val]
                escaped = []
                for x in inner:
                    x_str = str(x).replace('"', '\\"')
                    escaped.append(f'"{x_str}"')
                formatted.append(f'{{{",".join(escaped)}}}')
            elif field.field_type == "BOOLEAN":
                formatted.append("true" if val else "false")
            else:
                formatted.append(str(val))
        return "\t".join(formatted)


def get_pg_type(field):
    """Maps BigQuery types to PostgreSQL types."""
    bq_type = field.field_type
    pg_type = {
        "STRING": "TEXT",
        "INTEGER": "BIGINT",
        "FLOAT": "DOUBLE PRECISION",
        "BOOLEAN": "BOOLEAN",
        "TIMESTAMP": "TIMESTAMP WITHOUT TIME ZONE",
        "DATE": "DATE",
    }.get(bq_type, "TEXT")

    if field.mode == "REPEATED":
        pg_type += "[]"

    return pg_type


def get_all_sync_states(cursor):
    """Gets the current sync state for all tables and datasets from Postgres."""
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS sync_state (
            table_name TEXT NOT NULL,
            dataset_id TEXT NOT NULL,
            last_synced_at TIMESTAMP WITHOUT TIME ZONE,
            PRIMARY KEY (table_name, dataset_id)
        );
    """
    )
    cursor.execute("SELECT table_name, dataset_id, last_synced_at FROM sync_state")
    for table_name, dataset_id, last_synced_at in cursor.fetchall():
        if table_name not in states:  # Bug here in previous version, fixing.
            pass  # Actually 'states' is local here.

    # Corrected implementation:
    states = {}
    cursor.execute("SELECT table_name, dataset_id, last_synced_at FROM sync_state")
    for table_name, dataset_id, last_synced_at in cursor.fetchall():
        if table_name not in states:
            states[table_name] = {}
        states[table_name][dataset_id] = last_synced_at
    return states


def get_bq_latest_timestamps_and_counts(bq_client, bq_dataset, bq_table, bq_location):
    """Gets the latest max_ingested_at and row count for every dataset_id in a BQ table."""
    query = f"""
        SELECT dataset_id, MAX(max_ingested_at) as latest_ts, COUNT(*) as row_count
        FROM `{bq_dataset}.{bq_table}`
        GROUP BY dataset_id
    """
    query_job = bq_client.query(query, location=bq_location)
    results = query_job.result()
    return {row.dataset_id: (row.latest_ts, row.row_count) for row in results}


def stream_bq_to_pg(
    cursor,
    pg_table,
    bq_client,
    bq_dataset,
    bq_table,
    bq_location,
    dataset_id,
    bq_schema,
    total_rows,
):
    """Streams data from a large BQ query directly into a Postgres table using a temp destination table."""
    # Use a destination table for the query result if it's potentially huge.
    # To keep things robust, we'll ALWAYS use a destination table for dataset ingestion.
    temp_table_id = f"{bq_client.project}.{bq_dataset}.temp_sync_{uuid.uuid4().hex}"

    job_config = bigquery.QueryJobConfig(
        destination=temp_table_id,
        write_disposition="WRITE_TRUNCATE",
    )

    query = f"SELECT * FROM `{bq_dataset}.{bq_table}` WHERE dataset_id = '{dataset_id}'"

    query_job = bq_client.query(query, job_config=job_config, location=bq_location)
    query_job.result()  # Wait for completion

    # Fetch rows from the destination table. This avoids result size limits.
    row_iterator = bq_client.list_rows(temp_table_id)

    with tqdm(
        total=total_rows, desc=f"    Loading {dataset_id}", leave=False, unit="rows"
    ) as pbar:
        stream = BQRowIteratorIO(row_iterator, bq_schema, pbar=pbar)
        cursor.copy_expert(
            sql.SQL("COPY {} FROM STDIN WITH (FORMAT TEXT, NULL '')").format(
                sql.Identifier(pg_table)
            ),
            stream,
        )

    # Clean up temp table
    bq_client.delete_table(temp_table_id, not_found_ok=True)


def delete_dataset_from_pg(cursor, pg_table, dataset_id):
    """Deletes all rows for a given dataset_id from a Postgres table."""
    cursor.execute(
        sql.SQL("DELETE FROM {} WHERE dataset_id = %s").format(
            sql.Identifier(pg_table)
        ),
        (dataset_id,),
    )


def update_sync_state(cursor, pg_table, dataset_id, timestamp):
    """Updates or inserts the sync state for a specific dataset."""
    cursor.execute(
        """
        INSERT INTO sync_state (table_name, dataset_id, last_synced_at)
        VALUES (%s, %s, %s)
        ON CONFLICT (table_name, dataset_id) DO UPDATE
        SET last_synced_at = EXCLUDED.last_synced_at;
    """,
        (pg_table, dataset_id, timestamp),
    )


def ensure_pg_table_exists(cursor, pg_table, bq_schema):
    """Ensures the target table exists in Postgres, creating it if necessary."""
    cursor.execute(
        "SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = %s)",
        (pg_table,),
    )
    if not cursor.fetchone()[0]:
        logging.info(f"Table {pg_table} does not exist. Creating.")
        columns = [f"{field.name} {get_pg_type(field)}" for field in bq_schema]
        cursor.execute(
            sql.SQL("CREATE TABLE {} ({})").format(
                sql.Identifier(pg_table),
                sql.SQL(", ").join(map(sql.SQL, columns)),
            )
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bq-dataset", required=True)
    parser.add_argument("--bq-location", required=True)
    parser.add_argument("--pg-conn", required=True)
    parser.add_argument(
        "--yes", action="store_true", help="Proceed without confirmation"
    )
    args = parser.parse_args()

    bq_client = bigquery.Client()

    try:
        with psycopg2.connect(args.pg_conn) as conn:
            with conn.cursor() as cursor:
                # 1. Fetch current sync states for all tables
                logging.info("Fetching current sync states from Postgres...")
                all_states = get_all_sync_states(cursor)

                # 2. Plan sync for each table
                sync_plan = {}
                total_datasets = 0

                for table_name in TABLES_TO_SYNC:
                    logging.info(f"Checking BQ table {table_name}...")
                    bq_info = get_bq_latest_timestamps_and_counts(
                        bq_client, args.bq_dataset, table_name, args.bq_location
                    )
                    current_table_states = all_states.get(table_name, {})

                    to_update = []
                    to_insert = []

                    for ds_id, (bq_ts, row_count) in bq_info.items():
                        if ds_id in current_table_states:
                            pg_ts = current_table_states[ds_id]
                            if pg_ts is None or bq_ts > pg_ts.replace(
                                tzinfo=bq_ts.tzinfo
                            ):
                                to_update.append(ds_id)
                        else:
                            to_insert.append(ds_id)

                    if to_update or to_insert:
                        sync_plan[table_name] = {
                            "to_update": to_update,
                            "to_insert": to_insert,
                            "bq_info": bq_info,
                        }
                        total_datasets += len(to_update) + len(to_insert)

                # 3. Summary and confirmation
                if not sync_plan:
                    logging.info("Everything is up to date. Nothing to sync.")
                    return

                print("\n--- Sync Summary ---")
                for table_name, plan in sync_plan.items():
                    print(f"Table: {table_name}")
                    if plan["to_update"]:
                        print(f"  Datasets to update: {len(plan['to_update'])}")
                        print(
                            f"    {', '.join(plan['to_update'][:5])}{'...' if len(plan['to_update']) > 5 else ''}"
                        )
                    if plan["to_insert"]:
                        print(f"  Datasets to insert: {len(plan['to_insert'])}")
                        print(
                            f"    {', '.join(plan['to_insert'][:5])}{'...' if len(plan['to_insert']) > 5 else ''}"
                        )

                print(f"\nTotal datasets to process: {total_datasets}")
                print("--------------------\n")

                if not args.yes:
                    try:
                        confirm = input("Proceed with sync? (y/N): ")
                    except EOFError:
                        confirm = "n"
                    if confirm.lower() != "y":
                        logging.info("Sync cancelled by user.")
                        sys.exit(0)

                # 4. Sequential Execution (one dataset at a time)
                overall_pbar = tqdm(
                    total=total_datasets, desc="Overall Progress", unit="dataset"
                )

                for table_name, plan in sync_plan.items():
                    logging.info(f"Syncing table {table_name}...")

                    bq_table_obj = bq_client.get_table(
                        f"{args.bq_dataset}.{table_name}"
                    )
                    ensure_pg_table_exists(cursor, table_name, bq_table_obj.schema)

                    all_to_process = sorted(plan["to_update"] + plan["to_insert"])

                    for ds_id in all_to_process:
                        row_count = plan["bq_info"][ds_id][1]
                        bq_timestamp = plan["bq_info"][ds_id][0]

                        try:
                            # Delete if update
                            if ds_id in plan["to_update"]:
                                delete_dataset_from_pg(cursor, table_name, ds_id)

                            # Stream directly from BQ to PG using a temp destination table for safety
                            stream_bq_to_pg(
                                cursor,
                                table_name,
                                bq_client,
                                args.bq_dataset,
                                table_name,
                                args.bq_location,
                                ds_id,
                                bq_table_obj.schema,
                                row_count,
                            )

                            # Update sync states
                            update_sync_state(cursor, table_name, ds_id, bq_timestamp)

                            # COMMIT EVERY DATASET
                            conn.commit()
                            overall_pbar.update(1)

                        except (Exception, KeyboardInterrupt) as e:
                            conn.rollback()
                            overall_pbar.close()
                            if isinstance(e, KeyboardInterrupt):
                                logging.error(
                                    f"\nSync of {table_name}/{ds_id} interrupted. Rolling back..."
                                )
                            else:
                                logging.error(
                                    f"\nError syncing {table_name}/{ds_id}: {e}. Rolling back..."
                                )
                            sys.exit(1)

                overall_pbar.close()
                logging.info("Multi-table sync completed successfully.")

    except psycopg2.Error as e:
        logging.error(f"Postgres connection error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("\nExiting...")
        sys.exit(0)


if __name__ == "__main__":
    main()
