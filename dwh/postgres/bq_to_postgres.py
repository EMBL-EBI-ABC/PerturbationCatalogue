import argparse
import logging
import io
import sys
import uuid
import threading
import queue

from google.cloud import bigquery
import psycopg2
from psycopg2 import sql
from tqdm import tqdm
import pyarrow.csv as pacsv

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

TABLES_TO_SYNC = [
    "crispr_data",
    "mave_data",
    "perturb_seq_dea",
    "perturb_seq_gsea",
]


class BQToPGStream(io.RawIOBase):
    """
    High-performance stream that buffers Arrow batches from BQ (Producer)
    and serves them in TSV format to Postgres COPY (Consumer).
    """

    def __init__(self, bq_client, temp_table_id, bq_schema, pbar=None):
        self.bq_client = bq_client
        self.temp_table_id = temp_table_id
        self.bq_schema = bq_schema
        self.pbar = pbar

        self.queue = queue.Queue(maxsize=10)  # Buffer up to 10 batches
        self.buffer = b""
        self.finished = False
        self.error = None

        # Start the background producer thread
        self.thread = threading.Thread(target=self._producer)
        self.thread.daemon = True
        self.thread.start()

    def readable(self):
        return True

    def _producer(self):
        """Background thread to fetch data from BigQuery and serialize to TSV."""
        try:
            # Fetch data in Arrow record batches
            row_iterator = self.bq_client.list_rows(self.temp_table_id)
            arrow_batches = row_iterator.to_arrow_iterable()

            write_options = pacsv.WriteOptions(
                include_header=False, delimiter="\t", quoting_style="none"
            )

            for batch in arrow_batches:
                # Serialize the entire batch to TSV in-memory using C-based pyarrow
                out = io.BytesIO()
                pacsv.write_csv(batch, out, write_options=write_options)
                tsv_chunk = out.getvalue()

                self.queue.put(tsv_chunk)

                if self.pbar:
                    self.pbar.update(batch.num_rows)

            self.queue.put(None)  # Sentinel for completion
        except Exception as e:
            self.error = e
            self.queue.put(None)

    def read(self, n=-1):
        if self.error:
            raise self.error

        while not self.buffer:
            if self.finished:
                return b""

            chunk = self.queue.get()
            if chunk is None:
                self.finished = True
                continue

            self.buffer = chunk

        if n == -1 or n >= len(self.buffer):
            result = self.buffer
            self.buffer = b""
        else:
            result = self.buffer[:n]
            self.buffer = self.buffer[n:]
        return result


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
    """High-performance stream from BQ to PG using Arrow, threading, and bulk serialization."""
    temp_table_id = f"{bq_client.project}.{bq_dataset}.temp_sync_{uuid.uuid4().hex}"

    job_config = bigquery.QueryJobConfig(
        destination=temp_table_id,
        write_disposition="WRITE_TRUNCATE",
    )

    query = f"SELECT * FROM `{bq_dataset}.{bq_table}` WHERE dataset_id = '{dataset_id}'"

    # query_job = bq_client.query(query, job_config=job_config, location=bq_location)
    # query_job.result() # Wait for completion

    # Fetching the query results to a temp table is still needed for huge results scalability
    query_job = bq_client.query(query, job_config=job_config, location=bq_location)
    query_job.result()

    with tqdm(
        total=total_rows, desc=f"    Syncing {dataset_id}", leave=False, unit="rows"
    ) as pbar:
        # Use our high-performance asynchronous stream
        stream = BQToPGStream(bq_client, temp_table_id, bq_schema, pbar=pbar)
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

                            # High-performance streaming
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
