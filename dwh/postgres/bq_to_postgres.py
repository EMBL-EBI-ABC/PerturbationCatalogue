import argparse
import logging
import io
import sys

from google.cloud import bigquery
import psycopg2
from psycopg2 import sql
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


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
    states = {}
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


def stream_bq_to_pg(cursor, pg_table, bq_client, bq_query, bq_schema, total_rows):
    """Streams data from a BQ query directly into a Postgres table using COPY."""
    query_job = bq_client.query(bq_query)
    row_iterator = query_job.result()

    with tqdm(
        total=total_rows, desc=f"  Loading rows to {pg_table}", leave=False
    ) as pbar:
        stream = BQRowIteratorIO(row_iterator, bq_schema, pbar=pbar)
        cursor.copy_expert(
            sql.SQL("COPY {} FROM STDIN WITH (FORMAT TEXT, NULL '')").format(
                sql.Identifier(pg_table)
            ),
            stream,
        )


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
    parser.add_argument("--bq-table", required=True)
    parser.add_argument("--bq-location", required=True)
    parser.add_argument("--pg-conn", required=True)
    parser.add_argument("--pg-table", required=True)
    parser.add_argument(
        "--yes", action="store_true", help="Proceed without confirmation"
    )
    args = parser.parse_args()

    bq_client = bigquery.Client()

    try:
        with psycopg2.connect(args.pg_conn) as conn:
            with conn.cursor() as cursor:
                # 1. Fetch current sync states
                logging.info("Fetching current sync states from Postgres...")
                all_states = get_all_sync_states(cursor)
                current_table_states = all_states.get(args.pg_table, {})

                # 2. Fetch latest timestamps and counts from BQ
                logging.info(
                    f"Fetching latest dataset info from BQ table {args.bq_table}..."
                )
                bq_info = get_bq_latest_timestamps_and_counts(
                    bq_client, args.bq_dataset, args.bq_table, args.bq_location
                )

                # 3. Determine datasets to process
                to_update = []  # Exists in PG but timestamp in BQ is newer
                to_insert = []  # Does not exist in PG

                for ds_id, (bq_ts, _) in bq_info.items():
                    if ds_id in current_table_states:
                        pg_ts = current_table_states[ds_id]
                        if pg_ts is None or bq_ts > pg_ts.replace(tzinfo=bq_ts.tzinfo):
                            to_update.append(ds_id)
                    else:
                        to_insert.append(ds_id)

                # 4. Summary and confirmation
                total_to_process = len(to_update) + len(to_insert)
                if total_to_process == 0:
                    logging.info("Everything is up to date. Nothing to sync.")
                    return

                print("\n--- Sync Summary ---")
                print(f"Table: {args.pg_table}")
                print(f"Datasets to update (delete + re-ingest): {len(to_update)}")
                if to_update:
                    print(
                        f"  {', '.join(to_update[:10])}{'...' if len(to_update) > 10 else ''}"
                    )
                print(f"Datasets to insert (new): {len(to_insert)}")
                if to_insert:
                    print(
                        f"  {', '.join(to_insert[:10])}{'...' if len(to_insert) > 10 else ''}"
                    )
                print(f"Total datasets to process: {total_to_process}")
                print("--------------------\n")

                if not args.yes:
                    try:
                        confirm = input("Proceed with sync? (y/N): ")
                    except EOFError:
                        confirm = "n"
                    if confirm.lower() != "y":
                        logging.info("Sync cancelled by user.")
                        sys.exit(0)

                # 5. Execution
                bq_table_obj = bq_client.get_table(f"{args.bq_dataset}.{args.bq_table}")
                ensure_pg_table_exists(cursor, args.pg_table, bq_table_obj.schema)

                all_to_process = sorted(to_update + to_insert)
                chunk_size = 5

                overall_pbar = tqdm(
                    total=len(all_to_process), desc="Synchronizing datasets"
                )

                for i in range(0, len(all_to_process), chunk_size):
                    chunk = all_to_process[i : i + chunk_size]
                    chunk_rows = sum([bq_info[ds_id][1] for ds_id in chunk])

                    try:
                        # Delete if update
                        for ds_id in chunk:
                            if ds_id in to_update:
                                delete_dataset_from_pg(cursor, args.pg_table, ds_id)

                        # Build query for the chunk
                        ds_ids_str = ", ".join([f"'{ds}'" for ds in chunk])
                        bq_query = f"""
                            SELECT *
                            FROM `{args.bq_dataset}.{args.bq_table}`
                            WHERE dataset_id IN ({ds_ids_str})
                        """

                        # Stream directly from BQ to PG
                        stream_bq_to_pg(
                            cursor,
                            args.pg_table,
                            bq_client,
                            bq_query,
                            bq_table_obj.schema,
                            chunk_rows,
                        )

                        # Update sync states
                        for ds_id in chunk:
                            update_sync_state(
                                cursor, args.pg_table, ds_id, bq_info[ds_id][0]
                            )

                        # COMMIT EACH CHUNK
                        conn.commit()
                        overall_pbar.update(len(chunk))

                    except (Exception, KeyboardInterrupt) as e:
                        conn.rollback()
                        overall_pbar.close()
                        if isinstance(e, KeyboardInterrupt):
                            logging.error(
                                "\nSync interrupted by user. Rolling back current chunk..."
                            )
                        else:
                            logging.error(f"\nError during sync: {e}. Rolling back...")
                        sys.exit(1)

                overall_pbar.close()
                logging.info("Sync completed successfully.")

    except psycopg2.Error as e:
        logging.error(f"Postgres connection error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("\nExiting...")
        sys.exit(0)


if __name__ == "__main__":
    main()
