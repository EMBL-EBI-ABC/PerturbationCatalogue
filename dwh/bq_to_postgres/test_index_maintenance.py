import io
import unittest
from datetime import datetime, timezone
from unittest.mock import patch

import bq_to_postgres as sync
from psycopg2 import sql


def _render(statement):
    if isinstance(statement, str):
        return statement
    if isinstance(statement, sql.SQL):
        return statement.string
    if isinstance(statement, sql.Identifier):
        return ".".join(f'"{part}"' for part in statement.strings)
    if isinstance(statement, sql.Composed):
        return "".join(_render(part) for part in statement.seq)
    return str(statement)


class RecordingCursor:
    def __init__(self, relation_kinds=None):
        self.queries = []
        self.relation_kinds = iter(relation_kinds or ["m"] * 3)
        self.result = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def execute(self, statement, params=None):
        query = _render(statement)
        self.queries.append(query)
        if query.startswith("SELECT relkind"):
            relkind = next(self.relation_kinds)
            self.result = None if relkind is None else (relkind,)

    def fetchone(self):
        result = self.result
        self.result = None
        return result


class RecordingConnection:
    def __init__(self, cursor):
        self.recording_cursor = cursor

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def cursor(self):
        return self.recording_cursor


class IndexMaintenanceTest(unittest.TestCase):
    def _synchronizer(self, rebuild_indexes):
        synchronizer = sync.TableSynchronizer.__new__(sync.TableSynchronizer)
        synchronizer.drop_and_recreate_indexes = rebuild_indexes
        synchronizer.pg_conn_str = "unused in unit test"
        synchronizer._get_bq_schema = lambda _: []
        synchronizer._export_and_prepare_dataset = lambda *_: (
            "dataset_a",
            io.BytesIO(b"row\n"),
        )
        synchronizer._load_tsv_to_pg = lambda *_: None
        return synchronizer

    def _sync_one_dataset(self, rebuild_indexes):
        cursor = RecordingCursor()
        connection = RecordingConnection(cursor)
        synchronizer = self._synchronizer(rebuild_indexes)
        timestamp = datetime(2026, 1, 1, tzinfo=timezone.utc)
        plan = {
            "to_insert": [],
            "to_update": ["dataset_a"],
            "bq_info": {"dataset_a": (timestamp, 1)},
        }
        with (
            patch.object(sync.psycopg2, "connect", return_value=connection),
            patch.object(sync, "ensure_pg_table_exists"),
            patch.object(sync, "update_sync_state"),
        ):
            synchronizer._sync_unified("perturb_seq_dea", plan)
        return cursor.queries

    def test_default_selective_sync_preserves_indexes_and_refreshes_each_mv_once(self):
        queries = self._sync_one_dataset(rebuild_indexes=False)

        self.assertFalse(any(query.startswith("DROP INDEX") for query in queries))
        self.assertFalse(
            any(query.startswith("DROP MATERIALIZED VIEW") for query in queries)
        )
        self.assertEqual(
            sum(query.startswith("CREATE INDEX IF NOT EXISTS") for query in queries),
            len(sync.INDEX_DEFINITIONS["perturb_seq_dea"]),
        )
        for mv_name in sync.TABLE_MATERIALIZED_VIEWS["perturb_seq_dea"]:
            self.assertEqual(
                sum(
                    query == f'REFRESH MATERIALIZED VIEW "{mv_name}"'
                    for query in queries
                ),
                1,
            )
            self.assertFalse(
                any(
                    query.lstrip().startswith(f"CREATE MATERIALIZED VIEW {mv_name}")
                    for query in queries
                )
            )

    def test_explicit_index_rebuild_drops_only_configured_indexes(self):
        queries = self._sync_one_dataset(rebuild_indexes=True)
        drops = [query for query in queries if query.startswith("DROP INDEX IF EXISTS")]

        self.assertEqual(len(drops), len(sync.INDEX_DEFINITIONS["perturb_seq_dea"]))
        for index_name, _ in sync.INDEX_DEFINITIONS["perturb_seq_dea"]:
            self.assertIn(f'DROP INDEX IF EXISTS "{index_name}"', drops)
        self.assertFalse(
            any(query.startswith("DROP MATERIALIZED VIEW") for query in queries)
        )

    def test_missing_mvs_are_created_with_idempotent_indexes(self):
        cursor = RecordingCursor([None, None, None])
        sync.ensure_defined_materialized_views(
            cursor, "perturb_seq_dea", refresh_existing=True
        )

        self.assertEqual(
            sum(
                query.lstrip().startswith("CREATE MATERIALIZED VIEW")
                for query in cursor.queries
            ),
            3,
        )
        self.assertEqual(
            sum(
                query.startswith("CREATE UNIQUE INDEX IF NOT EXISTS")
                for query in cursor.queries
            ),
            3,
        )
        self.assertFalse(
            any(
                query.startswith("REFRESH MATERIALIZED VIEW")
                for query in cursor.queries
            )
        )

    def test_unexpected_relation_kind_is_not_overwritten(self):
        cursor = RecordingCursor(["r"])

        with self.assertRaisesRegex(RuntimeError, "Expected .* materialized view"):
            sync.ensure_defined_materialized_views(
                cursor, "perturb_seq_dea", refresh_existing=True
            )

        self.assertEqual(len(cursor.queries), 1)

    def test_existing_mvs_are_not_refreshed_when_no_dataset_changed(self):
        cursor = RecordingCursor()
        sync.ensure_defined_materialized_views(
            cursor, "perturb_seq_dea", refresh_existing=False
        )

        self.assertFalse(
            any(
                query.startswith("REFRESH MATERIALIZED VIEW")
                for query in cursor.queries
            )
        )
        self.assertEqual(
            sum(
                query.startswith("CREATE UNIQUE INDEX IF NOT EXISTS")
                for query in cursor.queries
            ),
            3,
        )


if __name__ == "__main__":
    unittest.main()
