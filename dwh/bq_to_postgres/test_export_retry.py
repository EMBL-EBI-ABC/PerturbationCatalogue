import io
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from google.api_core.exceptions import NotFound

import bq_to_postgres as sync


class Blob:
    def __init__(self, contents=None, missing=False):
        self.contents = contents
        self.missing = missing

    def download_to_file(self, target):
        if self.missing:
            raise NotFound("export shard disappeared")
        target.write(self.contents)

    def delete(self):
        pass


class Bucket:
    def __init__(self):
        self.blobs = {}

    def list_blobs(self, prefix):
        return self.blobs[prefix]


class ExportRetryTest(unittest.TestCase):
    def setUp(self):
        self.bucket = Bucket()
        self.synchronizer = sync.TableSynchronizer.__new__(sync.TableSynchronizer)
        self.synchronizer.bq_client = object()
        self.synchronizer.bq_dataset = "test"
        self.synchronizer.bq_location = "europe-west2"
        self.synchronizer.gcs_bucket = "tmp"
        self.synchronizer.gcs_client = SimpleNamespace(get_bucket=lambda _: self.bucket)

    def test_missing_shard_reexports_once_and_discards_partial_tsv(self):
        exports = []

        def export(*args):
            prefix = args[5]
            self.bucket.blobs[prefix] = (
                [Blob(b"partial"), Blob(missing=True)]
                if not exports
                else [Blob(b"fresh")]
            )
            exports.append(prefix)

        with (
            patch.object(sync, "export_dataset_to_gcs", side_effect=export),
            patch.object(sync.pq, "read_table", return_value=object()),
            patch.object(
                sync,
                "_prepare_table_for_copy",
                side_effect=[io.BytesIO(b"partial\n"), io.BytesIO(b"fresh\n")],
            ) as prepare,
        ):
            dataset_id, output = self.synchronizer._export_and_prepare_dataset(
                "perturb_seq_dea", "dataset_a", []
            )

        try:
            self.assertEqual(dataset_id, "dataset_a")
            self.assertEqual(output.read(), b"fresh\n")
            self.assertEqual(len(exports), 2)
            self.assertNotEqual(exports[0], exports[1])
            self.assertEqual(prepare.call_count, 2)
        finally:
            name = output.name
            output.close()
            sync.os.unlink(name)

    def test_second_missing_shard_fails_without_prepared_output(self):
        exports = []

        def export(*args):
            prefix = args[5]
            self.bucket.blobs[prefix] = [Blob(missing=True)]
            exports.append(prefix)

        with (
            patch.object(sync, "export_dataset_to_gcs", side_effect=export),
            patch.object(sync, "_prepare_table_for_copy") as prepare,
        ):
            with self.assertRaises(NotFound):
                self.synchronizer._export_and_prepare_dataset(
                    "perturb_seq_dea", "dataset_a", []
                )

        self.assertEqual(len(exports), 2)
        self.assertNotEqual(exports[0], exports[1])
        prepare.assert_not_called()


if __name__ == "__main__":
    unittest.main()
