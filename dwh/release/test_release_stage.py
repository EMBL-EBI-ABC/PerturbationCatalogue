import json
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from google.cloud import bigquery

import release_stage
from release_stage import _parse_dataset_ids, _stage, create_staging


class _Job:
    def result(self):
        return None


class _Client:
    def __init__(self):
        self.query_history = []

    def query(self, query, **kwargs):
        self.query_history.append(query)
        self.query_text = query
        self.query_kwargs = kwargs
        return _Job()

    def get_table(self, _):
        return SimpleNamespace(expires=None)

    def update_table(self, *_):
        return None


class _Rows:
    def __init__(self, rows=()):
        self.rows = rows

    def __iter__(self):
        return iter(self.rows)


class _StagingClient(_Client):
    def query(self, query, **kwargs):
        super().query(query, **kwargs)
        if query.startswith("CREATE OR REPLACE TABLE"):
            return _Job()
        if "release_run_perturb_seq_data" in query:
            return _JobWithRows(
                [
                    SimpleNamespace(dataset_id="pseq-1"),
                    SimpleNamespace(dataset_id="pseq-empty"),
                ]
            )
        return _JobWithRows()


class _JobWithRows(_Job):
    def __init__(self, rows=()):
        self.rows = _Rows(rows)

    def result(self):
        return self.rows


class _Blob:
    def upload_from_string(self, content, **_):
        self.content = content


class _StorageClient:
    def __init__(self):
        self.manifest_blob = _Blob()

    def bucket(self, _):
        return self

    def blob(self, _):
        return self.manifest_blob


class TestReleaseStage(unittest.TestCase):
    def test_dataset_ids_are_trimmed_nonempty_and_unique(self):
        self.assertEqual(_parse_dataset_ids(" ds1, ds2 "), ["ds1", "ds2"])
        with self.assertRaises(ValueError):
            _parse_dataset_ids(" , ")
        with self.assertRaises(ValueError):
            _parse_dataset_ids("ds1, ds1")

    def test_stage_uses_array_parameter_for_dataset_filter(self):
        client = _Client()
        _stage(
            client,
            "project",
            "dataset",
            "US",
            "stage",
            "SELECT dataset_id FROM `project.dataset.source`",
            ["ds1", "ds2"],
        )

        self.assertIn("dataset_id IN UNNEST(@dataset_ids)", client.query_text)
        parameter = client.query_kwargs["job_config"].query_parameters[0]
        self.assertIsInstance(parameter, bigquery.ArrayQueryParameter)
        self.assertEqual(parameter.name, "dataset_ids")
        self.assertEqual(parameter.values, ["ds1", "ds2"])

    def test_stage_keeps_full_release_query_unfiltered_by_default(self):
        client = _Client()
        _stage(
            client,
            "project",
            "dataset",
            "US",
            "stage",
            "SELECT dataset_id FROM `project.dataset.source`",
        )

        self.assertNotIn("UNNEST(@dataset_ids)", client.query_text)
        self.assertEqual(client.query_kwargs, {"location": "US"})

    def test_invalid_filter_is_rejected_before_cloud_access(self):
        with self.assertRaisesRegex(ValueError, "unique"):
            create_staging("p", "d", "US", "bucket", "run", ["ds1", "ds1"])

    def test_scoped_manifest_contains_only_requested_dataset_rows(self):
        client = _StagingClient()
        storage_client = _StorageClient()
        with (
            patch.object(release_stage.bigquery, "Client", return_value=client),
            patch.object(release_stage.storage, "Client", return_value=storage_client),
        ):
            task_count = create_staging(
                "project",
                "dataset",
                "US",
                "bucket",
                "run",
                ["pseq-1", "pseq-empty"],
            )

        manifest = json.loads(storage_client.manifest_blob.content)
        self.assertEqual(task_count, 2)
        self.assertEqual(
            [(item["modality"], item["dataset_id"]) for item in manifest["items"]],
            [("perturb-seq", "pseq-1"), ("perturb-seq", "pseq-empty")],
        )
        self.assertIn("gsea_table", manifest["items"][0])
        query = next(
            query
            for query in client.query_history
            if "release_run_perturb_seq_data" in query
            and query.startswith("SELECT DISTINCT")
        )
        self.assertIn("UNION DISTINCT", query)
        self.assertIn("release_run_dataset_metadata", query)
        self.assertIn("release_run_perturb_seq_gsea", query)


if __name__ == "__main__":
    unittest.main()
