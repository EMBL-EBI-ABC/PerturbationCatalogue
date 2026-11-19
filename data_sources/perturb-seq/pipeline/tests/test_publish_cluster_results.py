#!/usr/bin/env python3
"""Check the selected-ID BigQuery replacement commands."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "dea-gsea"))

import publish_cluster_results as publisher


class PublishClusterResultsTest(unittest.TestCase):
    def test_upload_delete_load_cleanup_is_scoped_to_manifest_ids(self):
        ids = [f"dataset_{index}" for index in range(20)]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "ids.json"
            manifest.write_text(
                json.dumps({"datasets": [{"dataset_id": item} for item in ids]})
            )
            self.assertEqual(publisher.read_ids(manifest), ids)
            for dataset_id in ids:
                result_dir = root / "pool" / "dea_gsea"
                result_dir.mkdir(parents=True, exist_ok=True)
                for kind in publisher.TABLES:
                    (result_dir / f"{dataset_id}.{kind}.parquet").touch()

            commands = []

            def run(command, check):
                self.assertTrue(check)
                commands.append(command)

            with (
                patch.object(publisher, "RESULTS_ROOT", root),
                patch.object(publisher.subprocess, "run", side_effect=run),
            ):
                prefix = publisher.publish(
                    ids, "dev-project", "europe-west2", "dev-tmp"
                )
                with self.assertRaisesRegex(ValueError, "Refusing non-development"):
                    publisher.publish(
                        ids, "production-project", "europe-west2", "dev-tmp"
                    )

        self.assertEqual(len(commands), 5)
        self.assertEqual(commands[0][:3], ["gcloud", "storage", "cp"])
        self.assertEqual(len(commands[0][3:-1]), 40)
        self.assertEqual(commands[0][-1], prefix + "/")
        delete_sql = commands[1][-1]
        self.assertEqual(delete_sql.count("DELETE FROM"), 2)
        self.assertTrue(all(f"'{dataset_id}'" in delete_sql for dataset_id in ids))
        self.assertNotIn("TRUNCATE", delete_sql)
        self.assertIn("pertpy_dea", commands[2][-2])
        self.assertIn("*.dea.parquet", commands[2][-1])
        self.assertIn("pertpy_gsea", commands[3][-2])
        self.assertIn("*.gsea.parquet", commands[3][-1])
        self.assertEqual(commands[4][:4], ["gcloud", "storage", "rm", "--recursive"])


if __name__ == "__main__":
    unittest.main()
