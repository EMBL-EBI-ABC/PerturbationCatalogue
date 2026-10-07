#!/usr/bin/env python3
"""Checks receipt validation, staging and scoped transactional replacement."""

from __future__ import annotations

import base64
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "dea-gsea"))

from load_results_to_bigquery import (
    DEA_SCHEMA,
    GcsDatasetFiles,
    GcsManifest,
    GcsParquet,
    GSEA_SCHEMA,
    EXPECTED_DATASET_COUNT,
    _load_stage,
    _verify_gcs_objects,
    apply_gcs_replacement,
    main,
    transaction_sql,
    validate_expected_datasets,
    validate_gcs_manifest,
)


def _receipt(dataset_ids: list[str]) -> dict:
    run_key = "0123456789abcdef0123456789abcdef"
    datasets = []
    for dataset_id in dataset_ids:
        item = {"dataset_id": dataset_id}
        for kind in ("dea", "gsea"):
            name = (
                f"perturb-seq-ingest/{run_key}/{dataset_id}/{dataset_id}.{kind}.parquet"
            )
            item[kind] = {
                "uri": f"gs://dev-tmp/{name}",
                "size_bytes": 123,
                "row_count": 1,
                "sha256": "a" * 64,
                "md5_hash": base64.b64encode(b"m" * 16).decode(),
                "crc32c": None,
                "generation": "987654321",
            }
        datasets.append(item)
    return {"run_key": run_key, "bucket": "dev-tmp", "datasets": datasets}


class BigQueryReplacementTest(unittest.TestCase):
    def test_receipt_is_bound_to_run_bucket_and_exact_object_names(self):
        receipt = _receipt(["first_2026"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "receipt.json"
            path.write_text(json.dumps(receipt))
            manifest = validate_gcs_manifest(path)
            self.assertEqual(manifest.datasets[0].dea.rows, 1)
            self.assertEqual(manifest.datasets[0].gsea.generation, 987654321)

            receipt["datasets"][0]["dea"]["uri"] = "gs://other-bucket/other-object"
            path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, "unexpected GCS URI"):
                validate_gcs_manifest(path)

    def test_selected_manifest_must_have_exactly_twenty_unique_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ids.json"
            entries = [{"dataset_id": f"dataset_{index}"} for index in range(20)]
            path.write_text(json.dumps({"datasets": entries}))
            self.assertEqual(
                len(validate_expected_datasets(path)), EXPECTED_DATASET_COUNT
            )

            path.write_text(json.dumps({"datasets": entries[:-1]}))
            with self.assertRaisesRegex(ValueError, "exactly 20"):
                validate_expected_datasets(path)

            entries[-1] = entries[0]
            path.write_text(json.dumps({"datasets": entries}))
            with self.assertRaisesRegex(ValueError, "must be unique"):
                validate_expected_datasets(path)

    def test_receipt_must_match_selected_twenty_ids_before_apply(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ids = [f"dataset_{index}" for index in range(EXPECTED_DATASET_COUNT)]
            receipt_path = root / "receipt.json"
            receipt_path.write_text(json.dumps(_receipt(ids)))
            selected_path = root / "selected.json"
            selected_ids = ids[:-1] + ["different_dataset"]
            selected_path.write_text(
                json.dumps(
                    {"datasets": [{"dataset_id": value} for value in selected_ids]}
                )
            )
            with self.assertRaisesRegex(ValueError, "do not match"):
                main(
                    [
                        "--gcs-manifest",
                        str(receipt_path),
                        "--id-manifest",
                        str(selected_path),
                    ]
                )

    def test_apply_rejects_a_receipt_outside_the_twenty_dataset_scope(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "receipt.json"
            path.write_text(json.dumps(_receipt(["dataset_0"])))
            manifest = validate_gcs_manifest(path)
            env = {
                "GCLOUD_PROJECT": "dev-project",
                "BQ_LOCATION": "europe-west2",
                "CLOUD_TMP_BUCKET": "dev-tmp",
            }
            with (
                patch.dict("os.environ", env),
                patch("load_results_to_bigquery._run_cli") as cli,
            ):
                with self.assertRaisesRegex(ValueError, "exactly 20 unique datasets"):
                    apply_gcs_replacement(manifest)
            cli.assert_not_called()

    def test_transaction_deletes_and_inserts_only_selected_ids_in_both_tables(self):
        sql = transaction_sql("dev-project", "ingest_run_dea", "ingest_run_gsea")
        self.assertEqual(sql.count("DELETE FROM"), 2)
        self.assertEqual(sql.count("INSERT INTO"), 2)
        self.assertEqual(sql.count("WHERE dataset_id IN UNNEST(@dataset_ids)"), 4)
        self.assertIn("BEGIN TRANSACTION", sql)
        self.assertIn("COMMIT TRANSACTION", sql)
        self.assertNotIn("TRUNCATE", sql)

    def test_cluster_receipt_objects_are_checked_by_generation_and_metadata(self):
        manifest_data = _receipt(["first_2026"])
        dataset = manifest_data["datasets"][0]
        files = {}
        for kind in ("dea", "gsea"):
            value = dataset[kind]
            name = value["uri"].split("/", 3)[-1]
            files[kind] = GcsParquet(
                uri=value["uri"],
                object_name=name,
                size_bytes=value["size_bytes"],
                rows=value["row_count"],
                sha256=value["sha256"],
                md5_hash=value["md5_hash"],
                crc32c=value["crc32c"],
                generation=int(value["generation"]),
            )
        manifest = GcsManifest(
            manifest_data["run_key"],
            "dev-tmp",
            [GcsDatasetFiles("first_2026", files["dea"], files["gsea"])],
        )
        responses = []
        for kind in ("dea", "gsea"):
            responses.append(
                json.dumps(
                    {
                        "generation": "987654321",
                        "size": "123",
                        "md5Hash": None,
                        "crc32c": None,
                        "metadata": {
                            "run_key": manifest.run_key,
                            "dataset_id": "first_2026",
                            "kind": kind,
                            "sha256": "a" * 64,
                            "row_count": "1",
                        },
                    }
                )
            )
        with patch("load_results_to_bigquery._run_cli", side_effect=responses):
            self.assertEqual(
                _verify_gcs_objects(manifest, "dev-project"),
                {item.object_name: item.generation for item in files.values()},
            )

    def test_load_stage_uses_parquet_list_inference_and_expires_table(self):
        files = [
            GcsParquet("gs://dev-tmp/one", "one", 10, 2, "a" * 64, "", None, 1),
            GcsParquet("gs://dev-tmp/two", "two", 10, 3, "b" * 64, "", None, 2),
        ]
        table_json = json.dumps({"schema": {"fields": DEA_SCHEMA}, "numRows": "5"})
        commands = []

        def fake_cli(args):
            commands.append(args)
            return table_json if args[0] == "bq" and args[2] == "show" else ""

        with patch("load_results_to_bigquery._run_cli", side_effect=fake_cli):
            stage = _load_stage(
                "dev-project", "europe-west2", "0" * 32, "dea", DEA_SCHEMA, files
            )

        self.assertEqual(stage, "dev-project.perturb_seq.ingest_" + "0" * 32 + "_dea")
        self.assertIn("--parquet_enable_list_inference=true", commands[0])
        self.assertIn("--replace=true", commands[0])
        self.assertEqual(commands[0][-1], "gs://dev-tmp/one,gs://dev-tmp/two")
        self.assertIn("--expiration=172800", commands[1])
        self.assertEqual(len(GSEA_SCHEMA), 13)


if __name__ == "__main__":
    unittest.main()
