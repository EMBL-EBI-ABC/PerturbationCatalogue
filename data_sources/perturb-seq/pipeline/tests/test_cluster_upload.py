#!/usr/bin/env python3
"""Checks cluster Parquet validation and the gcloud storage upload contract."""

from __future__ import annotations

import base64
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pyarrow as pa
import pyarrow.parquet as pq

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "dea-gsea"))

import cluster_upload
from io_schemas import DEA_SCHEMA


class ClusterUploadTest(unittest.TestCase):
    def test_upload_validates_data_and_uses_conditional_cli_upload(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset_id = "example_2026"
            source = root / f"{dataset_id}.dea.parquet"
            row = {field.name: None for field in DEA_SCHEMA}
            row["dataset_id"] = dataset_id
            pq.write_table(pa.Table.from_pylist([row], schema=DEA_SCHEMA), source)
            content = source.read_bytes()
            md5_hash = base64.b64encode(hashlib.md5(content).digest()).decode()
            sha256_hash = hashlib.sha256(content).hexdigest()
            run_key = "0123456789abcdef0123456789abcdef"
            object_name = (
                f"perturb-seq-ingest/{run_key}/{dataset_id}/{dataset_id}.dea.parquet"
            )
            metadata = {
                "run_key": run_key,
                "dataset_id": dataset_id,
                "kind": "dea",
                "sha256": sha256_hash,
                "row_count": "1",
            }
            described = {
                "name": object_name,
                "bucket": "dev-tmp",
                "generation": "987654321",
                "size": str(len(content)),
                "md5Hash": md5_hash,
                "crc32c": None,
                "metadata": metadata,
            }
            commands = []

            def fake_cli(args):
                commands.append(args)
                return json.dumps(described) if "describe" in args else ""

            source_entry = {
                "dea_parquet": str(source),
                "dea_size_bytes": len(content),
                "dea_rows": 1,
                "dea_md5_hash": md5_hash,
                "dea_sha256": sha256_hash,
                "dea_stat": {
                    "dev": source.stat().st_dev,
                    "ino": source.stat().st_ino,
                    "mtime_ns": source.stat().st_mtime_ns,
                },
            }
            with (
                patch.object(cluster_upload, "RESULTS_ROOT", str(root)),
                patch.object(cluster_upload, "_run_cli", side_effect=fake_cli),
            ):
                receipt = cluster_upload._upload_one(
                    run_key,
                    "dev-tmp",
                    "dev-project",
                    dataset_id,
                    "dea",
                    source_entry,
                )

            cp = commands[0]
            self.assertEqual(cp[:3], ["gcloud", "storage", "cp"])
            self.assertIn("--if-generation-match=0", cp)
            self.assertIn(f"--content-md5={md5_hash}", cp)
            self.assertIn(
                f"--custom-metadata="
                + ",".join(f"{k}={v}" for k, v in metadata.items()),
                cp,
            )
            self.assertEqual(cp[-1], f"gs://dev-tmp/{object_name}")
            self.assertEqual(receipt["row_count"], 1)
            self.assertEqual(receipt["generation"], "987654321")
            self.assertEqual(receipt["sha256"], sha256_hash)

    def test_wrong_dataset_rows_are_rejected_before_cloud_access(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "expected.dea.parquet"
            row = {field.name: None for field in DEA_SCHEMA}
            row["dataset_id"] = "other"
            pq.write_table(pa.Table.from_pylist([row], schema=DEA_SCHEMA), source)
            source_entry = {
                "dea_parquet": str(source),
                "dea_size_bytes": source.stat().st_size,
            }
            schema = [
                [field.name, str(field.type), field.nullable] for field in DEA_SCHEMA
            ]
            with (
                patch.object(cluster_upload, "RESULTS_ROOT", str(root)),
                patch.object(cluster_upload, "_run_cli") as cli,
            ):
                with self.assertRaisesRegex(
                    cluster_upload.TransferError, "Wrong dataset_id"
                ):
                    cluster_upload._verify_source(
                        str(source),
                        "expected",
                        "dea",
                        schema,
                    )
            cli.assert_not_called()


if __name__ == "__main__":
    unittest.main()
