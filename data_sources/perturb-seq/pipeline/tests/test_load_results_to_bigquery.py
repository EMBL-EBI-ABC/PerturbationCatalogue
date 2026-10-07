#!/usr/bin/env python3
"""Focused local checks for scoped DEA/GSEA replacement inputs and SQL."""

from __future__ import annotations

import json
import base64
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pyarrow as pa
import pyarrow.parquet as pq

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "dea-gsea"))

from io_schemas import DEA_SCHEMA, GSEA_SCHEMA
from load_results_to_bigquery import (
    GcsDatasetFiles,
    GcsManifest,
    GcsParquet,
    _expected_bq_schema,
    _stage_and_replace,
    _verify_gcs_objects,
    apply_replacement,
    transaction_sql,
    validate_manifest,
    validate_gcs_manifest,
)


def _write(path: Path, schema: pa.Schema, dataset_ids: list[str]) -> None:
    rows = []
    for dataset_id in dataset_ids:
        row = {field.name: None for field in schema}
        row["dataset_id"] = dataset_id
        if "leading_edge" in schema.names:
            row["leading_edge"] = []
        rows.append(row)
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path)


class BigQueryReplacementTest(unittest.TestCase):
    def test_apply_rejects_production_project_before_cloud_access(self):
        with self.assertRaisesRegex(ValueError, "production-like project"):
            apply_replacement([], "local_check", project="example-prod-project")

    def test_apply_project_must_match_selected_environment(self):
        with patch.dict("os.environ", {"GCLOUD_PROJECT": "approved-dev-project"}):
            with self.assertRaisesRegex(ValueError, "must match GCLOUD_PROJECT"):
                apply_replacement([], "local_check", project="other-dev-project")

    def test_manifest_validation_and_transaction_are_scoped(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = {"datasets": []}
            for dataset_id in ("first_2026", "second_2026"):
                dea = root / f"{dataset_id}.dea.parquet"
                gsea = root / f"{dataset_id}.gsea.parquet"
                _write(
                    dea, DEA_SCHEMA, [dataset_id] if dataset_id == "first_2026" else []
                )
                _write(gsea, GSEA_SCHEMA, [])
                manifest["datasets"].append(
                    {
                        "dataset_id": dataset_id,
                        "dea_parquet": dea.name,
                        "gsea_parquet": gsea.name,
                    }
                )
            manifest_path = root / "manifest.json"
            manifest_path.write_text(json.dumps(manifest))

            files = validate_manifest(manifest_path)
            self.assertEqual(
                [item.dataset_id for item in files], ["first_2026", "second_2026"]
            )
            self.assertEqual(
                [(item.dea_rows, item.gsea_rows) for item in files], [(1, 0), (0, 0)]
            )

            sql = transaction_sql(
                "dev-project",
                "ingest_run_dea",
                "ingest_run_gsea",
            )
            self.assertEqual(sql.count("DELETE FROM"), 2)
            self.assertEqual(sql.count("INSERT INTO"), 2)
            self.assertEqual(sql.count("WHERE dataset_id IN UNNEST(@dataset_ids)"), 4)
            self.assertIn("BEGIN TRANSACTION", sql)
            self.assertIn("COMMIT TRANSACTION", sql)
            self.assertNotIn("TRUNCATE", sql)

    def test_manifest_rejects_rows_for_a_different_id(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dea = root / "expected_2026.dea.parquet"
            gsea = root / "expected_2026.gsea.parquet"
            _write(dea, DEA_SCHEMA, ["other_2026"])
            _write(gsea, GSEA_SCHEMA, [])
            manifest_path = root / "manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "datasets": [
                            {
                                "dataset_id": "expected_2026",
                                "dea_parquet": dea.name,
                                "gsea_parquet": gsea.name,
                            }
                        ]
                    }
                )
            )
            with self.assertRaisesRegex(ValueError, "dataset_id values must all equal"):
                validate_manifest(manifest_path)

    def test_gcs_receipt_is_bound_to_exact_run_dataset_and_object(self):
        run_key = "0123456789abcdef0123456789abcdef"
        dataset_id = "first_2026"

        def file_receipt(kind):
            name = (
                f"perturb-seq-ingest/{run_key}/{dataset_id}/{dataset_id}.{kind}.parquet"
            )
            return {
                "uri": f"gs://dev-tmp/{name}",
                "size_bytes": 123,
                "row_count": 0,
                "sha256": "a" * 64,
                "md5_hash": base64.b64encode(b"m" * 16).decode(),
                "crc32c": base64.b64encode(b"crc!").decode(),
                "generation": "987654321",
            }

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "receipt.json"
            receipt = {
                "run_key": run_key,
                "bucket": "dev-tmp",
                "datasets": [
                    {
                        "dataset_id": dataset_id,
                        "dea": file_receipt("dea"),
                        "gsea": file_receipt("gsea"),
                    }
                ],
            }
            path.write_text(json.dumps(receipt))
            parsed = validate_gcs_manifest(path)
            self.assertEqual(parsed.datasets[0].dea.rows, 0)
            self.assertEqual(parsed.datasets[0].gsea.generation, 987654321)

            receipt["datasets"][0]["dea"]["crc32c"] = None
            path.write_text(json.dumps(receipt))
            self.assertIsNone(validate_gcs_manifest(path).datasets[0].dea.crc32c)

            receipt["datasets"][0]["dea"]["uri"] = "gs://other-bucket/other-object"
            path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, "unexpected GCS URI"):
                validate_gcs_manifest(path)

    def test_gcs_sources_are_rechecked_after_staging_before_replacement(self):
        from types import SimpleNamespace
        from unittest.mock import patch

        events = []

        class SchemaField:
            def __init__(self, name, field_type, mode):
                self.name = name
                self.field_type = field_type
                self.mode = mode

        class BigQuery:
            SourceFormat = SimpleNamespace(PARQUET="PARQUET")
            WriteDisposition = SimpleNamespace(WRITE_TRUNCATE="WRITE_TRUNCATE")
            LoadJobConfig = staticmethod(lambda **kwargs: kwargs)
            QueryJobConfig = staticmethod(lambda **kwargs: kwargs)
            ArrayQueryParameter = staticmethod(lambda *args: args)
            format_options = SimpleNamespace(
                ParquetOptions=lambda: SimpleNamespace(enable_list_inference=False)
            )

        BigQuery.SchemaField = SchemaField

        class Job:
            errors = None

            def result(self):
                return None

        class Client:
            def load_table_from_uri(self, _uris, table_id, **kwargs):
                self_outer.assertTrue(
                    kwargs["job_config"]["parquet_options"].enable_list_inference
                )
                events.append("load:" + table_id.rsplit("_", 1)[-1])
                return Job()

            def get_table(self, table_id):
                kind = "gsea" if table_id.endswith("_gsea") else "dea"
                schema = GSEA_SCHEMA if kind == "gsea" else DEA_SCHEMA
                return SimpleNamespace(
                    schema=_expected_bq_schema(schema, BigQuery), num_rows=1
                )

            def query(self, *_args, **_kwargs):
                events.append("replace")
                return Job()

        def verify_sources():
            events.append("verify")

        self_outer = self
        with (
            patch("load_results_to_bigquery._expire_stage"),
            patch("load_results_to_bigquery._counts", return_value={"demo": 1}),
        ):
            _stage_and_replace(
                Client(),
                BigQuery,
                "dev-project",
                "europe-west2",
                "0123456789abcdef0123456789abcdef",
                ["demo"],
                {"dea": ["gs://bucket/dea"], "gsea": ["gs://bucket/gsea"]},
                {"dea": {"demo": 1}, "gsea": {"demo": 1}},
                verify_sources=verify_sources,
            )

        self.assertEqual(events, ["load:dea", "load:gsea", "verify", "replace"])

    def test_cmek_objects_without_gcs_checksums_are_accepted(self):
        run_key = "0123456789abcdef0123456789abcdef"
        dataset_id = "first_2026"
        objects = {}
        parquet_files = {}
        for kind in ("dea", "gsea"):
            name = (
                f"perturb-seq-ingest/{run_key}/{dataset_id}/{dataset_id}.{kind}.parquet"
            )
            checksum_metadata = {
                "run_key": run_key,
                "dataset_id": dataset_id,
                "kind": kind,
                "sha256": "a" * 64,
                "row_count": "4",
            }
            objects[name] = SimpleNamespace(
                generation="987654321",
                size=123,
                md5_hash=None,
                crc32c=None,
                metadata=checksum_metadata,
            )
            parquet_files[kind] = GcsParquet(
                uri=f"gs://dev-tmp/{name}",
                object_name=name,
                size_bytes=123,
                rows=4,
                sha256="a" * 64,
                md5_hash=base64.b64encode(b"m" * 16).decode(),
                crc32c=None,
                generation=987654321,
            )

        class Bucket:
            def get_blob(self, name):
                return objects[name]

        class StorageClient:
            def bucket(self, _name):
                return Bucket()

        manifest = GcsManifest(
            run_key=run_key,
            bucket="dev-tmp",
            datasets=[
                GcsDatasetFiles(
                    dataset_id=dataset_id,
                    dea=parquet_files["dea"],
                    gsea=parquet_files["gsea"],
                )
            ],
        )
        self.assertEqual(
            _verify_gcs_objects(StorageClient(), manifest),
            {name: 987654321 for name in objects},
        )


if __name__ == "__main__":
    unittest.main()
