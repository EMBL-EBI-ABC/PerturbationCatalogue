#!/usr/bin/env python3
"""Focused local checks for scoped DEA/GSEA replacement inputs and SQL."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "dea-gsea"))

from io_schemas import DEA_SCHEMA, GSEA_SCHEMA
from load_results_to_bigquery import (
    apply_replacement,
    transaction_sql,
    validate_manifest,
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


if __name__ == "__main__":
    unittest.main()
