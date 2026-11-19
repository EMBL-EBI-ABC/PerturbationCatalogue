#!/usr/bin/env python3
"""Small Parquet parity checks for the bounded result merger."""

from __future__ import annotations

import hashlib
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "dea-gsea"))

from io_schemas import DEA_SCHEMA, GSEA_SCHEMA, dataframe_to_table, write_parquet
import merge_results
from merge_results import merge_parquet_files, normalize_leading_edge


INGESTED_AT = pd.Timestamp("2026-10-02T17:00:00Z")


def write_source(path: Path, frame: pd.DataFrame, schema: pa.Schema) -> None:
    table = dataframe_to_table(frame, schema)
    with pq.ParquetWriter(path, schema, compression="zstd") as writer:
        writer.write_table(table, row_group_size=2)


def old_merge(
    paths: list[Path],
    schema: pa.Schema,
    sort_columns: list[str],
    output: Path,
    dataset_id: str,
    ingested_at: pd.Timestamp,
    normalize_gsea: bool = False,
) -> pd.DataFrame:
    frames = [pd.read_parquet(path) for path in paths]
    frames = [frame for frame in frames if not frame.empty]
    if frames:
        merged = pd.concat(frames, ignore_index=True)
        merged["dataset_id"] = dataset_id
        merged["max_ingested_at"] = ingested_at
        if normalize_gsea:
            merged["leading_edge"] = merged["leading_edge"].apply(
                normalize_leading_edge
            )
        merged = merged.sort_values(sort_columns).reset_index(drop=True)
    else:
        merged = pd.DataFrame(columns=[field.name for field in schema])
    write_parquet(merged, output, schema)
    return pd.read_parquet(output)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class MergeResultsTest(unittest.TestCase):
    def test_sorted_merge_matches_pandas_reference_and_keeps_empty_schema(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dea_rows = [
                pd.DataFrame(
                    {
                        "dataset_id": ["old", "old", "old"],
                        "perturbed_target_symbol": ["B", "A", "B"],
                        "perturbed_target_ensg": ["E2", "E1", "E2"],
                        "effect_gene_symbol": ["zeta", "zeta", "alpha"],
                        "effect_gene_ensg": ["X2", "X1", "X2"],
                        "padj": [0.2, 0.1, 0.3],
                        "log2foldchange": [1.2, -1.0, 0.5],
                        "score_name": ["old"] * 3,
                        "score_value": [2.0, 1.0, 3.0],
                        "cell_type": [None] * 3,
                        "max_ingested_at": [INGESTED_AT] * 3,
                    }
                ),
                pd.DataFrame(
                    {
                        "dataset_id": ["old", "old", "old", "old"],
                        "perturbed_target_symbol": ["C", "A", None, None],
                        "perturbed_target_ensg": ["E3", "E1", None, None],
                        "effect_gene_symbol": ["beta", "zeta", "zeta", None],
                        "effect_gene_ensg": ["X3", "X1", "X4", None],
                        "padj": pd.Series([0.4, 0.5, np.nan, None], dtype=object),
                        "log2foldchange": [0.4, -0.7, 0.2, 0.0],
                        "score_name": ["old"] * 4,
                        "score_value": [4.0, 5.0, 6.0, 7.0],
                        "cell_type": [None] * 4,
                        "max_ingested_at": [INGESTED_AT] * 4,
                    }
                ),
            ]
            dea_paths = [root / "dea_0.parquet", root / "dea_1.parquet"]
            for path, frame in zip(dea_paths, dea_rows):
                write_source(path, frame, DEA_SCHEMA)

            expected_dea = old_merge(
                dea_paths,
                DEA_SCHEMA,
                ["perturbed_target_symbol", "effect_gene_symbol"],
                root / "expected.dea.parquet",
                "sample",
                INGESTED_AT,
            )
            dea_output = root / "actual.dea.parquet"
            dea_stats = merge_parquet_files(
                [str(path) for path in dea_paths],
                DEA_SCHEMA,
                ["perturbed_target_symbol", "effect_gene_symbol"],
                dea_output,
                "sample",
                INGESTED_AT,
                unique_columns=("perturbed_target_symbol",),
            )
            actual_dea = pd.read_parquet(dea_output)
            pd.testing.assert_frame_equal(actual_dea, expected_dea)
            self.assertIsInstance(actual_dea.index, pd.RangeIndex)
            same_key = actual_dea[
                (actual_dea["perturbed_target_symbol"] == "A")
                & (actual_dea["effect_gene_symbol"] == "zeta")
            ]
            self.assertEqual(same_key["padj"].tolist(), [0.1, 0.5])
            self.assertEqual(dea_stats["rows"], len(expected_dea))
            self.assertEqual(dea_stats["unique_counts"]["perturbed_target_symbol"], 3)
            self.assertEqual(
                pq.ParquetFile(dea_output).schema_arrow.equals(
                    DEA_SCHEMA, check_metadata=True
                ),
                True,
            )
            expected_padj = pq.read_table(root / "expected.dea.parquet")["padj"]
            actual_padj = pq.read_table(dea_output)["padj"]
            self.assertEqual(actual_padj.null_count, expected_padj.null_count)
            self.assertEqual(
                actual_padj.is_null().to_pylist(),
                expected_padj.is_null().to_pylist(),
            )
            self.assertEqual(
                [
                    index
                    for index, value in enumerate(actual_padj.to_pylist())
                    if value is not None and math.isnan(value)
                ],
                [
                    index
                    for index, value in enumerate(expected_padj.to_pylist())
                    if value is not None and math.isnan(value)
                ],
            )
            self.assertEqual(actual_padj.null_count, 0)

            gsea_rows = [
                pd.DataFrame(
                    {
                        "dataset_id": ["old", "old", "old"],
                        "term": ["ZETA", "ALPHA", None],
                        "perturbed_target_symbol": ["B", "A", None],
                        "perturbed_target_ensg": ["E2", "E1", None],
                        "es": [0.1, -0.2, 0.0],
                        "nes": [1.1, -1.2, 0.0],
                        "pval": [0.01, 0.02, 1.0],
                        "sidak": [0.03, 0.04, 1.0],
                        "fdr": [0.05, 0.06, 1.0],
                        "geneset_size": [20, 30, 0],
                        "leading_edge": [["G2"], "G1;G3", None],
                        "cell_type": [None] * 3,
                        "max_ingested_at": [INGESTED_AT] * 3,
                    }
                ),
                pd.DataFrame(
                    {
                        "dataset_id": ["old", "old"],
                        "term": ["ZETA", "BETA"],
                        "perturbed_target_symbol": ["A", "C"],
                        "perturbed_target_ensg": ["E1", "E3"],
                        "es": [0.3, 0.4],
                        "nes": [1.3, 1.4],
                        "pval": [0.03, 0.04],
                        "sidak": [0.05, 0.06],
                        "fdr": [0.07, 0.08],
                        "geneset_size": [25, 40],
                        "leading_edge": [["G4", "G5"], []],
                        "cell_type": [None] * 2,
                        "max_ingested_at": [INGESTED_AT] * 2,
                    }
                ),
            ]
            gsea_paths = [root / "gsea_0.parquet", root / "gsea_1.parquet"]
            for path, frame in zip(gsea_paths, gsea_rows):
                write_source(path, frame, GSEA_SCHEMA)

            expected_gsea = old_merge(
                gsea_paths,
                GSEA_SCHEMA,
                ["perturbed_target_symbol", "term"],
                root / "expected.gsea.parquet",
                "sample",
                INGESTED_AT,
                normalize_gsea=True,
            )
            gsea_output = root / "actual.gsea.parquet"
            gsea_stats = merge_parquet_files(
                [str(path) for path in gsea_paths],
                GSEA_SCHEMA,
                ["perturbed_target_symbol", "term"],
                gsea_output,
                "sample",
                INGESTED_AT,
                unique_columns=("term",),
            )
            actual_gsea = pd.read_parquet(gsea_output)
            pd.testing.assert_frame_equal(actual_gsea, expected_gsea)
            self.assertIsInstance(actual_gsea.index, pd.RangeIndex)
            self.assertEqual(gsea_stats["rows"], len(expected_gsea))
            self.assertEqual(gsea_stats["unique_counts"]["term"], 3)

            empty_output = root / "empty.dea.parquet"
            empty_stats = merge_parquet_files(
                [],
                DEA_SCHEMA,
                ["perturbed_target_symbol", "effect_gene_symbol"],
                empty_output,
                "sample",
                INGESTED_AT,
            )
            expected_empty = old_merge(
                [],
                DEA_SCHEMA,
                ["perturbed_target_symbol", "effect_gene_symbol"],
                root / "expected.empty.dea.parquet",
                "sample",
                INGESTED_AT,
            )
            actual_empty = pd.read_parquet(empty_output)
            pd.testing.assert_frame_equal(actual_empty, expected_empty)
            self.assertIsInstance(actual_empty.index, pd.RangeIndex)
            self.assertEqual(empty_stats["rows"], 0)
            self.assertEqual(
                pq.ParquetFile(empty_output).schema_arrow.equals(
                    DEA_SCHEMA, check_metadata=True
                ),
                True,
            )

            invalid = root / "invalid.parquet"
            pq.write_table(pa.table({"not_a_sort_key": ["x"]}), invalid)
            preserved_hash = sha256(dea_output)
            with self.assertRaises(KeyError):
                merge_parquet_files(
                    [str(invalid)],
                    DEA_SCHEMA,
                    ["perturbed_target_symbol", "effect_gene_symbol"],
                    dea_output,
                    "sample",
                    INGESTED_AT,
                )
            self.assertEqual(sha256(dea_output), preserved_hash)

            manifest_path = root / "manifest.json"
            metrics_path = root / "batch.metrics.json"
            manifest = {"fixture": "metadata"}
            batch_metrics = {"batch_id": "fixture_batch", "n_cells": 3}
            manifest_path.write_text(json.dumps(manifest))
            metrics_path.write_text(json.dumps(batch_metrics))
            output_dir = root / "main_output"
            argv = [
                "merge_results.py",
                "--dataset-id",
                "sample",
                "--outdir",
                str(output_dir),
                "--dea-files",
                *[str(path) for path in dea_paths],
                "--gsea-files",
                *[str(path) for path in gsea_paths],
                "--metrics-files",
                str(metrics_path),
                "--manifest",
                str(manifest_path),
            ]
            with patch.object(sys, "argv", argv):
                merge_results.main()
            summary = json.loads((output_dir / "sample.summary.json").read_text())
            self.assertEqual(summary["n_dea_rows"], len(expected_dea))
            self.assertEqual(summary["n_gsea_rows"], len(expected_gsea))
            self.assertEqual(summary["n_perturbations_in_dea"], 3)
            self.assertEqual(summary["n_terms_in_gsea"], 3)
            self.assertEqual(summary["manifest"], manifest)
            self.assertEqual(summary["batch_metrics"], [batch_metrics])
            cli_dea = pd.read_parquet(output_dir / "sample.dea.parquet")
            cli_gsea = pd.read_parquet(output_dir / "sample.gsea.parquet")
            shared_time = pd.Timestamp(summary["max_ingested_at"])
            self.assertEqual(set(cli_dea["max_ingested_at"]), {shared_time})
            self.assertEqual(set(cli_gsea["max_ingested_at"]), {shared_time})


if __name__ == "__main__":
    unittest.main()
