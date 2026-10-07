#!/usr/bin/env python3
"""Protocol check for cluster-to-GCS resumable Parquet upload."""

from __future__ import annotations

import base64
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.parse import quote

import pyarrow as pa
import pyarrow.parquet as pq

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "dea-gsea"))

import cluster_upload
from io_schemas import DEA_SCHEMA


class ResumableUploadTest(unittest.TestCase):
    def test_session_url_must_name_the_expected_bucket_and_object(self):
        dataset_id = "example_2026"
        object_name = (
            "perturb-seq-ingest/0123456789abcdef0123456789abcdef/"
            f"{dataset_id}/{dataset_id}.dea.parquet"
        )
        source_entry = {
            "object_name": object_name,
            "session_uri": (
                "https://storage.googleapis.com/upload/storage/v1/b/dev-tmp/o"
                "?uploadType=resumable&name=other.parquet&upload_id=fake"
            ),
        }
        with patch.object(cluster_upload, "_http_put") as http_put:
            with self.assertRaisesRegex(cluster_upload.TransferError, "object name"):
                cluster_upload._upload_one(
                    "0123456789abcdef0123456789abcdef",
                    "dev-tmp",
                    dataset_id,
                    "dea",
                    source_entry,
                    [],
                )
        http_put.assert_not_called()

    def test_chunk_upload_recovers_after_lost_response(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset_id = "example_2026"
            source = root / f"{dataset_id}.dea.parquet"
            rows = [
                {"dataset_id": dataset_id, "perturbed_target_symbol": "x" * 310_000},
                {"dataset_id": dataset_id, "perturbed_target_symbol": "y" * 310_000},
            ]
            pq.write_table(
                pa.Table.from_pylist(rows, schema=DEA_SCHEMA),
                source,
                compression="NONE",
            )
            content = source.read_bytes()
            self.assertGreater(len(content), 2 * 256 * 1024)

            accepted = bytearray()
            lost_once = False
            object_name = (
                "perturb-seq-ingest/0123456789abcdef0123456789abcdef/"
                f"{dataset_id}/{dataset_id}.dea.parquet"
            )
            session_uri = (
                "https://storage.googleapis.com/upload/storage/v1/b/dev-tmp/o"
                f"?uploadType=resumable&name={quote(object_name, safe='')}&upload_id=fake"
            )

            def fake_put(_url, body, headers):
                nonlocal lost_once
                content_range = headers["Content-Range"]
                if content_range.startswith("bytes */"):
                    response_headers = (
                        {"Range": f"bytes=0-{len(accepted) - 1}"} if accepted else {}
                    )
                    return 308, response_headers, b""
                start, end = map(
                    int,
                    content_range.removeprefix("bytes ").split("/")[0].split("-"),
                )
                self.assertEqual(start, len(accepted))
                self.assertEqual(end - start + 1, len(body))
                accepted.extend(body)
                if start > 0 and not lost_once:
                    lost_once = True
                    raise cluster_upload.RetryableTransferError("lost response")
                if len(accepted) == len(content):
                    self.assertEqual(
                        headers["X-Goog-Hash"],
                        "md5="
                        + base64.b64encode(hashlib.md5(content).digest()).decode(),
                    )
                    final = {
                        "bucket": "dev-tmp",
                        "name": object_name,
                        "size": str(len(content)),
                        "generation": "987654321",
                    }
                    return 201, {}, json.dumps(final).encode()
                return 308, {"Range": f"bytes=0-{len(accepted) - 1}"}, b""

            source_entry = {
                "source_path": str(source),
                "size_bytes": len(content),
                "object_name": object_name,
                "session_uri": session_uri,
            }
            schema = [
                [field.name, str(field.type), field.nullable] for field in DEA_SCHEMA
            ]
            with (
                patch.object(cluster_upload, "RESULTS_ROOT", str(root)),
                patch.object(cluster_upload, "CHUNK_SIZE", 256 * 1024),
                patch.object(cluster_upload, "_http_put", fake_put),
                patch.object(cluster_upload.time, "sleep"),
            ):
                receipt = cluster_upload._upload_one(
                    "0123456789abcdef0123456789abcdef",
                    "dev-tmp",
                    dataset_id,
                    "dea",
                    source_entry,
                    schema,
                )

            self.assertTrue(lost_once)
            self.assertEqual(bytes(accepted), content)
            self.assertEqual(receipt["size_bytes"], len(content))
            self.assertEqual(receipt["row_count"], 2)
            self.assertEqual(
                receipt["md5_hash"],
                base64.b64encode(hashlib.md5(content).digest()).decode(),
            )
            self.assertIsNone(receipt["crc32c"])
            self.assertEqual(receipt["generation"], "987654321")


if __name__ == "__main__":
    unittest.main()
