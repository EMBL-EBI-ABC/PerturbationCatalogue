#!/usr/bin/env python3
"""Check SRA locator pinning and checked-range delegation."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from stream_count import download_sra


class SraLocatorPinTest(unittest.TestCase):
    @staticmethod
    def locator_execute(accession, size, digest, url, commands):
        def execute(command, _log_path):
            commands.append(command)
            metadata_path = Path(command[command.index("--output") + 1])
            metadata_path.write_text(
                json.dumps(
                    {
                        "result": [
                            {
                                "bundle": accession,
                                "status": 200,
                                "files": [
                                    {
                                        "type": "sra",
                                        "accession": accession,
                                        "size": size,
                                        "md5": digest,
                                        "locations": [{"link": url}],
                                    }
                                ],
                            }
                        ]
                    }
                )
            )

        return execute

    def test_mismatch_stops_before_archive_transfer(self):
        accession = "SRR123456"
        archive = b"small test archive"
        archive_md5 = hashlib.md5(archive).hexdigest()

        for pins in (
            {"expected_bytes": len(archive) + 1, "expected_md5": archive_md5},
            {"expected_bytes": len(archive), "expected_md5": "0" * 32},
        ):
            with self.subTest(pins=pins), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                downloads, logs = root / "downloads", root / "logs"
                downloads.mkdir()
                logs.mkdir()
                commands = []
                execute = self.locator_execute(
                    accession,
                    len(archive),
                    archive_md5,
                    "https://download.example.invalid/archive.sra",
                    commands,
                )

                with patch("stream_count.download_ranges") as ranged_download:
                    with self.assertRaisesRegex(RuntimeError, "Pinned SRA metadata"):
                        download_sra(accession, downloads, logs, execute, **pins)
                    ranged_download.assert_not_called()

                self.assertEqual(
                    len(commands), 1, "archive ranges started after a pin mismatch"
                )
                target = downloads / accession
                self.assertFalse((target / "locator.json").exists())
                self.assertFalse(list(target.glob("part-*")))
                self.assertFalse(list(target.glob("*.sra")))

    def test_large_pinned_archive_uses_checked_range_helper(self):
        accession = "SRR987654"
        size = 111_880_000_000
        digest = "a" * 32
        url = "https://download.example.invalid/large.sra"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            downloads, logs = root / "downloads", root / "logs"
            downloads.mkdir()
            logs.mkdir()
            commands = []
            execute = self.locator_execute(accession, size, digest, url, commands)
            returned = {
                "bytes": size,
                "md5": digest,
                "ranges": 834,
                "workers": 32,
                "seconds": 1.0,
            }

            with patch(
                "stream_count.download_ranges", return_value=returned
            ) as ranged_download:
                metrics = download_sra(
                    accession,
                    downloads,
                    logs,
                    execute,
                    expected_bytes=size,
                    expected_md5=digest,
                )

            ranged_download.assert_called_once_with(
                url,
                size,
                digest,
                downloads / accession / f"{accession}.sra",
                workers=32,
            )
            self.assertEqual(
                metrics,
                {
                    "archive_bytes": size,
                    "archive_md5": digest,
                    "download_connections": 32,
                },
            )
            self.assertEqual(len(commands), 1, "only the tiny locator uses curl")
            self.assertFalse((downloads / accession / "locator.json").exists())


if __name__ == "__main__":
    unittest.main()
