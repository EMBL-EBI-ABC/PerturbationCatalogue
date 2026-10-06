import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import pyarrow.parquet as pq

from release_worker import generate


class Blob:
    def __init__(self, root, name):
        self.path = root / name
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def open(self, mode, **_):
        return self.path.open(mode)

    def delete(self):
        self.path.unlink(missing_ok=True)


class Bucket:
    def __init__(self, root):
        self.root = root

    def blob(self, name):
        return Blob(self.root, name)


class Result:
    def __init__(self, rows):
        self.rows = rows

    def to_arrow_iterable(self, **_):
        raise ImportError

    def __iter__(self):
        return iter(self.rows)


class Job:
    def __init__(self, rows):
        self.rows = rows

    def result(self, **_):
        return Result(self.rows)


class Client:
    def query(self, query, **_):
        if "metadata_table" in query:
            self.rows = [{"dataset_id": "demo", "sample_id": "s1"}]
        else:
            self.rows = [
                {
                    "perturbed_target_ensg": "ENSG1",
                    "perturbed_target_name": "G1",
                    "score_name": "score",
                    "score_value": 1.5,
                    "significant": True,
                    "significance_criteria": "p",
                }
            ]
        return Job(self.rows)

    def list_rows(self, _):
        return Result(self.rows)

    def delete_table(self, *_args, **_kwargs):
        pass


class TestReleaseWorker(unittest.TestCase):
    def test_generate(self):
        prefix = "release/crispr/demo"
        item = {
            "modality": "crispr",
            "dataset_id": "demo",
            "data_table": "project.dataset.data_table",
            "metadata_table": "project.dataset.metadata_table",
        }
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            generate(item, Bucket(tmp_path), prefix, Client(), "test", None)

            with gzip.open(tmp_path / f"{prefix}.csv.gz", "rt", newline="") as handle:
                self.assertEqual(next(csv.reader(handle))[0], "Perturbed Target ENSG")
            self.assertEqual(pq.read_table(tmp_path / f"{prefix}.parquet").num_rows, 1)
            self.assertEqual(
                json.loads((tmp_path / f"{prefix}.metadata.json").read_text())[
                    "dataset_id"
                ],
                "demo",
            )


if __name__ == "__main__":
    unittest.main()
