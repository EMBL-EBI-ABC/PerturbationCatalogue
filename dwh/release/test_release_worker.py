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
    def __init__(self):
        self.destinations = []

    def query(self, query, **_):
        self.destinations.append(_.get("job_config").destination)
        if "metadata_table" in query:
            self.rows = [{"dataset_id": "demo", "sample_id": "s1"}]
        elif "gsea_table" in query:
            self.rows = [
                {
                    "perturbed_target_ensg": "ENSG1",
                    "perturbed_target_symbol": "G1",
                    "term": "Pathway A",
                    "es": 0.5,
                    "nes": 1.2,
                    "pval": 0.01,
                    "sidak": 0.03,
                    "fdr": 0.02,
                    "geneset_size": 12,
                    "leading_edge": ["A", "B"],
                    "cell_type": "T cell",
                }
            ]
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

    def test_generate_writes_gsea_parquet_and_csv(self):
        prefix = "release/perturb-seq/demo"
        item = {
            "modality": "perturb-seq",
            "dataset_id": "demo",
            "data_table": "project.dataset.data_table",
            "gsea_table": "project.dataset.gsea_table",
            "metadata_table": "project.dataset.metadata_table",
        }
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            client = Client()
            generate(item, Bucket(tmp_path), prefix, client, "test", None)

            gsea = pq.read_table(tmp_path / f"{prefix}.gsea.parquet")
            self.assertTrue(
                str(gsea.schema.field("leading_edge").type).startswith("list")
            )
            self.assertEqual(gsea.column("leading_edge").to_pylist(), [["A", "B"]])
            with gzip.open(
                tmp_path / f"{prefix}.gsea.csv.gz", "rt", newline=""
            ) as handle:
                header, row = list(csv.reader(handle))
            self.assertEqual(header[9], "Leading Edge")
            self.assertEqual(row[9], "A;B")
            data_destination, gsea_destination = client.destinations[:2]
            self.assertNotEqual(data_destination.table_id, gsea_destination.table_id)
            self.assertTrue(data_destination.table_id.endswith("_data"))
            self.assertTrue(gsea_destination.table_id.endswith("_gsea"))


if __name__ == "__main__":
    unittest.main()
