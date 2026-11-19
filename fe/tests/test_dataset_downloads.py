"""Run with the frontend environment: python -m unittest discover -s tests."""

import sys
import unittest
from pathlib import Path

import dash

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
app = dash.Dash(__name__, use_pages=True, pages_folder="")
from pages.dataset import _dataset_download_controls


class DatasetDownloadsTest(unittest.TestCase):
    def test_h5ad_links_only_appear_for_reprocessed_perturb_seq(self):
        for modality, reprocessed, expected_count in (
            ("perturb-seq", True, 4),
            ("perturb-seq", False, 2),
            ("mave", True, 2),
            ("crispr-screen", True, 2),
        ):
            with self.subTest(modality=modality, reprocessed=reprocessed):
                controls = _dataset_download_controls("demo", modality, reprocessed)
                self.assertEqual(len(controls.children), expected_count)
                if expected_count == 4:
                    self.assertTrue(controls.children[2].href.endswith("?format=h5ad"))
                    self.assertTrue(
                        controls.children[3].href.endswith("?format=filtered.h5ad")
                    )


if __name__ == "__main__":
    unittest.main()
