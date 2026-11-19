"""Run with the frontend environment: python -m unittest discover -s tests."""

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import dash

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
app = dash.Dash(__name__, use_pages=True, pages_folder="")
from pages.dataset import render_dataset


class DatasetDownloadsTest(unittest.TestCase):
    def test_h5ad_links_share_metadata_row_only_for_reprocessed_perturb_seq(self):
        for modality, reprocessed, expected_count in (
            ("Perturb-seq", True, 3),
            ("Perturb-seq", False, 1),
            ("MAVE", True, 1),
            ("CRISPR screen", True, 1),
        ):
            with self.subTest(modality=modality, reprocessed=reprocessed):
                metadata = {
                    "data_modalities": [modality],
                    "perturb_seq_reprocessed": reprocessed,
                }
                with patch(
                    "pages.dataset.fetch_dataset", return_value=(metadata, None)
                ):
                    content, _, _ = render_dataset({"dataset_id": "demo"})
                top_row = content.children[1]
                self.assertEqual(len(top_row.children), expected_count)
                self.assertTrue(top_row.children[0].href.endswith("?format=metadata"))
                if expected_count == 3:
                    self.assertTrue(top_row.children[1].href.endswith("?format=h5ad"))
                    self.assertTrue(
                        top_row.children[2].href.endswith("?format=filtered.h5ad")
                    )
                data_section = content.children[3]
                self.assertEqual(
                    data_section.children[1].children,
                    (
                        "Differential Expression (DEA)"
                        if modality == "Perturb-seq"
                        else "Data"
                    ),
                )
                data_links = data_section.children[3].children
                self.assertEqual(len(data_links), 2)
                label = "DEA" if modality == "Perturb-seq" else "full data"
                self.assertEqual(
                    data_links[0].children.children[-1], f"Download {label} (Parquet)"
                )


if __name__ == "__main__":
    unittest.main()
