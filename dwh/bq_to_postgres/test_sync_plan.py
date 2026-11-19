import argparse
import unittest
from datetime import datetime, timezone

from sync_plan import parse_dataset_ids, plan_datasets, validate_force_dataset_ids


class SyncPlanTest(unittest.TestCase):
    def test_forced_ids_are_validated_and_reloaded_without_timestamp_change(self):
        force_ids = parse_dataset_ids("dataset_a")
        timestamp = datetime(2026, 1, 1, tzinfo=timezone.utc)
        bq_info = {
            "dataset_a": (timestamp, 2),
            "dataset_b": (timestamp, 3),
        }

        validate_force_dataset_ids(force_ids, {"dataset_a", "dataset_b"})
        inserts, updates = plan_datasets(
            bq_info,
            {"dataset_a": timestamp, "dataset_b": timestamp},
            force_ids,
            force_updates=True,
        )

        self.assertEqual(inserts, [])
        self.assertEqual(updates, ["dataset_a"])
        inserts, updates = plan_datasets(
            {}, {"dataset_a": timestamp}, force_ids, force_updates=True
        )
        self.assertEqual(inserts, [])
        self.assertEqual(updates, ["dataset_a"])
        with self.assertRaises(ValueError):
            validate_force_dataset_ids(force_ids, {"dataset_b"})
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_dataset_ids("")
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_dataset_ids("dataset_a,,dataset_b")


if __name__ == "__main__":
    unittest.main()
