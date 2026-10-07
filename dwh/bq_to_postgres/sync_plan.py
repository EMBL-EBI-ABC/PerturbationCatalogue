import argparse
import re
from datetime import timezone


_DATASET_ID = re.compile(r"[A-Za-z0-9_-]+\Z")


def parse_dataset_ids(value):
    dataset_ids = value.split(",") if value else []
    if not dataset_ids or any(not _DATASET_ID.fullmatch(item) for item in dataset_ids):
        raise argparse.ArgumentTypeError(
            "must be a non-empty comma-separated list of valid dataset IDs"
        )
    return set(dataset_ids)


def validate_force_dataset_ids(dataset_ids, metadata_dataset_ids):
    missing = dataset_ids - set(metadata_dataset_ids)
    if missing:
        raise ValueError(
            "Force dataset IDs must exist in perturb_seq.metadata; missing: "
            + ", ".join(sorted(missing))
        )


def plan_datasets(bq_info, pg_info, force_dataset_ids, force_updates=False):
    to_insert, to_update = [], []
    for dataset_id, (bq_ts, _row_count) in bq_info.items():
        if dataset_id not in pg_info:
            to_insert.append(dataset_id)
            continue

        pg_ts = pg_info[dataset_id]
        if pg_ts and pg_ts.tzinfo is None:
            pg_ts = pg_ts.replace(tzinfo=timezone.utc)
        if (force_updates and dataset_id in force_dataset_ids) or bq_ts > pg_ts:
            to_update.append(dataset_id)
    if force_updates:
        for dataset_id in sorted(force_dataset_ids - set(bq_info)):
            (to_update if dataset_id in pg_info else to_insert).append(dataset_id)
    return to_insert, to_update
