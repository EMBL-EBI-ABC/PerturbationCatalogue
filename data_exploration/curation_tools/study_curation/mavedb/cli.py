"""Command-line entry points for MaveDB metadata post-processing."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from curation_tools.study_curation.mavedb.workflow import (
    DEFAULT_MAVEDB_CSV_DIR,
    publish_to_bigquery,
    run_pipeline,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Post-process curated MaveDB metadata")
    commands = parser.add_subparsers(dest="command", required=True)

    run = commands.add_parser("run", help="Create local curated parquet outputs")
    run.add_argument("--llm-metadata", type=Path, required=True)
    run.add_argument("--manual-metadata", type=Path)
    run.add_argument("--mavedb-csv-dir", type=Path, default=DEFAULT_MAVEDB_CSV_DIR)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument(
        "--save-joint-artifact",
        action="store_true",
        help="Also save score data with validated metadata as a joint parquet",
    )

    publish = commands.add_parser(
        "publish", help="Merge already-reviewed final parquet outputs into BigQuery"
    )
    publish.add_argument("--metadata-parquet", type=Path, required=True)
    publish.add_argument("--data-parquet", type=Path, required=True)
    publish.add_argument("--project-id")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        manifest = run_pipeline(
            llm_metadata_path=args.llm_metadata,
            manual_metadata_path=args.manual_metadata,
            mavedb_csv_dir=args.mavedb_csv_dir,
            output_dir=args.output_dir,
            save_joint_artifact=args.save_joint_artifact,
        )
        print(json.dumps(manifest, indent=2))
    else:
        publish_to_bigquery(
            metadata_path=args.metadata_parquet,
            data_path=args.data_parquet,
            project_id=args.project_id,
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
