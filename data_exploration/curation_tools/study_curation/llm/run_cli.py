"""Supported command-line interface for native SQLite curation runs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from curation_tools.study_curation.llm.source_ingestion import build_create_run_request
from curation_tools.study_curation.llm.workflow import CurationWorkflow
from curation_tools.study_curation.paths import (
    CURATION_RUNS_DIR,
    DEFAULT_SCHEMA_PATH,
    DEFAULT_STEP1_PROMPT_PATH,
    DEFAULT_STEP2_PROMPT_PATH,
    DEFAULT_STEP3_PROMPT_PATH,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Native SQLite curation runs")
    parser.add_argument("--runs-root", type=Path, default=CURATION_RUNS_DIR)
    commands = parser.add_subparsers(dest="command", required=True)

    create = commands.add_parser("create", help="Create and snapshot a fresh run")
    create.add_argument("--run-name", required=True)
    create.add_argument("--publication-dir", type=Path, required=True)
    create.add_argument("--mavedb-metadata-dir", type=Path, required=True)
    create.add_argument("--urn-to-dois-file", type=Path, required=True)
    create.add_argument(
        "--step1-prompt",
        type=Path,
        default=DEFAULT_STEP1_PROMPT_PATH,
    )
    create.add_argument(
        "--step2-prompt",
        type=Path,
        default=DEFAULT_STEP2_PROMPT_PATH,
    )
    create.add_argument(
        "--step3-prompt",
        type=Path,
        default=DEFAULT_STEP3_PROMPT_PATH,
    )
    create.add_argument("--schema-path", type=Path, default=DEFAULT_SCHEMA_PATH)
    create.add_argument("--model", default="google/gemini-3.7-flash")
    create.add_argument("--max-workers", type=int, default=8)
    create.add_argument(
        "--urn",
        dest="selected_urns",
        action="append",
        help="Restrict the new run to this MaveDB URN; repeat for multiple URNs.",
    )

    step = commands.add_parser("step", help="Run missing work for one persisted step")
    step.add_argument("--run", required=True)
    step.add_argument(
        "--step", choices=("step1", "step2", "step3", "step4", "step5"), required=True
    )
    step.add_argument("--item", dest="items", action="append")
    step.add_argument("--retry-failed", action="store_true")

    status = commands.add_parser("status", help="Inspect a persisted run")
    status.add_argument("--run", required=True)

    seal = commands.add_parser(
        "seal", help="Make final artifacts immutable and exportable"
    )
    seal.add_argument("--run", required=True)

    export = commands.add_parser(
        "export", help="Write explicit final JSON and CSV deliverables"
    )
    export.add_argument("--run", required=True)
    export.add_argument("--output-dir", type=Path, required=True)
    export.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    workflow = CurationWorkflow(args.runs_root)
    if args.command == "create":
        summary = workflow.create_run(
            build_create_run_request(
                run_name=args.run_name,
                publication_dir=args.publication_dir,
                mavedb_metadata_dir=args.mavedb_metadata_dir,
                urn_to_dois_file=args.urn_to_dois_file,
                step1_prompt=args.step1_prompt,
                step2_prompt=args.step2_prompt,
                step3_prompt=args.step3_prompt,
                schema_path=args.schema_path,
                model=args.model,
                max_workers=args.max_workers,
                selected_urns=args.selected_urns,
            )
        )
        print(json.dumps({"run": summary.run_name, "state": summary.state}))
    elif args.command == "step":
        summary = workflow.run_step(args.run, args.step, args.items, args.retry_failed)
        print(
            json.dumps(
                {
                    "step": summary.step,
                    "state": summary.state,
                    "execution_id": summary.execution_id,
                }
            )
        )
    elif args.command == "status":
        status = workflow.get_status(args.run)
        print(
            json.dumps(
                {
                    "run": status.summary.run_name,
                    "state": status.summary.state,
                    "revision": status.summary.revision,
                    "items": len(status.items),
                    "executions": status.executions,
                },
                default=str,
            )
        )
    elif args.command == "seal":
        summary = workflow.seal_run(args.run)
        print(json.dumps({"run": summary.run_name, "state": summary.state}))
    else:
        result = workflow.export_final(args.run, args.output_dir, args.overwrite)
        print(
            json.dumps(
                {"export_id": result.export_id, "output_dir": str(result.output_dir)}
            )
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
