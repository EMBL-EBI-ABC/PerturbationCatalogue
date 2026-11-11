#!/usr/bin/env python3
"""Fit MaveDB GMM candidates in whole-study batches."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

PRODUCTION_SCORES = (
    REPOSITORY_ROOT / "data_sources" / "mavedb" / "output" / "processed" / "scores.tsv"
)
PRODUCTION_METADATA = (
    REPOSITORY_ROOT
    / "data_sources"
    / "mavedb"
    / "output"
    / "processed"
    / "metadata.tsv"
)
SMOKE_FIXTURES = (
    REPOSITORY_ROOT / "data_sources" / "mavedb" / "tests" / "fixtures" / "gmm" / "smoke"
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fit flexible or modes GMM candidates for MaveDB score sets."
    )
    parser.add_argument(
        "--profile", choices=("production", "smoke"), default="production"
    )
    parser.add_argument("--study_type", default="maves")
    parser.add_argument("--model", choices=("flexible", "modes"), default="modes")
    parser.add_argument("--index", type=int, default=0)
    parser.add_argument("--variants_per_batch", type=int, default=50_000)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--scores", type=Path)
    parser.add_argument("--meta", type=Path)
    parser.add_argument(
        "--input-format", choices=("catalogue", "legacy"), default="catalogue"
    )
    parser.add_argument(
        "--seed",
        type=int,
        help=(
            "Override the default dataset-derived seed and pass this exact "
            "seed to every Stan fit; also seed row subsampling."
        ),
    )
    parser.add_argument("--list-batches", action="store_true")
    parser.add_argument("--study-id", action="append", dest="study_ids")
    parser.add_argument("--max-studies", type=int)
    parser.add_argument("--max-study-rows", type=int)
    parser.add_argument("--max-fit-rows", type=int)
    parser.add_argument("--chains", type=int)
    parser.add_argument("--parallel-chains", type=int)
    parser.add_argument("--iter-warmup", type=int)
    parser.add_argument("--iter-sampling", type=int)
    parser.add_argument("--memory-budget-mb", type=int)
    parser.add_argument(
        "--summary-json", type=Path, help="Write a portable plan or batch result JSON."
    )
    parser.add_argument(
        "--expected-plan",
        type=Path,
        help="Verify inputs, code and settings against this plan before fitting.",
    )
    parser.add_argument("--fail-on-study-error", action="store_true")
    parser.add_argument(
        "--collect-plan",
        type=Path,
        help="Collect successful batch bundles for this plan, without fitting.",
    )
    parser.add_argument("--batch-results", type=Path, nargs="+")
    return parser


def _option_was_set(argv: Sequence[str], action: argparse.Action) -> bool:
    return any(
        argument == option or argument.startswith(f"{option}=")
        for argument in argv
        for option in action.option_strings
    )


def _explicit_overrides(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    argv: Sequence[str],
) -> dict[str, Any]:
    overrides: dict[str, Any] = {}
    for action in parser._actions:
        if not action.option_strings or action.dest == "help":
            continue
        if not _option_was_set(argv, action):
            continue
        value = getattr(args, action.dest)
        if isinstance(value, Path):
            value = str(value)
        overrides[action.dest] = value
    return overrides


def _resolve_path(path: Path | None, default: Path) -> Path:
    return (path if path is not None else default).expanduser().resolve()


def _build_config(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    argv: Sequence[str],
) -> Any:
    """Build typed settings and BatchConfig without importing the Stan runner."""

    from data_sources.mavedb.analysis._gmm_batch import BatchConfig
    from data_sources.mavedb.analysis.gmm import SamplerSettings

    if args.profile == "smoke":
        defaults = {
            "chains": 2,
            "parallel_chains": 1,
            "iter_warmup": 200,
            "iter_sampling": 200,
            "adapt_delta": 0.9,
            "max_treedepth": 10,
            "refresh": 4000,
            "memory_budget_mb": 2048,
            "seed": None,
            "max_fit_rows": 500,
        }
        max_studies = 1
        max_study_rows = 500
    else:
        defaults = asdict(SamplerSettings())
        max_studies = None
        max_study_rows = None

    setting_overrides = {
        "chains": args.chains,
        "parallel_chains": args.parallel_chains,
        "iter_warmup": args.iter_warmup,
        "iter_sampling": args.iter_sampling,
        "memory_budget_mb": args.memory_budget_mb,
        "seed": args.seed,
        "max_fit_rows": args.max_fit_rows,
    }
    defaults.update(
        {name: value for name, value in setting_overrides.items() if value is not None}
    )
    settings = SamplerSettings(**defaults)

    if args.max_studies is not None:
        max_studies = args.max_studies
    if args.max_study_rows is not None:
        max_study_rows = args.max_study_rows
    if args.profile == "smoke":
        if max_studies is not None and max_studies > 1:
            raise ValueError("smoke profile allows at most one study")
        if max_study_rows is not None and max_study_rows > 500:
            raise ValueError("smoke profile caps selected rows at 500 per study")
        if settings.max_fit_rows > 500:
            raise ValueError("smoke profile caps fitted rows at 500")
        if settings.memory_budget_mb > 2048:
            raise ValueError("smoke profile memory budget cannot exceed 2048 MiB")

    if args.input_format == "legacy" and (args.scores is None or args.meta is None):
        raise ValueError("legacy input requires explicit --scores and --meta paths")

    default_scores = (
        SMOKE_FIXTURES / "scores.tsv" if args.profile == "smoke" else PRODUCTION_SCORES
    )
    default_metadata = (
        SMOKE_FIXTURES / "metadata.tsv"
        if args.profile == "smoke"
        else PRODUCTION_METADATA
    )
    scores_path = _resolve_path(args.scores, default_scores)
    meta_path = _resolve_path(args.meta, default_metadata)
    out_path = _resolve_path(args.out, REPOSITORY_ROOT / "intermediates")

    explicit_overrides = _explicit_overrides(parser, args, argv)
    return BatchConfig(
        scores_path=scores_path,
        meta_path=meta_path,
        input_format=args.input_format,
        profile=args.profile,
        study_type=args.study_type,
        model=args.model,
        index=args.index,
        variants_per_batch=args.variants_per_batch,
        out=out_path,
        study_ids=tuple(args.study_ids or ()),
        max_studies=max_studies,
        max_study_rows=max_study_rows,
        settings=settings,
        list_batches=args.list_batches,
        explicit_overrides=explicit_overrides,
    )


def _settings_dict(settings: Any) -> dict[str, Any]:
    return asdict(settings)


def _effective_config(config: Any) -> dict[str, Any]:
    return {
        "profile": config.profile,
        "study_type": config.study_type,
        "model": config.model,
        "input_format": config.input_format,
        "scores": str(config.scores_path),
        "metadata": str(config.meta_path),
        "out": str(config.out),
        "index": config.index,
        "variants_per_batch": config.variants_per_batch,
        "study_ids": list(config.study_ids),
        "max_studies": config.max_studies,
        "max_study_rows": config.max_study_rows,
        "settings": _settings_dict(config.settings),
        "explicit_overrides": config.explicit_overrides,
    }


def _summary_dict(summary: Any, list_batches: bool) -> dict[str, Any]:
    fields = (
        "profile",
        "model",
        "batch_index",
        "batch_count",
        "eligible_study_count",
        "processed_study_count",
        "successful_study_count",
        "failed_study_count",
        "eligible_rows",
        "score_rows_written",
        "selected_study_ids",
        "assignments",
        "failures",
        "output_paths",
    )
    result = {
        field: getattr(summary, field) for field in fields if hasattr(summary, field)
    }
    result["status"] = "batch_plan" if list_batches else "completed"
    return result


def _json_default(value: Any) -> str:
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "item"):
        return value.item()
    return str(value)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(raw_argv)

    # Reject unsupported cohorts and invalid selection controls before loading
    # input tables or importing the batch implementation.
    if args.study_type != "maves":
        parser.error("only --study_type maves is supported")
    if args.index < 0:
        parser.error("--index must be zero or a positive batch number")
    if args.variants_per_batch <= 0:
        parser.error("--variants_per_batch must be positive")
    if args.max_studies is not None and args.max_studies <= 0:
        parser.error("--max-studies must be positive")
    if args.max_study_rows is not None and args.max_study_rows <= 0:
        parser.error("--max-study-rows must be positive")

    try:
        if args.collect_plan is not None:
            from data_sources.mavedb.analysis._gmm_pipeline import collect_batches

            if args.out is None or not args.batch_results:
                raise ValueError("--collect-plan requires --out and --batch-results")
            collect_batches(args.collect_plan, args.batch_results, args.out)
            return 0
        if (
            args.summary_json is not None
            and not args.list_batches
            and args.expected_plan is None
        ):
            raise ValueError("A batch --summary-json requires --expected-plan")
        config = _build_config(parser, args, raw_argv)
        if args.expected_plan is not None:
            from data_sources.mavedb.analysis._gmm_pipeline import verify_plan

            if args.list_batches:
                raise ValueError(
                    "--expected-plan cannot be combined with --list-batches"
                )
            verify_plan(config, args.expected_plan)
        print(
            json.dumps(
                {"event": "effective_config", **_effective_config(config)},
                sort_keys=True,
                default=_json_default,
            )
        )
        from data_sources.mavedb.analysis._gmm_batch import run_batch

        summary = run_batch(config)
        if args.summary_json is not None:
            from data_sources.mavedb.analysis._gmm_pipeline import (
                batch_record,
                make_plan,
                write_json,
            )

            record = (
                make_plan(config, summary)
                if args.list_batches
                else batch_record(config, summary, args.expected_plan)
            )
            write_json(args.summary_json, record)
        print(
            json.dumps(
                _summary_dict(summary, args.list_batches),
                sort_keys=True,
                default=_json_default,
            )
        )
        if not args.list_batches:
            print(
                f"Processed {summary.processed_study_count} studies: "
                f"{summary.successful_study_count} succeeded, "
                f"{summary.failed_study_count} failed."
            )
        return 3 if args.fail_on_study_error and summary.failed_study_count else 0
    except (ImportError, OSError, RuntimeError, ValueError) as error:
        print(f"fit_gmm_batch: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
