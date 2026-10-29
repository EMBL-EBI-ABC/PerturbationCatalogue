"""Batch configuration, staged output writing, and study orchestration."""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from .gmm import (
    BatchPlan,
    SamplerSettings,
    StudyFit,
    family_schemas,
    fit_study,
    plan_batches,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
PRODUCTION_SCORES = REPOSITORY_ROOT / "data_sources/mavedb/output/processed/scores.tsv"
PRODUCTION_METADATA = (
    REPOSITORY_ROOT / "data_sources/mavedb/output/processed/metadata.tsv"
)
SMOKE_FIXTURES = REPOSITORY_ROOT / "data_sources/mavedb/tests/fixtures/gmm/smoke"
SMOKE_SCORES = SMOKE_FIXTURES / "scores.tsv"
SMOKE_METADATA = SMOKE_FIXTURES / "metadata.tsv"

OUTPUT_FILENAMES = {
    "scores": "score{index}.tsv",
    "densities": "density{index}.tsv",
    "params": "fitted_params{index}.tsv",
    "posterior": "posterior{index}.tsv",
    "fit_summary": "fit_summary{index}.tsv",
    "loo": "loo{index}.tsv",
}

QC_STUDY_COLUMNS = (
    "study_id",
    "input_rows",
    "eligible_rows",
    "selected_rows",
    "fit_rows",
    "excluded_rows",
    "invalid_rows",
    "subset_excluded_rows",
    "success",
    "top_model",
    "failure_reason",
    "candidate_warnings",
    "pareto_k_summary",
    "estimated_peak_memory_mb",
    "selected_text_bytes",
    "profile",
    "max_study_rows",
    "smoke_subset_applied",
)
QC_SMOKE_ID_COLUMNS = ("study_id", "dataset_id", "source_row", "fit_sample")


@dataclass(frozen=True)
class BatchConfig:
    """Validated paths, selection, and sampler settings for a GMM batch."""

    scores_path: Path | str | None = None
    meta_path: Path | str | None = None
    input_format: str = "catalogue"
    profile: str = "production"
    study_type: str = "maves"
    model: str = "modes"
    index: int = 0
    variants_per_batch: int = 50_000
    out: Path | str | None = None
    study_ids: tuple[str, ...] = ()
    max_studies: int | None = None
    max_study_rows: int | None = None
    settings: SamplerSettings | None = None
    list_batches: bool = False
    explicit_overrides: Mapping[str, Any] = field(default_factory=dict)
    score_chunk_size: int = 10_000

    def __post_init__(self) -> None:
        if self.input_format not in {"catalogue", "legacy"}:
            raise ValueError("input_format must be 'catalogue' or 'legacy'")
        if self.profile not in {"production", "smoke"}:
            raise ValueError("profile must be 'production' or 'smoke'")
        if self.study_type != "maves":
            raise ValueError(
                "study_type must be 'maves'; other study types are unsupported"
            )
        if self.model not in {"flexible", "modes"}:
            raise ValueError("model must be 'flexible' or 'modes'")
        if (
            not isinstance(self.index, int)
            or isinstance(self.index, bool)
            or self.index < 0
        ):
            raise ValueError("index must be a nonnegative integer")
        if (
            not isinstance(self.variants_per_batch, int)
            or isinstance(self.variants_per_batch, bool)
            or self.variants_per_batch <= 0
        ):
            raise ValueError("variants_per_batch must be a positive integer")
        if (
            not isinstance(self.score_chunk_size, int)
            or isinstance(self.score_chunk_size, bool)
            or self.score_chunk_size <= 0
        ):
            raise ValueError("score_chunk_size must be a positive integer")

        study_ids = (
            (self.study_ids,)
            if isinstance(self.study_ids, str)
            else tuple(self.study_ids)
        )
        normalized_study_ids = tuple(str(study_id) for study_id in study_ids)
        if any(not study_id for study_id in normalized_study_ids):
            raise ValueError("study IDs must be nonempty")
        if len(set(normalized_study_ids)) != len(normalized_study_ids):
            raise ValueError("study_ids cannot contain duplicates")
        object.__setattr__(self, "study_ids", normalized_study_ids)

        if self.input_format == "legacy" and (
            self.scores_path is None or self.meta_path is None
        ):
            raise ValueError("legacy input requires explicit scores_path and meta_path")
        if self.profile == "smoke":
            default_scores = SMOKE_SCORES
            default_metadata = SMOKE_METADATA
        else:
            default_scores = PRODUCTION_SCORES
            default_metadata = PRODUCTION_METADATA
        scores_path = _resolve_path(
            default_scores if self.scores_path is None else self.scores_path
        )
        meta_path = _resolve_path(
            default_metadata if self.meta_path is None else self.meta_path
        )
        out_path = _resolve_path(
            REPOSITORY_ROOT / "intermediates" if self.out is None else self.out
        )
        object.__setattr__(self, "scores_path", scores_path)
        object.__setattr__(self, "meta_path", meta_path)
        object.__setattr__(self, "out", out_path)

        max_studies = self.max_studies
        max_study_rows = self.max_study_rows
        settings = self.settings
        if max_studies is not None and (
            not isinstance(max_studies, int)
            or isinstance(max_studies, bool)
            or max_studies <= 0
        ):
            raise ValueError("max_studies must be a positive integer when supplied")
        if max_study_rows is not None and (
            not isinstance(max_study_rows, int)
            or isinstance(max_study_rows, bool)
            or max_study_rows <= 0
        ):
            raise ValueError("max_study_rows must be a positive integer when supplied")
        if self.profile == "smoke":
            if max_studies is None:
                max_studies = 1
            if max_study_rows is None:
                max_study_rows = 500
            if settings is None:
                settings = _smoke_sampler_settings()
            if not isinstance(settings, SamplerSettings):
                raise TypeError("settings must be a SamplerSettings instance")
            if max_studies > 1:
                raise ValueError("smoke profile supports at most one study")
            if max_study_rows > 500:
                raise ValueError("smoke profile supports at most 500 study rows")
            if settings.memory_budget_mb > 2048:
                raise ValueError("smoke profile memory_budget_mb cannot exceed 2048")
            if settings.max_fit_rows > 500:
                raise ValueError("smoke profile max_fit_rows cannot exceed 500")
            if (
                settings.chains > 2
                or settings.parallel_chains > 1
                or settings.iter_warmup > 200
                or settings.iter_sampling > 200
            ):
                raise ValueError("smoke sampler overrides cannot exceed smoke defaults")
        elif settings is None:
            settings = SamplerSettings()
        if max_studies is not None and len(normalized_study_ids) > max_studies:
            raise ValueError("study_ids cannot exceed max_studies")
        if not isinstance(settings, SamplerSettings):
            raise TypeError("settings must be a SamplerSettings instance")
        if settings.max_fit_rows > 256_000:
            raise ValueError("max_fit_rows cannot exceed the production cap of 256000")
        object.__setattr__(self, "max_studies", max_studies)
        object.__setattr__(self, "max_study_rows", max_study_rows)
        object.__setattr__(self, "settings", settings)
        object.__setattr__(self, "explicit_overrides", dict(self.explicit_overrides))

    @property
    def output_dir(self) -> Path:
        if self.profile == "smoke":
            return Path(self.out) / "smoke" / "maves" / self.model
        return Path(self.out) / "maves" / self.model

    def validate_input_paths(self) -> None:
        for label, path in (("scores", self.scores_path), ("metadata", self.meta_path)):
            if not Path(path).is_file():
                raise FileNotFoundError(f"{label} input does not exist: {path}")


@dataclass(frozen=True)
class StudyFailure:
    study_id: str
    reason: str


@dataclass(frozen=True)
class BatchSummary:
    profile: str
    model: str
    batch_index: int
    batch_count: int
    eligible_study_count: int
    processed_study_count: int
    successful_study_count: int
    failed_study_count: int
    eligible_rows: int
    score_rows_written: int
    selected_study_ids: tuple[str, ...]
    assignments: tuple[tuple[int, str, int], ...]
    failures: tuple[StudyFailure, ...]
    output_paths: Mapping[str, Path]
    effective_settings: Mapping[str, Any]
    list_batches: bool = False
    plan: BatchPlan | None = None


def estimate_study_memory_mb(
    row_count: int, settings: SamplerSettings, selected_text_bytes: int = 0
) -> int:
    """Conservative peak-RAM estimate for one fit and streamed postprocessing.

    The CmdStan adapter stores generated likelihoods in a disk-backed memmap
    and computes PSIS in bounded observation blocks. Budget two block-sized
    allowances for mapped pages and the copied workspace, then add the
    per-observation PSIS arrays, selected study frames, and a conservative
    Python/CmdStan baseline.
    """

    if not isinstance(row_count, int) or isinstance(row_count, bool) or row_count < 0:
        raise ValueError("row_count must be a nonnegative integer")
    if (
        not isinstance(selected_text_bytes, int)
        or isinstance(selected_text_bytes, bool)
        or selected_text_bytes < 0
    ):
        raise ValueError("selected_text_bytes must be a nonnegative integer")
    python_and_stan_mb = 384 + 128 * settings.parallel_chains
    study_frames_mb = math.ceil(row_count * 2048 / (1024**2))
    # _observation_block_size caps each mapped/copy allowance at 64 MiB and
    # otherwise limits it to one eighth of the configured memory budget.
    likelihood_block_mb = min(64, settings.memory_budget_mb / 8)
    likelihood_mb = math.ceil(2 * likelihood_block_mb)
    selected_text_mb = math.ceil(selected_text_bytes * 4 / (1024**2))
    loo_workspace_mb = math.ceil(
        settings.chains * settings.iter_sampling * 8 * 16 / (1024**2)
    )
    scalar_draw_mb = math.ceil(
        settings.chains * settings.iter_sampling * 8 * 20 / (1024**2)
    )
    return (
        python_and_stan_mb
        + study_frames_mb
        + selected_text_mb
        + likelihood_mb
        + loo_workspace_mb
        + scalar_draw_mb
    )


def run_batch(config: BatchConfig, runner: Any | None = None) -> BatchSummary:
    """Read selected score sets, fit them, and publish six staged TSV files."""

    if not isinstance(config, BatchConfig):
        raise TypeError("config must be a BatchConfig instance")
    config.validate_input_paths()
    prepared_context = _prepare_inputs(config)
    with prepared_context as prepared:
        plan = plan_batches(prepared.study_counts, config.variants_per_batch)
        selected_study_ids = _select_studies(plan, config.index)
        selected_rows = sum(
            plan.study_counts[study_id] for study_id in selected_study_ids
        )
        if config.list_batches:
            return BatchSummary(
                profile=config.profile,
                model=config.model,
                batch_index=config.index,
                batch_count=plan.num_batches,
                eligible_study_count=len(plan.study_counts),
                processed_study_count=0,
                successful_study_count=0,
                failed_study_count=0,
                eligible_rows=selected_rows,
                score_rows_written=0,
                selected_study_ids=selected_study_ids,
                assignments=plan.assignments,
                failures=(),
                output_paths={},
                effective_settings=_effective_settings(config),
                list_batches=True,
                plan=plan,
            )
        if not selected_study_ids:
            raise ValueError("no eligible studies were selected")

        memory_estimates = {
            study_id: estimate_study_memory_mb(
                min(
                    plan.study_counts[study_id],
                    config.max_study_rows or plan.study_counts[study_id],
                ),
                config.settings,
                selected_text_bytes=int(
                    _first_value(
                        _source_study_qc(prepared, study_id),
                        ("selected_text_bytes",),
                        0,
                    )
                ),
            )
            for study_id in selected_study_ids
        }
        over_budget = {
            study_id: estimated
            for study_id, estimated in memory_estimates.items()
            if estimated > config.settings.memory_budget_mb
        }
        if over_budget:
            study_id, estimated = max(over_budget.items(), key=lambda item: item[1])
            raise MemoryError(
                f"estimated peak memory for {study_id} is {estimated} MiB, exceeding "
                f"the {config.settings.memory_budget_mb} MiB budget"
            )

        identity_columns = (
            ("dataset_id", "source_row") if config.input_format == "catalogue" else ()
        )
        schemas = family_schemas(config.model, identity_columns)
        return _fit_and_publish(
            config,
            prepared,
            plan,
            selected_study_ids,
            memory_estimates,
            schemas,
            runner,
        )


def _prepare_inputs(config: BatchConfig) -> Any:
    from ._gmm_inputs import prepare_inputs

    def stage_selector(counts: Mapping[str, int]) -> tuple[str, ...]:
        plan = plan_batches(counts, config.variants_per_batch)
        return _select_studies(plan, config.index)

    return prepare_inputs(
        config.scores_path,
        config.meta_path,
        input_format=config.input_format,
        study_ids=config.study_ids,
        max_studies=config.max_studies,
        max_study_rows=config.max_study_rows,
        seed=config.settings.seed,
        chunk_size=5000,
        stage_selector=None if config.list_batches else stage_selector,
        counts_only=config.list_batches,
    )


def _fit_and_publish(
    config: BatchConfig,
    prepared: Any,
    plan: BatchPlan,
    selected_study_ids: tuple[str, ...],
    memory_estimates: Mapping[str, int],
    schemas: Mapping[str, tuple[str, ...]],
    runner: Any | None,
) -> BatchSummary:
    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    failures: list[StudyFailure] = []
    processed_count = 0
    successful_count = 0
    failed_count = 0
    score_rows_written = 0
    selected_seen: set[str] = set()
    output_paths = {
        table: output_dir / filename.format(index=config.index)
        for table, filename in OUTPUT_FILENAMES.items()
    }

    with tempfile.TemporaryDirectory(prefix=".gmm-stage-", dir=output_dir) as temp_dir:
        stage_dir = Path(temp_dir)
        stage_paths = {
            table: stage_dir / path.name for table, path in output_paths.items()
        }
        for table, stage_path in stage_paths.items():
            _append_frame(
                pd.DataFrame(columns=schemas[table]), stage_path, schemas[table]
            )

        qc_stage_dir = stage_dir / "qc" / f"batch{config.index}"
        qc_stage_dir.mkdir(parents=True, exist_ok=True)
        qc_study_path = qc_stage_dir / "study_diagnostics.tsv"
        _append_frame(
            pd.DataFrame(columns=QC_STUDY_COLUMNS), qc_study_path, QC_STUDY_COLUMNS
        )
        qc_smoke_path = None
        if config.profile == "smoke" and config.input_format == "catalogue":
            qc_smoke_path = qc_stage_dir / "smoke_selected_rows.tsv"
            _append_frame(
                pd.DataFrame(columns=QC_SMOKE_ID_COLUMNS),
                qc_smoke_path,
                QC_SMOKE_ID_COLUMNS,
            )

        for study_id, rows in prepared.iter_studies(selected_study_ids):
            study_id = str(study_id)
            if study_id not in selected_study_ids:
                raise ValueError(f"input adapter returned unselected study {study_id}")
            if study_id in selected_seen:
                raise ValueError(
                    f"input adapter returned study {study_id} more than once"
                )
            if not isinstance(rows, pd.DataFrame):
                raise TypeError(f"rows for {study_id} must be a pandas DataFrame")
            if config.input_format == "catalogue":
                missing_identity = {"dataset_id", "source_row"} - set(rows.columns)
                if missing_identity:
                    raise ValueError(
                        f"catalogue rows for {study_id} lack identities: "
                        + ", ".join(sorted(missing_identity))
                    )

            fit = fit_study(
                study_id, rows, config.model, config.settings, runner=runner
            )
            _append_fit_outputs(fit, stage_paths, schemas, config.score_chunk_size)
            processed_count += 1
            score_rows_written += fit.score_row_count
            if fit.success:
                successful_count += 1
            else:
                failed_count += 1
                failures.append(
                    StudyFailure(study_id, fit.failure_reason or "study fit failed")
                )

            study_qc = _study_qc_record(
                prepared,
                study_id,
                rows,
                fit,
                plan.study_counts[study_id],
                memory_estimates[study_id],
                config,
            )
            _append_frame(pd.DataFrame([study_qc]), qc_study_path, QC_STUDY_COLUMNS)
            if qc_smoke_path is not None:
                _append_smoke_identities(qc_smoke_path, study_id, rows, fit)
            selected_seen.add(study_id)
            del fit, rows

        if selected_seen != set(selected_study_ids):
            missing_studies = sorted(set(selected_study_ids) - selected_seen)
            raise ValueError(
                "input adapter omitted selected studies: " + ", ".join(missing_studies)
            )

        qc_stage_files = _write_run_qc(prepared, qc_stage_dir, config, plan)
        _validate_staged_outputs(stage_paths, schemas)
        _validate_tsv(qc_study_path, QC_STUDY_COLUMNS)
        if qc_smoke_path is not None:
            _validate_tsv(qc_smoke_path, QC_SMOKE_ID_COLUMNS)

        final_files: dict[Path, Path] = {
            output_paths[table]: stage_path for table, stage_path in stage_paths.items()
        }
        qc_output_dir = output_dir / "qc" / f"batch{config.index}"
        for qc_name, stage_path in qc_stage_files.items():
            final_files[qc_output_dir / qc_name] = stage_path
        final_files[qc_output_dir / qc_study_path.name] = qc_study_path
        if qc_smoke_path is not None:
            final_files[qc_output_dir / qc_smoke_path.name] = qc_smoke_path
        _publish_with_rollback(final_files, stage_dir / ".backups")

    published_paths = dict(output_paths)
    qc_output_dir = output_dir / "qc" / f"batch{config.index}"
    published_paths.update(
        {
            "study_diagnostics": qc_output_dir / "study_diagnostics.tsv",
            "run_metadata": qc_output_dir / "run_metadata.json",
        }
    )
    if (qc_output_dir / "smoke_selected_rows.tsv").exists():
        published_paths["smoke_selected_rows"] = (
            qc_output_dir / "smoke_selected_rows.tsv"
        )
    for name in ("eligibility.tsv", "exclusions.tsv"):
        if (qc_output_dir / name).exists():
            published_paths[name.removesuffix(".tsv")] = qc_output_dir / name

    return BatchSummary(
        profile=config.profile,
        model=config.model,
        batch_index=config.index,
        batch_count=plan.num_batches,
        eligible_study_count=len(plan.study_counts),
        processed_study_count=processed_count,
        successful_study_count=successful_count,
        failed_study_count=failed_count,
        eligible_rows=sum(
            plan.study_counts[study_id] for study_id in selected_study_ids
        ),
        score_rows_written=score_rows_written,
        selected_study_ids=selected_study_ids,
        assignments=plan.assignments,
        failures=tuple(failures),
        output_paths=published_paths,
        effective_settings=_effective_settings(config),
        list_batches=False,
        plan=plan,
    )


def _append_fit_outputs(
    fit: StudyFit,
    stage_paths: Mapping[str, Path],
    schemas: Mapping[str, tuple[str, ...]],
    score_chunk_size: int,
) -> None:
    for score_chunk in fit.iter_score_chunks(score_chunk_size):
        _append_frame(score_chunk, stage_paths["scores"], schemas["scores"])
    for table, frame in (
        ("densities", fit.densities),
        ("params", fit.params),
        ("posterior", fit.posterior),
        ("fit_summary", fit.fit_summary),
        ("loo", fit.loo_tbl),
    ):
        _append_frame(frame, stage_paths[table], schemas[table])


def _study_qc_record(
    prepared: Any,
    study_id: str,
    rows: pd.DataFrame,
    fit: StudyFit,
    full_eligible_rows: int,
    estimated_peak_memory_mb: int,
    config: BatchConfig,
) -> dict[str, Any]:
    source_qc = _source_study_qc(prepared, study_id)
    input_rows = _first_value(
        source_qc,
        ("raw_rows", "raw_count", "input_rows", "row_count", "n_rows"),
        full_eligible_rows,
    )
    eligible_rows = _first_value(
        source_qc, ("eligible_rows", "eligible_count"), full_eligible_rows
    )
    excluded_rows = _first_value(
        source_qc, ("excluded_rows", "excluded_count", "n_excluded"), 0
    )
    invalid_rows = _first_value(source_qc, ("invalid_rows",), 0)
    subset_excluded_rows = _first_value(source_qc, ("subset_excluded_rows",), 0)
    warnings = {name: list(value) for name, value in fit.candidate_warnings.items()}
    pareto = dict(fit.candidate_pareto_summaries)
    if not pareto:
        pareto = {
            candidate: _pareto_summary(candidate_fit.pareto_k)
            for candidate, candidate_fit in fit.candidate_fits.items()
        }
    return {
        "study_id": study_id,
        "input_rows": input_rows,
        "eligible_rows": eligible_rows,
        "selected_rows": len(rows),
        "fit_rows": len(fit.sample_row_positions),
        "excluded_rows": excluded_rows,
        "invalid_rows": invalid_rows,
        "subset_excluded_rows": subset_excluded_rows,
        "success": fit.success,
        "top_model": fit.top_model,
        "failure_reason": fit.failure_reason,
        "candidate_warnings": json.dumps(warnings, sort_keys=True),
        "pareto_k_summary": json.dumps(pareto, sort_keys=True),
        "estimated_peak_memory_mb": estimated_peak_memory_mb,
        "selected_text_bytes": _first_value(source_qc, ("selected_text_bytes",), 0),
        "profile": config.profile,
        "max_study_rows": config.max_study_rows,
        "smoke_subset_applied": bool(
            config.profile == "smoke" and len(rows) < eligible_rows
        ),
    }


def _append_smoke_identities(
    path: Path, study_id: str, rows: pd.DataFrame, fit: StudyFit
) -> None:
    if not {"dataset_id", "source_row"}.issubset(rows.columns):
        return
    selected_pairs = list(zip(rows["dataset_id"].tolist(), rows["source_row"].tolist()))
    fit_pairs: set[tuple[Any, Any]] = set()
    if len(fit.sample_row_positions):
        sample_rows = fit.eligible_rows.iloc[fit.sample_row_positions]
        if {"dataset_id", "source_row"}.issubset(sample_rows.columns):
            fit_pairs = set(zip(sample_rows["dataset_id"], sample_rows["source_row"]))
    frame = pd.DataFrame(
        [
            {
                "study_id": study_id,
                "dataset_id": dataset_id,
                "source_row": source_row,
                "fit_sample": (dataset_id, source_row) in fit_pairs,
            }
            for dataset_id, source_row in selected_pairs
        ],
        columns=QC_SMOKE_ID_COLUMNS,
    )
    _append_frame(frame, path, QC_SMOKE_ID_COLUMNS)


def _write_run_qc(
    prepared: Any, qc_stage_dir: Path, config: BatchConfig, plan: BatchPlan
) -> dict[str, Path]:
    metadata_path = qc_stage_dir / "run_metadata.json"
    metadata = {
        "profile": config.profile,
        "input_format": config.input_format,
        "study_type": config.study_type,
        "model": config.model,
        "batch_index": config.index,
        "variants_per_batch": config.variants_per_batch,
        "batch_count": plan.num_batches,
        "selected_study_ids": list(_select_studies(plan, config.index)),
        "scores_path": str(config.scores_path),
        "meta_path": str(config.meta_path),
        "output_dir": str(config.output_dir),
        "effective_settings": _json_safe(_effective_settings(config)),
        "explicit_overrides": _json_safe(dict(config.explicit_overrides)),
        "max_studies": config.max_studies,
        "max_study_rows": config.max_study_rows,
        "score_chunk_size": config.score_chunk_size,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    qc_files = {metadata_path.name: metadata_path}
    for attribute, filename in (
        ("eligibility_path", "eligibility.tsv"),
        ("exclusions_path", "exclusions.tsv"),
    ):
        source = getattr(prepared, attribute, None)
        if source is not None and Path(source).is_file():
            destination = qc_stage_dir / filename
            shutil.copyfile(source, destination)
            qc_files[filename] = destination
    return qc_files


def _publish_with_rollback(final_files: Mapping[Path, Path], backup_dir: Path) -> None:
    """Replace validated files after staging; filesystems lack group atomicity."""

    backup_dir.mkdir(parents=True, exist_ok=True)
    existing: dict[Path, Path] = {}
    for final_path in final_files:
        if final_path.exists():
            backup = backup_dir / f"{len(existing):04d}-{final_path.name}"
            shutil.copy2(final_path, backup)
            existing[final_path] = backup

    replaced: list[Path] = []
    try:
        for final_path, staged_path in final_files.items():
            final_path.parent.mkdir(parents=True, exist_ok=True)
            os.replace(staged_path, final_path)
            replaced.append(final_path)
    except OSError:
        for final_path in reversed(replaced):
            backup = existing.get(final_path)
            if backup is None:
                final_path.unlink(missing_ok=True)
            else:
                shutil.copy2(backup, final_path)
        raise


def _validate_staged_outputs(
    paths: Mapping[str, Path], schemas: Mapping[str, tuple[str, ...]]
) -> None:
    for table, path in paths.items():
        _validate_tsv(path, schemas[table])


def _validate_tsv(path: Path, expected_columns: Iterable[str]) -> None:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t", quotechar='"')
        header = next(reader, None)
        if tuple(header or ()) != tuple(expected_columns):
            raise ValueError(f"unexpected TSV header in {path}")
        width = len(tuple(expected_columns))
        for row_number, row in enumerate(reader, start=2):
            if len(row) != width:
                raise ValueError(f"invalid TSV row width at {path}:{row_number}")


def _append_frame(frame: pd.DataFrame, path: Path, columns: Iterable[str]) -> None:
    columns = tuple(columns)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    writable = frame.reindex(columns=columns).copy()
    for column in writable.columns:
        if pd.api.types.is_bool_dtype(writable[column].dtype):
            writable[column] = writable[column].astype("string").str.upper()
    writable.to_csv(
        path,
        sep="\t",
        index=False,
        mode="a",
        header=write_header,
        na_rep="NA",
        float_format="%.17g",
        quoting=csv.QUOTE_MINIMAL,
        lineterminator="\n",
    )


def _select_studies(plan: BatchPlan, index: int) -> tuple[str, ...]:
    if not plan.num_batches:
        if index == 0:
            return ()
        raise ValueError("no batches are available")
    if index == 0:
        return tuple(study_id for batch in plan.batches for study_id in batch)
    if index > plan.num_batches:
        raise ValueError(
            f"batch index {index} exceeds the {plan.num_batches} planned batches"
        )
    return plan.studies_for_batch(index)


def _effective_settings(config: BatchConfig) -> dict[str, Any]:
    return {
        **asdict(config.settings),
        "profile": config.profile,
        "max_studies": config.max_studies,
        "max_study_rows": config.max_study_rows,
        "input_chunk_size": 5000,
        "score_chunk_size": config.score_chunk_size,
    }


def _smoke_sampler_settings() -> SamplerSettings:
    return SamplerSettings(
        chains=2,
        parallel_chains=1,
        iter_warmup=200,
        iter_sampling=200,
        seed=123,
        memory_budget_mb=2048,
        max_fit_rows=500,
    )


def _resolve_path(path: Path | str) -> Path:
    if isinstance(path, str) and not path.strip():
        raise ValueError("paths cannot be empty")
    return Path(path).expanduser().resolve()


def _source_study_qc(prepared: Any, study_id: str) -> Mapping[str, Any]:
    study_qc = getattr(prepared, "study_qc", {})
    if isinstance(study_qc, Mapping):
        value = study_qc.get(study_id, {})
        return value if isinstance(value, Mapping) else {}
    return {}


def _first_value(
    source: Mapping[str, Any], names: tuple[str, ...], default: Any
) -> Any:
    for name in names:
        if name in source and source[name] is not None:
            return source[name]
    return default


def _pareto_summary(values: np.ndarray) -> dict[str, float | int | None]:
    values = np.asarray(values, dtype=float)
    if not len(values):
        return {"count": 0, "max": None, "q90": None, "over_0_7": 0}
    finite = values[np.isfinite(values)]
    if not len(finite):
        return {"count": len(values), "max": None, "q90": None, "over_0_7": 0}
    return {
        "count": int(len(values)),
        "max": float(np.max(finite)),
        "q90": float(np.quantile(finite, 0.9)),
        "over_0_7": int(np.count_nonzero(finite > 0.7)),
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value
