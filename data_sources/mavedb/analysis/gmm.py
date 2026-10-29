"""Numerical helpers and one-study workflow for the MaveDB GMM analysis.

The numerical functions in this module intentionally follow the post-fit R
helpers in ``ot_perturbation_catalog_mave_experiments/src/gmms.R``.  The
Stan likelihood has its own constants and remains the runner's responsibility.
"""

from __future__ import annotations

import hashlib
import importlib
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Iterator, Protocol, runtime_checkable

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm, skewnorm

MODE_NAMES = ("ko", "del", "neut", "pos", "noise")
MODE_CANDIDATES = ("ko", "del", "neut", "pos")
FLEXIBLE_CANDIDATES = ("flexible",)

DEFAULT_DENSITY_PARAMETERS = {
    "theta_noise": 0.01,
    "mu_ko": 0.1,
    "sigma_ko": 0.05,
    "xi_del": 0.8,
    "omega_del": 0.7,
    "alpha_del": -20.0,
    "mu_neut": 1.0,
    "sigma_neut": 0.1,
    "xi_pos": 1.2,
    "omega_pos": 0.7,
    "alpha_pos": 20.0,
}

PARAMETER_NAME_MAP = {
    "frac_non_neut[1]": "frac_ko",
    "frac_non_neut[2]": "frac_del",
    "frac_non_neut[3]": "frac_pos",
    "y_median": "median",
    "y_mad": "mad",
}

PARAMETER_COLUMNS = (
    "frac_neut",
    "frac_ko",
    "frac_del",
    "frac_pos",
    "sigma_neut",
    "beta",
    "delta",
    "theta_ko",
    "theta_del",
    "theta_neut",
    "theta_pos",
    "lp__",
    "theta_noise",
    "median",
    "mad",
)

POSTERIOR_COLUMNS = (
    "iteration",
    "chain",
    "draw",
    "lp__",
    "frac_neut",
    "frac_ko",
    "frac_del",
    "frac_pos",
    "theta_noise",
    "theta_ko",
    "theta_del",
    "theta_neut",
    "theta_pos",
    "beta",
    "delta",
    "sigma_neut",
)

FIT_SUMMARY_COLUMNS = (
    "variable",
    "mean",
    "median",
    "sd",
    "mad",
    "q5",
    "q95",
    "rhat",
    "ess_bulk",
    "ess_tail",
)

LOO_COLUMNS = (
    "elpd_diff",
    "se_diff",
    "elpd_loo",
    "se_elpd_loo",
    "p_loo",
    "se_p_loo",
    "looic",
    "se_looic",
)

_SCORE_BASE_COLUMNS = (
    "score_urn",
    "hgvs_nt",
    "hgvs_pro",
    "s",
    "w",
    "mode_density",
    "mode_sd",
    "model",
    "top_model",
)
_SCORE_IDENTITY_COLUMNS = ("dataset_id", "source_row")
_DENSITY_FLEXIBLE_COLUMNS = (
    "score_urn",
    "mode",
    "name",
    "x",
    "density",
    "scaled_density",
    "model",
    "top_model",
)
_DENSITY_MODES_COLUMNS = (
    "score_urn",
    "model",
    "top_model",
    "mode",
    "name",
    "x",
    "density",
    "scaled_density",
)
_PARAMETER_OUTPUT_COLUMNS = {
    "flexible": (
        "lp__",
        "frac_neut",
        "frac_ko",
        "frac_del",
        "frac_pos",
        "sigma_neut",
        "beta",
        "delta",
        "theta_noise",
        "theta_ko",
        "theta_del",
        "theta_neut",
        "theta_pos",
        "median",
        "mad",
    ),
    "modes": (
        "lp__",
        "frac_neut",
        "sigma_neut",
        "beta",
        "delta",
        "theta_noise",
        "theta_ko",
        "theta_del",
        "theta_neut",
        "theta_pos",
        "median",
        "mad",
    ),
}
_POSTERIOR_FLEXIBLE_COLUMNS = (
    "iteration",
    "chain",
    "draw",
    "lp__",
    "frac_neut",
    "frac_ko",
    "frac_del",
    "frac_pos",
    "theta_noise",
    "theta_ko",
    "theta_del",
    "theta_neut",
    "theta_pos",
    "beta",
    "delta",
    "sigma_neut",
)
_POSTERIOR_MODES_COLUMNS = (
    "iteration",
    "chain",
    "draw",
    "lp__",
    "frac_neut",
    "theta_noise",
    "theta_ko",
    "theta_del",
    "theta_neut",
    "theta_pos",
    "beta",
    "delta",
    "sigma_neut",
)
_FIT_SUMMARY_OUTPUT_COLUMNS = (
    "score_urn",
    "model",
    "top_model",
    *FIT_SUMMARY_COLUMNS,
)
_LOO_OUTPUT_COLUMNS = ("score_urn", "model", *LOO_COLUMNS)


def family_schemas(
    model: str, identity_columns: Iterable[str] = ()
) -> dict[str, tuple[str, ...]]:
    """Return frozen output headers for a model family.

    ``score_urn`` is the legacy compatibility key. Catalogue writers may add
    ``dataset_id`` and ``source_row`` to score outputs without changing the
    legacy headers when those columns are absent.
    """

    _validate_model(model)
    available_identities = set(identity_columns)
    identities = tuple(
        column for column in _SCORE_IDENTITY_COLUMNS if column in available_identities
    )
    if model == "flexible":
        score_columns = (*_SCORE_BASE_COLUMNS, *identities)
    else:
        score_columns = (
            "model",
            "top_model",
            "score_urn",
            "hgvs_nt",
            "hgvs_pro",
            "s",
            "w",
            "mode_density",
            "mode_sd",
            *identities,
        )
    parameter_columns = _PARAMETER_OUTPUT_COLUMNS[model]
    posterior_columns = (
        _POSTERIOR_FLEXIBLE_COLUMNS if model == "flexible" else _POSTERIOR_MODES_COLUMNS
    )
    return {
        "scores": score_columns,
        "densities": (
            _DENSITY_FLEXIBLE_COLUMNS if model == "flexible" else _DENSITY_MODES_COLUMNS
        ),
        "params": ("score_urn", "model", "top_model", *parameter_columns),
        "posterior": ("score_urn", "model", "top_model", *posterior_columns),
        "fit_summary": _FIT_SUMMARY_OUTPUT_COLUMNS,
        "loo": _LOO_OUTPUT_COLUMNS,
    }


def g(x: Any, median: float, mad: float) -> Any:
    """Standardize a score with the sampled median and scaled MAD."""

    values = np.asarray(x, dtype=float)
    result = (values - median) / mad
    return _scalar_or_array(result, x)


def g_inv(x: Any, median: float, mad: float) -> Any:
    """Invert :func:`g`."""

    values = np.asarray(x, dtype=float)
    result = values * mad + median
    return _scalar_or_array(result, x)


def f(x: Any, beta: float, delta: float) -> Any:
    """Map latent fitness ``w`` back to the transformed score scale."""

    values = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        result = -(np.log(2.0 / values - 1.0)) / beta - delta
    return _scalar_or_array(result, x)


def f_inv(x: Any, beta: float, delta: float) -> Any:
    """Map transformed scores to latent fitness using a stable logistic."""

    values = np.asarray(x, dtype=float)
    result = 2.0 * expit(beta * (values + delta))
    return _scalar_or_array(result, x)


def dflexgmm(
    x: Any,
    theta_ko: float,
    theta_del: float,
    theta_neut: float,
    theta_pos: float,
    theta_noise: float = DEFAULT_DENSITY_PARAMETERS["theta_noise"],
    mu_ko: float = DEFAULT_DENSITY_PARAMETERS["mu_ko"],
    sigma_ko: float = DEFAULT_DENSITY_PARAMETERS["sigma_ko"],
    xi_del: float = DEFAULT_DENSITY_PARAMETERS["xi_del"],
    omega_del: float = DEFAULT_DENSITY_PARAMETERS["omega_del"],
    alpha_del: float = DEFAULT_DENSITY_PARAMETERS["alpha_del"],
    mu_neut: float = DEFAULT_DENSITY_PARAMETERS["mu_neut"],
    sigma_neut: float = DEFAULT_DENSITY_PARAMETERS["sigma_neut"],
    xi_pos: float = DEFAULT_DENSITY_PARAMETERS["xi_pos"],
    omega_pos: float = DEFAULT_DENSITY_PARAMETERS["omega_pos"],
    alpha_pos: float = DEFAULT_DENSITY_PARAMETERS["alpha_pos"],
    name: str = "Flex GMM",
) -> pd.DataFrame:
    """Return weighted component densities in R's fixed mode order."""

    values = np.asarray(x, dtype=float).reshape(-1)
    component_densities = (
        norm.pdf(values, loc=mu_ko, scale=sigma_ko) * theta_ko,
        skewnorm.pdf(values, a=alpha_del, loc=xi_del, scale=omega_del) * theta_del,
        norm.pdf(values, loc=mu_neut, scale=sigma_neut) * theta_neut,
        skewnorm.pdf(values, a=alpha_pos, loc=xi_pos, scale=omega_pos) * theta_pos,
        np.where((values >= 0.0) & (values <= 2.0), 0.5 * theta_noise, 0.0),
    )
    return pd.concat(
        [
            pd.DataFrame(
                {
                    "mode": mode,
                    "name": name,
                    "x": values,
                    "density": density,
                }
            )
            for mode, density in zip(MODE_NAMES, component_densities, strict=True)
        ],
        ignore_index=True,
    )


def transformed_density(
    beta: float,
    delta: float,
    median: float,
    mad: float,
    theta_ko: float,
    theta_del: float,
    theta_neut: float,
    theta_pos: float,
    theta_noise: float = DEFAULT_DENSITY_PARAMETERS["theta_noise"],
    mu_ko: float = DEFAULT_DENSITY_PARAMETERS["mu_ko"],
    sigma_ko: float = DEFAULT_DENSITY_PARAMETERS["sigma_ko"],
    xi_del: float = DEFAULT_DENSITY_PARAMETERS["xi_del"],
    omega_del: float = DEFAULT_DENSITY_PARAMETERS["omega_del"],
    alpha_del: float = DEFAULT_DENSITY_PARAMETERS["alpha_del"],
    mu_neut: float = DEFAULT_DENSITY_PARAMETERS["mu_neut"],
    sigma_neut: float = DEFAULT_DENSITY_PARAMETERS["sigma_neut"],
    xi_pos: float = DEFAULT_DENSITY_PARAMETERS["xi_pos"],
    omega_pos: float = DEFAULT_DENSITY_PARAMETERS["omega_pos"],
    alpha_pos: float = DEFAULT_DENSITY_PARAMETERS["alpha_pos"],
    length_out: int = 500,
) -> pd.DataFrame:
    """Build density curves on latent-fitness and score axes.

    The score-axis ordinates deliberately remain unchanged: the legacy R
    helper does not apply a change-of-variable Jacobian.
    """

    if (
        not isinstance(length_out, int)
        or isinstance(length_out, bool)
        or length_out < 1
    ):
        raise ValueError("length_out must be a positive integer")
    w = np.linspace(0.001, 2.0, length_out)
    density_w = dflexgmm(
        w,
        theta_ko,
        theta_del,
        theta_neut,
        theta_pos,
        theta_noise,
        mu_ko,
        sigma_ko,
        xi_del,
        omega_del,
        alpha_del,
        mu_neut,
        sigma_neut,
        xi_pos,
        omega_pos,
        alpha_pos,
        name="w",
    )
    density_s = density_w.loc[(density_w["x"] > 0.03) & (density_w["x"] < 1.97)].copy()
    density_s["x"] = f(density_s["x"].to_numpy(), beta, delta) * mad + median
    density_s["name"] = "s"
    return pd.concat([density_w, density_s], ignore_index=True)[
        ["mode", "name", "x", "density"]
    ]


def transformed_scores(
    scores: pd.DataFrame, params: Mapping[str, Any] | pd.Series | pd.DataFrame
) -> pd.DataFrame:
    """Return score rows with transformed fitness and both classifications."""

    if "s" not in scores.columns:
        raise ValueError("scores must contain the raw score column 's'")
    values = pd.to_numeric(scores["s"], errors="coerce").to_numpy(dtype=float)
    median = float(_param_value(params, "median"))
    mad = float(_param_value(params, "mad"))
    beta = float(_param_value(params, "beta"))
    delta = float(_param_value(params, "delta"))
    w = np.asarray(f_inv(g(values, median, mad), beta, delta), dtype=float)
    out = scores.copy()
    out["w"] = w
    out["mode_density"] = classify_point_density(
        w,
        _param_value(params, "theta_ko"),
        _param_value(params, "theta_del"),
        _param_value(params, "theta_neut"),
        _param_value(params, "theta_pos"),
        _param_value(params, "sigma_neut"),
    )
    out["mode_sd"] = classify_point_sd(w, _param_value(params, "sigma_neut"), 0.99)
    return out


def classify_point_density(
    x: Any,
    theta_ko: float,
    theta_del: float,
    theta_neut: float,
    theta_pos: float,
    sigma_neut: float,
) -> np.ndarray:
    """Choose the largest weighted density; ties use R's mode order."""

    values = np.asarray(x, dtype=float).reshape(-1)
    densities = dflexgmm(
        values,
        theta_ko,
        theta_del,
        theta_neut,
        theta_pos,
        sigma_neut=sigma_neut,
    )
    density_matrix = densities["density"].to_numpy().reshape(len(MODE_NAMES), -1).T
    labels: list[str | None] = []
    for row in density_matrix:
        valid = ~np.isnan(row)
        if not valid.any():
            labels.append(None)
        else:
            labels.append(MODE_NAMES[int(np.argmax(np.where(valid, row, -np.inf)))])
    return np.asarray(labels, dtype=object)


def classify_point_sd(x: Any, sigma_neut: float, threshold: float = 0.95) -> np.ndarray:
    """Classify by distance from neutral with the historical KO boundary."""

    values = np.asarray(x, dtype=float).reshape(-1)
    t_neut = norm.ppf(threshold, loc=0.0, scale=sigma_neut)
    t_ko = norm.ppf(threshold, loc=0.1, scale=0.05)
    out = np.full(values.shape, "neut", dtype=object)
    out[values > 1.0 + t_neut] = "pos"
    out[values < 1.0 - t_neut] = "del"
    out[values < 0.1 + t_ko] = "ko"
    return out


@dataclass(frozen=True)
class SamplerSettings:
    """Sampler and study-memory settings shared with the Stan adapter."""

    chains: int = 4
    parallel_chains: int = 4
    iter_warmup: int = 1000
    iter_sampling: int = 2000
    adapt_delta: float = 0.9
    max_treedepth: int = 10
    refresh: int = 4000
    memory_budget_mb: int = 16384
    seed: int | None = None
    max_fit_rows: int = 256_000

    def __post_init__(self) -> None:
        positive_integer_fields = (
            "chains",
            "parallel_chains",
            "iter_warmup",
            "iter_sampling",
            "max_treedepth",
            "memory_budget_mb",
            "max_fit_rows",
        )
        for field_name in positive_integer_fields:
            value = getattr(self, field_name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise ValueError(f"{field_name} must be a positive integer")
        if (
            not isinstance(self.refresh, int)
            or isinstance(self.refresh, bool)
            or self.refresh < 0
        ):
            raise ValueError("refresh must be a nonnegative integer")
        if self.parallel_chains > self.chains:
            raise ValueError("parallel_chains cannot exceed chains")
        if not math.isfinite(self.adapt_delta) or not 0.0 < self.adapt_delta < 1.0:
            raise ValueError("adapt_delta must be finite and between 0 and 1")
        if self.seed is not None and (
            not isinstance(self.seed, int)
            or isinstance(self.seed, bool)
            or not 0 <= self.seed <= 2**32 - 1
        ):
            raise ValueError("seed must be None or an integer from 0 to 2**32 - 1")


@dataclass(frozen=True)
class BatchPlan:
    """Deterministic one-based whole-study batch assignment."""

    variants_per_batch: int
    batches: tuple[tuple[str, ...], ...]
    batch_counts: tuple[int, ...]
    study_counts: Mapping[str, int]
    study_to_batch: Mapping[str, int]

    @property
    def num_batches(self) -> int:
        return len(self.batches)

    @property
    def assignments(self) -> tuple[tuple[int, str, int], ...]:
        """Return ``(batch, study_id, row_count)`` tuples in planned order."""

        return tuple(
            (batch_index, study_id, self.study_counts[study_id])
            for batch_index, studies in enumerate(self.batches, start=1)
            for study_id in studies
        )

    def studies_for_batch(self, index: int) -> tuple[str, ...]:
        if not isinstance(index, int) or isinstance(index, bool):
            raise ValueError("batch index must be an integer")
        if not 1 <= index <= self.num_batches:
            raise ValueError(f"batch index must be between 1 and {self.num_batches}")
        return self.batches[index - 1]


@dataclass(frozen=True)
class StudyFitError(Exception):
    """Expected data or numerical failure that produces placeholder outputs."""

    message: str

    def __str__(self) -> str:
        return self.message


@runtime_checkable
class CandidateResult(Protocol):
    """Structural view of ``_gmm_stan.CandidateFit`` without eager imports."""

    candidate: str
    posterior: pd.DataFrame
    fit_summary: pd.DataFrame
    params: Mapping[str, float]
    loo: Mapping[str, float]
    pointwise: np.ndarray
    pareto_k: np.ndarray
    warnings: list[str]


Runner = Any


@dataclass
class StudyFit:
    """Compact per-study result with score rows generated in bounded chunks."""

    study_id: str
    model: str
    settings: SamplerSettings
    success: bool
    top_model: str | None
    eligible_rows: pd.DataFrame = field(repr=False)
    sample_row_positions: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=int)
    )
    candidate_fits: Mapping[str, CandidateResult] = field(
        default_factory=dict, repr=False
    )
    params: pd.DataFrame = field(default_factory=pd.DataFrame)
    posterior: pd.DataFrame = field(default_factory=pd.DataFrame)
    fit_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    densities: pd.DataFrame = field(default_factory=pd.DataFrame)
    loo_tbl: pd.DataFrame = field(default_factory=pd.DataFrame)
    failure_reason: str | None = None
    candidate_warnings: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    candidate_pareto_summaries: Mapping[str, Mapping[str, Any]] = field(
        default_factory=dict
    )

    @property
    def score_row_count(self) -> int:
        candidate_count = len(
            FLEXIBLE_CANDIDATES if self.model == "flexible" else MODE_CANDIDATES
        )
        return (
            len(self.eligible_rows) * candidate_count
            if self.success
            else len(self.eligible_rows)
        )

    def iter_score_chunks(self, chunk_size: int = 50_000) -> Iterator[pd.DataFrame]:
        """Yield final score outputs without materializing all candidates."""

        if (
            not isinstance(chunk_size, int)
            or isinstance(chunk_size, bool)
            or chunk_size < 1
        ):
            raise ValueError("chunk_size must be a positive integer")
        schemas = family_schemas(self.model, self.eligible_rows.columns)
        score_columns = schemas["scores"]
        if not self.success:
            for start in range(0, len(self.eligible_rows), chunk_size):
                rows = self.eligible_rows.iloc[start : start + chunk_size]
                chunk = _base_score_rows(rows, self.study_id)
                chunk["w"] = np.nan
                chunk["mode_density"] = np.nan
                chunk["mode_sd"] = np.nan
                chunk["model"] = np.nan
                chunk["top_model"] = pd.array([pd.NA] * len(chunk), dtype="boolean")
                yield chunk.reindex(columns=score_columns)
            return

        output_candidates = _output_candidate_order(self.model, self.top_model)
        for candidate in output_candidates:
            candidate_fit = self.candidate_fits[candidate]
            params = _normalise_params(candidate_fit.params)
            for start in range(0, len(self.eligible_rows), chunk_size):
                rows = self.eligible_rows.iloc[start : start + chunk_size]
                chunk = _base_score_rows(rows, self.study_id)
                transformed = transformed_scores(chunk, params)
                transformed["model"] = candidate
                transformed["top_model"] = candidate == self.top_model
                yield transformed.reindex(columns=score_columns)


def plan_batches(
    study_counts: Mapping[str, int] | Iterable[tuple[str, int]],
    variants_per_batch: int,
) -> BatchPlan:
    """Plan whole-study batches by decreasing count, then full ID.

    A study larger than the target starts and occupies its own batch. Batch
    numbers are one-based; an empty plan contains no batch zero.
    """

    if (
        not isinstance(variants_per_batch, int)
        or isinstance(variants_per_batch, bool)
        or variants_per_batch <= 0
    ):
        raise ValueError("variants_per_batch must be a positive integer")
    items = list(
        study_counts.items() if isinstance(study_counts, Mapping) else study_counts
    )
    normalized: dict[str, int] = {}
    for study_id, count in items:
        study_key = str(study_id)
        if not study_key:
            raise ValueError("study IDs must be nonempty")
        if study_key in normalized:
            raise ValueError(f"duplicate study ID: {study_key}")
        if (
            not isinstance(count, (int, np.integer))
            or isinstance(count, bool)
            or count < 0
        ):
            raise ValueError(f"row count for {study_key} must be a nonnegative integer")
        normalized[study_key] = int(count)

    ordered = sorted(normalized, key=lambda study_id: (-normalized[study_id], study_id))
    batches: list[list[str]] = []
    batch_counts: list[int] = []
    current: list[str] = []
    current_count = 0

    def close_current() -> None:
        nonlocal current, current_count
        if current:
            batches.append(current)
            batch_counts.append(current_count)
            current = []
            current_count = 0

    for study_id in ordered:
        count = normalized[study_id]
        if count > variants_per_batch:
            close_current()
            batches.append([study_id])
            batch_counts.append(count)
            continue
        if current and current_count + count > variants_per_batch:
            close_current()
        current.append(study_id)
        current_count += count
    close_current()
    study_to_batch = {
        study_id: batch_index
        for batch_index, studies in enumerate(batches, start=1)
        for study_id in studies
    }
    return BatchPlan(
        variants_per_batch=variants_per_batch,
        batches=tuple(tuple(studies) for studies in batches),
        batch_counts=tuple(batch_counts),
        study_counts=normalized,
        study_to_batch=study_to_batch,
    )


def fit_study(
    study_id: str,
    rows: pd.DataFrame,
    model: str,
    settings: SamplerSettings,
    runner: Runner | None = None,
) -> StudyFit:
    """Fit one score set, compare candidates, and retain chunkable outputs.

    Fitting and transformed scores use finite raw scores with nonmissing
    protein HGVS. A single bounded sample and its row order are shared by all
    candidates. Only explicitly expected runner/numerical errors become
    placeholder results; programming errors propagate to the caller.
    """

    _validate_model(model)
    if not isinstance(settings, SamplerSettings):
        raise TypeError("settings must be a SamplerSettings instance")
    if not isinstance(rows, pd.DataFrame):
        raise TypeError("rows must be a pandas DataFrame")
    score_column = (
        "score" if "score" in rows.columns else "s" if "s" in rows.columns else None
    )
    if score_column is None:
        raise ValueError("rows must contain a numeric 'score' or 's' column")
    if "hgvs_pro" not in rows.columns:
        raise ValueError("rows must contain the protein HGVS column 'hgvs_pro'")

    eligible = _select_eligible_rows(rows, score_column)
    if "s" not in eligible.columns:
        eligible["s"] = pd.to_numeric(eligible[score_column], errors="coerce").to_numpy(
            dtype=float
        )
    elif score_column != "s":
        eligible["s"] = pd.to_numeric(eligible[score_column], errors="coerce").to_numpy(
            dtype=float
        )

    candidate_names = FLEXIBLE_CANDIDATES if model == "flexible" else MODE_CANDIDATES
    from ._gmm_stan import SamplerError

    sample_positions = np.array([], dtype=int)
    candidate_warnings: dict[str, tuple[str, ...]] = {}
    candidate_pareto_summaries: dict[str, Mapping[str, Any]] = {}
    try:
        if eligible.empty:
            raise StudyFitError("no rows have both finite scores and protein HGVS")
        sample_positions, sample_scores = _sample_fit_rows(
            eligible,
            settings.max_fit_rows,
            _stable_seed(settings.seed, study_id, "__sample__"),
        )
        scaled_mad = _scaled_mad(sample_scores)
        if not math.isfinite(scaled_mad) or scaled_mad == 0.0:
            raise StudyFitError("sampled scores have zero or nonfinite scaled MAD")

        active_runner = runner if runner is not None else _default_runner()
        candidate_fits: dict[str, CandidateResult] = {}
        for candidate in candidate_names:
            seed = _stable_seed(settings.seed, study_id, candidate)
            try:
                result = active_runner(
                    sample_scores.copy(), candidate, settings, seed=seed
                )
            except (
                StudyFitError,
                SamplerError,
                FloatingPointError,
                OverflowError,
                np.linalg.LinAlgError,
            ) as error:
                failure_candidate = getattr(error, "candidate", None) or candidate
                failure_warnings = tuple(getattr(error, "warnings", ()))
                if failure_warnings:
                    prior_warnings = candidate_warnings.get(failure_candidate, ())
                    candidate_warnings[failure_candidate] = tuple(
                        dict.fromkeys((*prior_warnings, *failure_warnings))
                    )
                raise StudyFitError(
                    f"{failure_candidate} candidate failed: "
                    f"{type(error).__name__}: {error}"
                ) from error
            if not isinstance(result, CandidateResult):
                raise TypeError(
                    f"runner returned an invalid result for candidate {candidate}"
                )
            if result.candidate != candidate:
                raise ValueError(
                    f"runner returned candidate {result.candidate!r} for requested {candidate!r}"
                )
            candidate_fits[candidate] = result
            candidate_warnings[candidate] = tuple(result.warnings)
            candidate_pareto_summaries[candidate] = _pareto_k_summary(result.pareto_k)

        if model == "flexible":
            top_model = "flexible"
            loo_tbl = _flexible_loo_table(
                candidate_fits["flexible"], expected_observations=len(sample_scores)
            )
        else:
            loo_tbl, top_model = _compare_mode_candidates(
                candidate_fits, expected_observations=len(sample_scores)
            )

        return _successful_study_fit(
            study_id,
            model,
            settings,
            eligible,
            sample_positions,
            candidate_fits,
            top_model,
            loo_tbl,
            candidate_warnings,
            candidate_pareto_summaries,
        )
    except (
        StudyFitError,
        SamplerError,
        FloatingPointError,
        OverflowError,
        np.linalg.LinAlgError,
    ) as error:
        return _failed_study_fit(
            study_id,
            model,
            settings,
            eligible,
            sample_positions,
            candidate_warnings,
            candidate_pareto_summaries,
            f"{type(error).__name__}: {error}",
        )


def _successful_study_fit(
    study_id: str,
    model: str,
    settings: SamplerSettings,
    eligible_rows: pd.DataFrame,
    sample_row_positions: np.ndarray,
    candidate_fits: Mapping[str, CandidateResult],
    top_model: str,
    loo_tbl: pd.DataFrame,
    candidate_warnings: Mapping[str, tuple[str, ...]],
    candidate_pareto_summaries: Mapping[str, Mapping[str, Any]],
) -> StudyFit:
    schemas = family_schemas(model)
    output_candidates = _output_candidate_order(model, top_model)
    params_rows: list[dict[str, Any]] = []
    posterior_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []
    density_frames: list[pd.DataFrame] = []

    for candidate in output_candidates:
        result = candidate_fits[candidate]
        is_top = candidate == top_model
        params = _normalise_params(result.params)
        params_rows.append(
            {
                "score_urn": study_id,
                "model": candidate,
                "top_model": is_top,
                **{column: params.get(column, np.nan) for column in PARAMETER_COLUMNS},
            }
        )

        posterior = _normalise_posterior(result.posterior, candidate)
        posterior.insert(0, "top_model", is_top)
        posterior.insert(0, "model", candidate)
        posterior.insert(0, "score_urn", study_id)
        posterior_frames.append(posterior.reindex(columns=schemas["posterior"]))

        fit_summary = _normalise_fit_summary(result.fit_summary)
        fit_summary.insert(0, "top_model", is_top)
        fit_summary.insert(0, "model", candidate)
        fit_summary.insert(0, "score_urn", study_id)
        summary_frames.append(fit_summary.reindex(columns=schemas["fit_summary"]))

        density = _transformed_density_for_params(params)
        if not np.isfinite(density[["x", "density"]].to_numpy(dtype=float)).all():
            raise StudyFitError(f"candidate {candidate} produced nonfinite densities")
        if model == "flexible":
            density["scaled_density"] = density["density"] / 100.0
        else:
            max_density = float(density["density"].max())
            if not math.isfinite(max_density) or max_density == 0.0:
                raise StudyFitError(
                    f"candidate {candidate} has no finite positive density"
                )
            density["scaled_density"] = density["density"] / max_density
        density["model"] = candidate
        density["top_model"] = is_top
        density.insert(0, "score_urn", study_id)
        density_frames.append(density.reindex(columns=schemas["densities"]))

    params_frame = pd.DataFrame(params_rows, columns=schemas["params"])
    posterior_frame = pd.concat(posterior_frames, ignore_index=True)
    summary_frame = pd.concat(summary_frames, ignore_index=True)
    density_frame = pd.concat(density_frames, ignore_index=True)
    loo_output = loo_tbl.copy()
    loo_output.insert(0, "score_urn", study_id)
    loo_output = loo_output.reindex(columns=schemas["loo"])
    return StudyFit(
        study_id=study_id,
        model=model,
        settings=settings,
        success=True,
        top_model=top_model,
        eligible_rows=eligible_rows,
        sample_row_positions=sample_row_positions,
        candidate_fits=dict(candidate_fits),
        params=params_frame,
        posterior=posterior_frame,
        fit_summary=summary_frame,
        densities=density_frame,
        loo_tbl=loo_output,
        candidate_warnings=candidate_warnings,
        candidate_pareto_summaries=dict(candidate_pareto_summaries),
    )


def _failed_study_fit(
    study_id: str,
    model: str,
    settings: SamplerSettings,
    eligible_rows: pd.DataFrame,
    sample_row_positions: np.ndarray,
    candidate_warnings: Mapping[str, tuple[str, ...]],
    candidate_pareto_summaries: Mapping[str, Mapping[str, Any]],
    failure_reason: str,
) -> StudyFit:
    schemas = family_schemas(model, eligible_rows.columns)
    return StudyFit(
        study_id=study_id,
        model=model,
        settings=settings,
        success=False,
        top_model=None,
        eligible_rows=eligible_rows,
        sample_row_positions=sample_row_positions,
        candidate_fits={},
        params=_study_sentinel_frame(schemas["params"], study_id),
        posterior=_study_sentinel_frame(schemas["posterior"], study_id),
        fit_summary=_study_sentinel_frame(schemas["fit_summary"], study_id),
        densities=_study_sentinel_frame(schemas["densities"], study_id),
        loo_tbl=_study_sentinel_frame(schemas["loo"], study_id),
        failure_reason=failure_reason,
        candidate_warnings=dict(candidate_warnings),
        candidate_pareto_summaries=dict(candidate_pareto_summaries),
    )


def _compare_mode_candidates(
    candidates: Mapping[str, CandidateResult],
    expected_observations: int,
) -> tuple[pd.DataFrame, str]:
    comparison_rows: list[dict[str, Any]] = []
    comparison_order = ("ko", "del", "neut")
    for candidate in comparison_order:
        result = candidates[candidate]
        _validate_loo_result(result, candidate, expected_observations)

    original_winner = max(
        comparison_order,
        key=lambda candidate: (
            float(candidates[candidate].loo["elpd_loo"]),
            -comparison_order.index(candidate),
        ),
    )
    best_elpd = float(candidates[original_winner].loo["elpd_loo"])
    best_pointwise = np.asarray(candidates[original_winner].pointwise, dtype=float)

    for candidate in comparison_order:
        result = candidates[candidate]
        pointwise = np.asarray(result.pointwise, dtype=float)
        difference = pointwise - best_pointwise
        se_diff = (
            0.0 if candidate == original_winner else _standard_error_of_sum(difference)
        )
        comparison_rows.append(
            {
                "model": candidate,
                "elpd_diff": float(result.loo["elpd_loo"]) - best_elpd,
                "se_diff": se_diff,
                **_loo_metrics(result.loo),
            }
        )
    comparison_rows.sort(
        key=lambda row: (
            -float(row["elpd_loo"]),
            comparison_order.index(str(row["model"])),
        )
    )
    retained = [row for row in comparison_rows if float(row["p_loo"]) < 12.0]
    if not retained:
        raise StudyFitError(
            "all ko/del/neut candidates were rejected by the p_loo < 12 filter"
        )
    top_model = str(retained[0]["model"])
    return pd.DataFrame(retained, columns=("model", *LOO_COLUMNS)), top_model


def _flexible_loo_table(
    candidate: CandidateResult, expected_observations: int
) -> pd.DataFrame:
    _validate_loo_result(candidate, "flexible", expected_observations)
    return pd.DataFrame(
        [
            {
                "model": "flexible",
                "elpd_diff": np.nan,
                "se_diff": np.nan,
                **_loo_metrics(candidate.loo),
            }
        ],
        columns=("model", *LOO_COLUMNS),
    )


def _validate_loo_result(
    candidate: CandidateResult, candidate_name: str, expected_observations: int
) -> None:
    metric_names = ("elpd_loo", "se_elpd_loo", "p_loo", "se_p_loo", "looic", "se_looic")
    missing = [name for name in metric_names if name not in candidate.loo]
    if missing:
        raise ValueError(
            f"candidate {candidate_name} LOO metrics are missing: {', '.join(missing)}"
        )
    for name in metric_names:
        if not math.isfinite(float(candidate.loo[name])):
            raise StudyFitError(
                f"candidate {candidate_name} has nonfinite LOO metric {name}"
            )
    pointwise = np.asarray(candidate.pointwise, dtype=float)
    if (
        pointwise.ndim != 1
        or len(pointwise) != expected_observations
        or not np.isfinite(pointwise).all()
    ):
        raise StudyFitError(
            f"candidate {candidate_name} has invalid pointwise LOO values"
        )
    if len(candidate.pareto_k) != len(pointwise):
        raise ValueError(
            f"candidate {candidate_name} Pareto-k length does not match pointwise LOO"
        )


def _loo_metrics(metrics: Mapping[str, float]) -> dict[str, float]:
    return {name: float(metrics[name]) for name in LOO_COLUMNS if name in metrics}


def _standard_error_of_sum(values: np.ndarray) -> float:
    if len(values) <= 1:
        return 0.0
    return float(np.sqrt(len(values) * np.var(values, ddof=1)))


def _pareto_k_summary(values: np.ndarray) -> dict[str, float | int | None]:
    pareto_k = np.asarray(values, dtype=float)
    finite = pareto_k[np.isfinite(pareto_k)]
    if not len(finite):
        return {
            "count": int(len(pareto_k)),
            "max": None,
            "q90": None,
            "over_0_7": 0,
        }
    return {
        "count": int(len(pareto_k)),
        "max": float(np.max(finite)),
        "q90": float(np.quantile(finite, 0.9)),
        "over_0_7": int(np.count_nonzero(finite > 0.7)),
    }


def _transformed_density_for_params(params: Mapping[str, Any]) -> pd.DataFrame:
    density_parameters = {
        name: params.get(name, DEFAULT_DENSITY_PARAMETERS.get(name))
        for name in (
            "beta",
            "delta",
            "median",
            "mad",
            "theta_ko",
            "theta_del",
            "theta_neut",
            "theta_pos",
            "mu_ko",
            "sigma_ko",
            "xi_del",
            "omega_del",
            "alpha_del",
            "mu_neut",
            "sigma_neut",
            "xi_pos",
            "omega_pos",
            "alpha_pos",
        )
    }
    missing = [name for name, value in density_parameters.items() if value is None]
    if missing:
        raise ValueError(
            f"candidate params are missing density values: {', '.join(missing)}"
        )
    return transformed_density(**density_parameters)


def _normalise_params(params: Mapping[str, Any]) -> dict[str, Any]:
    normalized = {
        PARAMETER_NAME_MAP.get(name, name): value for name, value in params.items()
    }
    return normalized


def _normalise_posterior(posterior: pd.DataFrame, candidate: str) -> pd.DataFrame:
    if not isinstance(posterior, pd.DataFrame):
        raise TypeError(f"candidate {candidate} posterior must be a pandas DataFrame")
    normalized = posterior.rename(columns=PARAMETER_NAME_MAP).copy()
    normalized = normalized.rename(
        columns={".iteration": "iteration", ".chain": "chain", ".draw": "draw"}
    )
    unknown = set(normalized.columns) - set(POSTERIOR_COLUMNS)
    if unknown:
        raise ValueError(
            f"candidate {candidate} posterior has unsupported columns: {', '.join(sorted(unknown))}"
        )
    if candidate == "neut" and "frac_neut" in normalized.columns:
        normalized["frac_neut"] = np.nan
    return normalized.reindex(columns=POSTERIOR_COLUMNS)


def _normalise_fit_summary(fit_summary: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(fit_summary, pd.DataFrame):
        raise TypeError("candidate fit_summary must be a pandas DataFrame")
    normalized = fit_summary.rename(columns=PARAMETER_NAME_MAP).copy()
    if "variable" in normalized.columns:
        normalized["variable"] = normalized["variable"].map(
            lambda name: PARAMETER_NAME_MAP.get(name, name)
        )
    unknown = set(normalized.columns) - set(FIT_SUMMARY_COLUMNS)
    if unknown:
        raise ValueError(
            "candidate fit_summary has unsupported columns: "
            + ", ".join(sorted(unknown))
        )
    return normalized.reindex(columns=FIT_SUMMARY_COLUMNS)


def _select_eligible_rows(rows: pd.DataFrame, score_column: str) -> pd.DataFrame:
    numeric_scores = pd.to_numeric(rows[score_column], errors="coerce").to_numpy(
        dtype=float
    )
    finite_score = np.isfinite(numeric_scores)
    protein_present = rows["hgvs_pro"].notna().to_numpy()
    selected = rows.iloc[np.flatnonzero(finite_score & protein_present)].copy()
    return selected


def _sample_fit_rows(
    eligible: pd.DataFrame, max_fit_rows: int, seed: int | None
) -> tuple[np.ndarray, np.ndarray]:
    if len(eligible) <= max_fit_rows:
        positions = np.arange(len(eligible), dtype=int)
    else:
        rng = np.random.default_rng(seed)
        positions = np.sort(rng.choice(len(eligible), size=max_fit_rows, replace=False))
    fit_rows = eligible.iloc[positions]
    scores = pd.to_numeric(fit_rows["s"], errors="coerce").to_numpy(dtype=float)
    if not len(scores) or not np.isfinite(scores).all():
        raise StudyFitError("sampled score vector is empty or nonfinite")
    return positions, scores


def _scaled_mad(scores: np.ndarray) -> float:
    median = float(np.median(scores))
    return float(np.median(np.abs(scores - median)) * 1.4826)


def _stable_seed(base_seed: int | None, study_id: str, name: str) -> int | None:
    if base_seed is None:
        return None
    payload = f"mavedb-gmm\0{base_seed}\0{study_id}\0{name}".encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:4], "big") % (2**32 - 1) + 1


def _default_runner() -> Runner:
    module = importlib.import_module("._gmm_stan", package=__package__)
    runner_type = getattr(module, "CmdStanRunner")
    return runner_type()


def _base_score_rows(rows: pd.DataFrame, study_id: str) -> pd.DataFrame:
    score_values = pd.to_numeric(rows["s"], errors="coerce").to_numpy(dtype=float)
    out = pd.DataFrame(index=pd.RangeIndex(len(rows)))
    out["score_urn"] = study_id
    out["hgvs_nt"] = rows["hgvs_nt"].to_numpy() if "hgvs_nt" in rows.columns else np.nan
    out["hgvs_pro"] = rows["hgvs_pro"].to_numpy()
    out["s"] = score_values
    for identity in _SCORE_IDENTITY_COLUMNS:
        if identity in rows.columns:
            out[identity] = rows[identity].to_numpy()
    return out


def _output_candidate_order(model: str, top_model: str | None) -> tuple[str, ...]:
    candidates = FLEXIBLE_CANDIDATES if model == "flexible" else MODE_CANDIDATES
    if top_model is None:
        return tuple(sorted(candidates))
    return (
        top_model,
        *tuple(candidate for candidate in sorted(candidates) if candidate != top_model),
    )


def _sentinel_frame(columns: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame([{column: np.nan for column in columns}], columns=columns)


def _study_sentinel_frame(columns: Sequence[str], study_id: str) -> pd.DataFrame:
    sentinel = _sentinel_frame(columns)
    if "score_urn" in sentinel.columns:
        sentinel["score_urn"] = study_id
    return sentinel


def _param_value(
    params: Mapping[str, Any] | pd.Series | pd.DataFrame, name: str
) -> Any:
    if isinstance(params, pd.DataFrame):
        if len(params) != 1:
            raise ValueError("params DataFrame must contain exactly one row")
        return params.iloc[0][name]
    return params[name]


def _scalar_or_array(result: np.ndarray, original: Any) -> Any:
    if np.ndim(original) == 0:
        return float(result)
    return result


def _validate_model(model: str) -> None:
    if model not in {"flexible", "modes"}:
        raise ValueError("model must be 'flexible' or 'modes'")
