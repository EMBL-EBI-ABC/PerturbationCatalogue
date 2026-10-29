"""CmdStan execution and compact posterior/PSIS-LOO adaptation for MaveDB GMMs.

The large generated ``log_lik`` matrix is streamed one CSV row at a time to a
temporary binary file and read back in budget-sized memory-mapped blocks. This
module intentionally never calls CmdStanPy's
``draws()``, ``stan_variable()`` or ``draws_pd()`` methods, each of which can
materialize the full likelihood array in memory.

Array conventions are explicit:

* ``summarize_draws`` accepts scalar variables with shape ``(chain, draw)``.
* ``compute_loo`` accepts ``log_lik`` with shape ``(draw, chain, observation)``.
* ``LOOResult.pointwise`` is the one-dimensional ``elpd_loo_i`` vector used
  for paired model comparison; ``pareto_k`` and ``relative_eff`` have the
  same observation order.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import mmap
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.special import logsumexp

LOGGER = logging.getLogger(__name__)
STAN_DIR = Path(__file__).resolve().parent / "stan"
_CANDIDATE_STAN_FILES = {
    "flexible": "fixed_gmm_flexible.stan",
    "ko": "fixed_gmm_ko.stan",
    "del": "fixed_gmm_del.stan",
    "neut": "fixed_gmm_neut.stan",
    "pos": "fixed_gmm_pos.stan",
}
_DEFAULT_CACHE_DIR = Path(__file__).resolve().parents[1] / "tmp" / "gmm_stan_cache"
_ARVIZ_VERSION = "0.22.0"
_SUMMARY_COLUMNS = [
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
]
_SUMMARY_VARIABLE_MAP = {
    "frac_non_neut[1]": "frac_ko",
    "frac_non_neut[2]": "frac_del",
    "frac_non_neut[3]": "frac_pos",
    "y_median": "median",
    "y_mad": "mad",
}
_POSTERIOR_VARIABLE_ORDER = [
    "lp__",
    "frac_neut",
    "frac_non_neut[1]",
    "frac_non_neut[2]",
    "frac_non_neut[3]",
    "theta_noise",
    "theta_ko",
    "theta_del",
    "theta_neut",
    "theta_pos",
    "beta",
    "delta",
    "sigma_neut",
]
_SCALAR_VARIABLES = {
    "lp__",
    "frac_neut",
    "frac_non_neut[1]",
    "frac_non_neut[2]",
    "frac_non_neut[3]",
    "theta_noise",
    "theta_ko",
    "theta_del",
    "theta_neut",
    "theta_pos",
    "beta",
    "delta",
    "sigma_neut",
    "y_median",
    "y_mad",
}
_INCLUDE_RE = re.compile(r'^\s*#include\s+(?:"([^"]+)"|<([^>]+)>|([^\s]+))')
_LOG_LIK_RE = re.compile(r"^log_lik(?:\[(\d+)\]|\.(\d+))$")
_FRAC_NON_NEUT_RE = re.compile(r"^frac_non_neut(?:\[(\d+)\]|\.(\d+))$")


class SamplerError(RuntimeError):
    """Expected candidate-level sampling or numerical failure."""

    def __init__(
        self,
        message: str,
        *,
        candidate: str | None = None,
        warnings: Sequence[str] | None = None,
    ) -> None:
        super().__init__(message)
        self.candidate = candidate
        if isinstance(warnings, str):
            self.warnings = (warnings,)
        else:
            self.warnings = tuple(warnings or ())


@dataclass
class CandidateFit:
    """Compact candidate result, with no CmdStan object or likelihood cube."""

    candidate: str
    posterior: pd.DataFrame
    fit_summary: pd.DataFrame
    params: dict[str, float]
    loo: dict[str, float]
    pointwise: np.ndarray
    pareto_k: np.ndarray
    warnings: list[str]


@dataclass(frozen=True)
class LOOResult:
    """Pointwise PSIS-LOO values and aggregate statistics."""

    summary: dict[str, float]
    pointwise: np.ndarray
    lpd_i: np.ndarray
    p_loo_i: np.ndarray
    pareto_k: np.ndarray
    relative_eff: np.ndarray


@dataclass(frozen=True)
class DiskLogLikelihood:
    """Sequential binary likelihood store ordered as ``(chain, draw, obs)``."""

    path: Path
    chains: int
    draws: int
    observations: int

    def __post_init__(self) -> None:
        expected_bytes = (
            self.chains * self.draws * self.observations * np.dtype(np.float64).itemsize
        )
        if self.path.stat().st_size != expected_bytes:
            raise ValueError(
                f"log_lik store has {self.path.stat().st_size} bytes; "
                f"expected {expected_bytes}"
            )

    @property
    def shape(self) -> tuple[int, int, int]:
        """Expose the public ``(draw, chain, observation)`` dimension order."""

        return self.draws, self.chains, self.observations

    def read_observation_block(self, start: int, stop: int) -> np.ndarray:
        """Copy one bounded observation block and evict its mapped file pages."""

        if not 0 <= start < stop <= self.observations:
            raise ValueError("observation block must be a nonempty in-range slice")
        mapped = np.memmap(
            self.path,
            mode="r",
            dtype=np.float64,
            shape=(self.chains, self.draws, self.observations),
            order="C",
        )
        mapping = mapped._mmap
        try:
            block = np.array(mapped[:, :, start:stop], dtype=np.float64, copy=True)
        finally:
            del mapped
            _evict_and_close_mmap(mapping)
        return block


def _evict_and_close_mmap(mapping: mmap.mmap) -> None:
    """Drop file-backed resident pages, then unmap portably."""

    dont_need = getattr(mmap, "MADV_DONTNEED", None)
    advise = getattr(mapping, "madvise", None)
    if dont_need is not None and advise is not None:
        try:
            advise(dont_need)
        except (OSError, ValueError):
            # Closing the mapping below remains a portable RSS-release fallback.
            pass
    mapping.close()


def _import_arviz() -> tuple[Any, Any]:
    """Import the pinned ArviZ stack after redirecting optional caches."""

    cache_root = Path(tempfile.gettempdir()) / "mavedb-gmm-caches"
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    os.environ.setdefault("NUMBA_CACHE_DIR", str(cache_root / "numba"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root / "xdg"))
    for variable in ("MPLCONFIGDIR", "NUMBA_CACHE_DIR", "XDG_CACHE_HOME"):
        Path(os.environ[variable]).mkdir(parents=True, exist_ok=True)

    import arviz as az
    import xarray as xr

    if az.__version__ != _ARVIZ_VERSION:
        raise RuntimeError(
            f"MaveDB GMM requires arviz=={_ARVIZ_VERSION}; found {az.__version__}"
        )
    return az, xr


def _diagnostic_value(result: Any) -> float:
    """Extract one scalar from ArviZ's one-variable xarray Dataset result."""

    if hasattr(result, "to_array"):
        result = result.to_array().values
    else:
        result = result.values
    return float(np.asarray(result).item())


def _canonical_csv_name(name: str) -> str:
    """Normalize CmdStan's dotted array columns to R-style indexed names."""

    name = name.strip()
    match = _FRAC_NON_NEUT_RE.fullmatch(name)
    if match:
        index = match.group(1) or match.group(2)
        return f"frac_non_neut[{index}]"
    return name


def _summary_variable_name(name: str) -> str:
    canonical = _canonical_csv_name(name)
    return _SUMMARY_VARIABLE_MAP.get(canonical, canonical)


def _arrays_from_draws(
    draws: pd.DataFrame | Mapping[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Return variables as equal-shaped ``(chain, draw)`` float arrays."""

    if isinstance(draws, pd.DataFrame):
        missing = {"chain", "iteration"}.difference(draws.columns)
        if missing:
            raise ValueError(
                f"draw table is missing identifier columns: {sorted(missing)}"
            )
        identifiers = draws[["chain", "iteration"]]
        if identifiers.isna().any().any() or identifiers.duplicated().any():
            raise ValueError(
                "draw table must have unique, nonmissing chain/iteration keys"
            )

        chain_ids = np.sort(draws["chain"].unique())
        iterations = np.sort(draws["iteration"].unique())
        if chain_ids.size == 0 or iterations.size == 0:
            raise ValueError("draw table is empty")
        if not np.array_equal(chain_ids, np.arange(1, chain_ids.size + 1)):
            raise ValueError("chain identifiers must be contiguous and one-based")
        if not np.array_equal(iterations, np.arange(1, iterations.size + 1)):
            raise ValueError("iteration identifiers must be contiguous and one-based")
        expected_rows = int(chain_ids.size * iterations.size)
        if len(draws) != expected_rows:
            raise ValueError(
                "each chain must contain the same complete set of iterations"
            )

        variables: dict[str, np.ndarray] = {}
        for raw_name in draws.columns:
            if raw_name in {"chain", "iteration", "draw", "model"}:
                continue
            name = _canonical_csv_name(str(raw_name))
            if name in variables:
                raise ValueError(
                    f"duplicate draw variable after name normalization: {name}"
                )
            ordered = draws.sort_values(["chain", "iteration"], kind="stable")
            values = pd.to_numeric(ordered[raw_name], errors="raise").to_numpy(
                dtype=np.float64
            )
            variables[name] = values.reshape(chain_ids.size, iterations.size)
        return variables

    variables = {}
    for raw_name, raw_values in draws.items():
        name = _canonical_csv_name(str(raw_name))
        if name in variables:
            raise ValueError(
                f"duplicate draw variable after name normalization: {name}"
            )
        values = np.asarray(raw_values, dtype=np.float64)
        if values.ndim == 1:
            values = values[np.newaxis, :]
        if values.ndim != 2:
            raise ValueError(
                f"draw variable {name!r} must have (chain, draw) dimensions"
            )
        variables[name] = values

    shapes = {value.shape for value in variables.values()}
    if len(shapes) > 1:
        raise ValueError("all draw variables must have identical chain/draw dimensions")
    return variables


def summarize_draws(
    draws: pd.DataFrame | Mapping[str, np.ndarray],
) -> pd.DataFrame:
    """Summarize scalar posterior draws with R-compatible columns.

    ``draws`` may be a wide table with one-based ``chain`` and ``iteration``
    columns, or a mapping of Stan variable names to arrays shaped
    ``(chain, draw)``. Quantiles use NumPy's linear/type-7 convention, ``sd``
    is the sample standard deviation, and ``mad`` uses R's 1.4826 normalizing
    constant. Rank-normalized diagnostics use the pinned ArviZ 0.22.0 API.
    """

    variables = _arrays_from_draws(draws)
    if not variables:
        return pd.DataFrame(columns=_SUMMARY_COLUMNS)

    az, xr = _import_arviz()
    first_shape = next(iter(variables.values())).shape
    chains, iterations = first_shape
    rows: list[dict[str, float | str]] = []

    for name, chain_draws in variables.items():
        if chain_draws.shape != first_shape:
            raise ValueError(
                "all draw variables must have identical chain/draw dimensions"
            )
        values = chain_draws.reshape(-1)
        if not np.isfinite(values).all():
            raise SamplerError(f"posterior variable {name!r} contains non-finite draws")

        constant = bool(np.all(values == values[0]))
        rhat = np.nan
        ess_bulk = np.nan
        ess_tail = np.nan
        if not constant and chains >= 2 and iterations >= 4:
            sample = xr.DataArray(chain_draws, dims=("chain", "draw"))
            rhat = _diagnostic_value(az.rhat(sample, method="rank"))
            ess_bulk = _diagnostic_value(az.ess(sample, method="bulk"))
            ess_tail = _diagnostic_value(az.ess(sample, method="tail"))

        rows.append(
            {
                "variable": _summary_variable_name(name),
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "sd": float(np.std(values, ddof=1)) if values.size > 1 else np.nan,
                "mad": float(np.median(np.abs(values - np.median(values))) * 1.4826),
                "q5": float(np.quantile(values, 0.05, method="linear")),
                "q95": float(np.quantile(values, 0.95, method="linear")),
                "rhat": rhat,
                "ess_bulk": ess_bulk,
                "ess_tail": ess_tail,
            }
        )

    return pd.DataFrame(rows, columns=_SUMMARY_COLUMNS)


def _relative_efficiency(
    az: Any,
    xr: Any,
    log_likelihood: np.ndarray,
) -> float:
    """Compute CmdStanR/loo relative efficiency for one observation."""

    # loo::relative_eff.array uses posterior::ess_mean(exp(log_lik)) / S.
    # Shifting by the maximum leaves ESS unchanged and avoids exponent overflow.
    if np.all(np.ptp(log_likelihood, axis=1) == 0.0):
        return np.nan
    maximum = float(np.max(log_likelihood))
    likelihood = np.exp(log_likelihood - maximum)
    sample = xr.DataArray(likelihood, dims=("chain", "draw"))
    effective_sample_size = _diagnostic_value(az.ess(sample, method="mean"))
    return effective_sample_size / float(likelihood.size)


def _psis_log_weights_fixed_tail(
    log_ratios: np.ndarray,
    relative_efficiency: float,
) -> tuple[np.ndarray, float]:
    """Smooth exactly R loo's fixed-size tail, including cutoff ties.

    ArviZ 0.22's high-level ``psislw`` selects tail samples with the strict
    predicate ``x > cutoff``. The R loo 2.10.1 implementation instead selects
    the last ``ceil(min(0.2*S, 3*sqrt(S/reff)))`` order statistics, which can
    include any observations tied at the cutoff. Use ArviZ's pinned GPD fit
    and inverse-CDF primitives while preserving the R tail membership rule.
    """

    from arviz.stats.stats import _gpdfit, _gpinv

    centered = np.asarray(log_ratios, dtype=np.float64).copy()
    if centered.ndim != 1 or centered.size == 0:
        raise ValueError("log importance ratios must be a nonempty vector")
    centered -= float(np.max(centered))
    sample_count = centered.size
    tail_length = int(
        np.ceil(
            min(
                0.2 * sample_count,
                3.0 * np.sqrt(sample_count / relative_efficiency),
            )
        )
    )
    pareto_k = np.inf

    if tail_length >= 5:
        order = np.argsort(centered, kind="stable")
        tail_indices = order[-tail_length:]
        tail = centered[tail_indices]
        if np.ptp(tail) >= np.finfo(np.float64).eps / 100.0:
            cutoff = centered[order[-tail_length - 1]]
            exp_cutoff = float(np.exp(cutoff))
            excesses = np.exp(tail) - exp_cutoff
            fitted_k, sigma = _gpdfit(excesses)
            if np.isfinite(fitted_k) and np.isfinite(sigma) and sigma > 0:
                probabilities = (
                    np.arange(tail_length, dtype=np.float64) + 0.5
                ) / tail_length
                smoothed_excesses = _gpinv(probabilities, fitted_k, sigma)
                centered[tail_indices] = np.log(smoothed_excesses + exp_cutoff)
                pareto_k = float(fitted_k)

    centered[centered > 0.0] = 0.0
    centered -= float(logsumexp(centered))
    return centered, pareto_k


def _observation_block_size(
    draws: int,
    chains: int,
    observations: int,
    memory_budget_mb: int,
) -> int:
    """Choose a conservative likelihood block bounded by the RAM budget."""

    if memory_budget_mb <= 0:
        raise ValueError("memory_budget_mb must be positive")
    bytes_per_observation = draws * chains * np.dtype(np.float64).itemsize
    budget_bytes = memory_budget_mb * 1024**2
    # Reserve most of the process budget for CmdStan, Python frames and the
    # per-observation PSIS work arrays. The mapped pages and copied block each
    # stay within this smaller allowance.
    block_budget_bytes = min(64 * 1024**2, budget_bytes // 8)
    if budget_bytes < bytes_per_observation * 8:
        raise SamplerError(
            "memory budget is too small for one observation's PSIS workspace"
        )
    return max(1, min(observations, block_budget_bytes // bytes_per_observation))


def compute_loo(
    log_lik: np.ndarray | DiskLogLikelihood,
    memory_budget_mb: int = 2048,
) -> LOOResult:
    """Compute chain-aware PSIS-LOO for ``log_lik[draw, chain, observation]``.

    In-memory arrays are read one observation at a time. ``DiskLogLikelihood``
    stores rows sequentially as ``(chain, draw, observation)``; it is read in
    bounded observation blocks, copied into a workspace sized from
    ``memory_budget_mb``, then unmapped with ``MADV_DONTNEED`` where available.
    No operation creates a boolean mask or a second resident copy of the full
    likelihood cube.
    """

    disk_store = log_lik if isinstance(log_lik, DiskLogLikelihood) else None
    if disk_store is not None:
        draws, chains, observations = disk_store.shape
    else:
        values = np.asanyarray(log_lik)
        if values.ndim != 3:
            raise ValueError("log_lik must have shape (draw, chain, observation)")
        draws, chains, observations = values.shape
    if draws * chains < 4:
        raise SamplerError("PSIS-LOO requires at least four retained draws")
    if observations == 0:
        raise SamplerError("PSIS-LOO requires at least one observation")

    az, xr = _import_arviz()
    lpd_i = np.empty(observations, dtype=np.float64)
    elpd_loo_i = np.empty(observations, dtype=np.float64)
    p_loo_i = np.empty(observations, dtype=np.float64)
    pareto_k = np.empty(observations, dtype=np.float64)
    relative_eff = np.empty(observations, dtype=np.float64)

    if disk_store is None:
        likelihood_blocks = (
            (start, values[:, :, start : start + 1].transpose(1, 0, 2))
            for start in range(observations)
        )
    else:
        block_size = _observation_block_size(
            draws, chains, observations, memory_budget_mb
        )
        likelihood_blocks = _disk_likelihood_blocks(disk_store, block_size)

    for block_start, block in likelihood_blocks:
        if block.shape[0] != chains or block.shape[1] != draws:
            raise ValueError("likelihood block has unexpected chain/draw dimensions")
        block_observations = block.shape[2]
        for local_observation in range(block_observations):
            observation = block_start + local_observation
            chain_draws = np.asarray(block[:, :, local_observation], dtype=np.float64)
            if not np.isfinite(chain_draws).all():
                raise SamplerError(
                    f"generated log_lik contains non-finite values at observation {observation + 1}"
                )
            flat_log_lik = chain_draws.reshape(-1)
            reff = _relative_efficiency(az, xr, chain_draws)
            if not np.isfinite(reff) or reff <= 0:
                raise SamplerError(
                    "could not compute positive finite relative efficiency for "
                    f"observation {observation + 1}"
                )
            relative_eff[observation] = reff

            log_ratios = -flat_log_lik
            smoothed_log_weights, pareto = _psis_log_weights_fixed_tail(
                log_ratios, reff
            )
            lpd_i[observation] = float(
                logsumexp(flat_log_lik) - np.log(flat_log_lik.size)
            )
            elpd_loo_i[observation] = float(
                logsumexp(flat_log_lik + smoothed_log_weights)
            )
            p_loo_i[observation] = lpd_i[observation] - elpd_loo_i[observation]
            pareto_k[observation] = pareto
        if disk_store is not None:
            del block

    def standard_error(pointwise: np.ndarray) -> float:
        if pointwise.size < 2:
            return np.nan
        return float(np.sqrt(pointwise.size * np.var(pointwise, ddof=1)))

    elpd_loo = float(np.sum(elpd_loo_i))
    summary = {
        "elpd_loo": elpd_loo,
        "se_elpd_loo": standard_error(elpd_loo_i),
        "p_loo": float(np.sum(p_loo_i)),
        "se_p_loo": standard_error(p_loo_i),
        "looic": -2.0 * elpd_loo,
        "se_looic": 2.0 * standard_error(elpd_loo_i),
    }
    return LOOResult(
        summary=summary,
        pointwise=elpd_loo_i,
        lpd_i=lpd_i,
        p_loo_i=p_loo_i,
        pareto_k=pareto_k,
        relative_eff=relative_eff,
    )


def _disk_likelihood_blocks(
    store: DiskLogLikelihood,
    block_size: int,
) -> Iterator[tuple[int, np.ndarray]]:
    """Yield copied chain/draw blocks from a disk-backed likelihood store."""

    for start in range(0, store.observations, block_size):
        stop = min(start + block_size, store.observations)
        yield start, store.read_observation_block(start, stop)


def _initial_values(candidate: str) -> dict[str, Any]:
    initial = {
        "flexible": {
            "sigma_neut": 0.1,
            "beta": 1.5,
            "delta": 0.0,
            "frac_neut": 0.8,
            "frac_non_neut": [0.8, 0.19, 0.01],
        },
        "ko": {"sigma_neut": 0.1, "beta": 2.0, "delta": -0.4, "frac_neut": 0.65},
        "del": {"sigma_neut": 0.1, "beta": 1.25, "delta": -0.1, "frac_neut": 0.7},
        "neut": {"beta": 1.01, "delta": 0.0, "frac_neut": 0.9},
        "pos": {"sigma_neut": 0.1, "beta": 1.25, "delta": 0.0, "frac_neut": 0.75},
    }
    return initial[candidate]


def _included_stan_files(stan_file: Path) -> list[Path]:
    """Collect a model and its transitive ``#include`` sources."""

    pending = [stan_file.resolve()]
    included: set[Path] = set()
    while pending:
        current = pending.pop()
        if current in included:
            continue
        included.add(current)
        text = current.read_text(encoding="utf-8")
        for line in text.splitlines():
            match = _INCLUDE_RE.match(line)
            if not match:
                continue
            include_name = match.group(1) or match.group(2) or match.group(3)
            include_path = (current.parent / include_name).resolve()
            if not include_path.is_file():
                include_path = (STAN_DIR / include_name).resolve()
            if not include_path.is_file():
                raise FileNotFoundError(
                    f"Stan include {include_name!r} referenced from {current} was not found"
                )
            pending.append(include_path)
    return sorted(included)


def _compiler_identity() -> dict[str, str]:
    configured_cxx = os.environ.get("CXX", "c++")
    command = shlex.split(configured_cxx)
    compiler = shutil.which(command[0]) if command else None
    version = "unavailable"
    if compiler:
        result = subprocess.run(
            [compiler, "--version"],
            capture_output=True,
            check=False,
            text=True,
            timeout=10,
        )
        version = (
            (result.stdout or result.stderr).splitlines()[0]
            if (result.stdout or result.stderr)
            else "unknown"
        )
    return {
        "configured_cxx": configured_cxx,
        "resolved_cxx": compiler or "unavailable",
        "version": version,
        "cxxflags": os.environ.get("CXXFLAGS", ""),
    }


def model_cache_key(candidate: str) -> str:
    """Return the content/toolchain digest used for this candidate executable."""

    if candidate not in _CANDIDATE_STAN_FILES:
        raise ValueError(f"unknown GMM candidate: {candidate!r}")
    import cmdstanpy

    stan_file = STAN_DIR / _CANDIDATE_STAN_FILES[candidate]
    options = {"cpp_options": {}, "stanc_options": {}}
    cmdstan_path = str(cmdstanpy.cmdstan_path())
    cmdstan_version = str(cmdstanpy.cmdstan_version())
    local_make = Path(cmdstan_path) / "make" / "local"
    fingerprint: dict[str, Any] = {
        "candidate": candidate,
        "options": options,
        "cmdstanpy": cmdstanpy.__version__,
        "cmdstan_path": cmdstan_path,
        "cmdstan_version": cmdstan_version,
        "python": os.sys.version,
        "compiler": _compiler_identity(),
        "make_local": (
            local_make.read_text(encoding="utf-8") if local_make.exists() else ""
        ),
    }
    digest = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode("utf-8"))
    for source in _included_stan_files(stan_file):
        relative = source.relative_to(STAN_DIR.resolve()).as_posix()
        digest.update(relative.encode("utf-8"))
        digest.update(source.read_bytes())
    return digest.hexdigest()


@contextmanager
def _exclusive_file_lock(lock_path: Path) -> Iterator[None]:
    """Serialize compilation across processes on Unix and Windows."""

    lock_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import fcntl
    except ImportError:
        fcntl = None

    with lock_path.open("a+b") as handle:
        if fcntl is not None:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        else:
            import msvcrt

            handle.seek(0, os.SEEK_END)
            if handle.tell() == 0:
                handle.write(b"0")
                handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)


def _diagnostic_warnings(stdout_files: Sequence[str]) -> list[str]:
    warnings: list[str] = []
    for stdout_file in stdout_files:
        path = Path(stdout_file)
        if not path.is_file():
            continue
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if "Warning:" in line or "WARNING:" in line:
                    message = line.strip()
                    if message and message not in warnings:
                        warnings.append(message)
    return warnings


def _read_cmdstan_csvs(
    csv_files: Sequence[str],
    candidate: str,
    chains: int,
    iterations: int,
    observations: int,
    scratch: Path,
    max_treedepth: int,
) -> tuple[pd.DataFrame, dict[str, np.ndarray], DiskLogLikelihood, list[str]]:
    """Stream retained draw rows into compact scalar arrays and a disk file."""

    if len(csv_files) != chains:
        raise SamplerError(
            f"CmdStan returned {len(csv_files)} chain CSVs; expected {chains}"
        )

    likelihood_path = scratch / "log_lik.float64.bin"
    scalar_draws: dict[str, np.ndarray] = {}
    expected_scalar_names: set[str] | None = None
    warning_messages: list[str] = []
    divergent_count = 0
    treedepth_count = 0

    with likelihood_path.open("wb") as likelihood_handle:
        for chain_index, csv_path in enumerate(csv_files):
            with Path(csv_path).open(newline="", encoding="utf-8") as handle:
                reader = csv.reader(handle)
                header: list[str] | None = None
                for row in reader:
                    if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                        continue
                    header = [column.strip() for column in row]
                    break
                if header is None:
                    raise SamplerError(f"CmdStan CSV has no header: {csv_path}")

                scalar_columns: dict[str, int] = {}
                log_lik_columns: dict[int, int] = {}
                diagnostic_columns: dict[str, int] = {}
                for column_index, raw_name in enumerate(header):
                    name = _canonical_csv_name(raw_name)
                    if name in _SCALAR_VARIABLES:
                        if name in scalar_columns:
                            raise SamplerError(f"duplicate scalar CSV column {name!r}")
                        scalar_columns[name] = column_index
                    likelihood_match = _LOG_LIK_RE.fullmatch(raw_name)
                    if likelihood_match:
                        observation = int(
                            likelihood_match.group(1) or likelihood_match.group(2)
                        )
                        log_lik_columns[observation] = column_index
                    if raw_name in {"divergent__", "treedepth__"}:
                        diagnostic_columns[raw_name] = column_index

                if sorted(log_lik_columns) != list(range(1, observations + 1)):
                    raise SamplerError(
                        f"CmdStan CSV contains {len(log_lik_columns)} log_lik columns; "
                        f"expected observations 1..{observations}"
                    )
                likelihood_order = [
                    log_lik_columns[observation]
                    for observation in range(1, observations + 1)
                ]
                required = {
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
                    "y_median",
                    "y_mad",
                }
                missing = required.difference(scalar_columns)
                if missing:
                    raise SamplerError(
                        f"CmdStan CSV is missing scalar columns: {sorted(missing)}"
                    )
                if candidate == "flexible":
                    expected_fractions = {
                        "frac_non_neut[1]",
                        "frac_non_neut[2]",
                        "frac_non_neut[3]",
                    }
                    missing_fractions = expected_fractions.difference(scalar_columns)
                    if missing_fractions:
                        raise SamplerError(
                            "flexible CmdStan CSV is missing simplex columns: "
                            f"{sorted(missing_fractions)}"
                        )

                current_scalar_names = set(scalar_columns)
                if expected_scalar_names is None:
                    expected_scalar_names = current_scalar_names
                elif current_scalar_names != expected_scalar_names:
                    raise SamplerError(
                        "CmdStan chain CSVs have different scalar headers"
                    )

                for name in scalar_columns:
                    scalar_draws.setdefault(
                        name, np.empty((chains, iterations), dtype=np.float64)
                    )

                row_count = 0
                for row in reader:
                    if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                        continue
                    if row_count >= iterations:
                        raise SamplerError(
                            f"CmdStan CSV has excess retained draws: {csv_path}"
                        )
                    if len(row) != len(header):
                        raise SamplerError(
                            f"CmdStan CSV row {row_count + 1} has {len(row)} fields; "
                            f"expected {len(header)}: {csv_path}"
                        )
                    try:
                        for name, column_index in scalar_columns.items():
                            scalar_draws[name][chain_index, row_count] = float(
                                row[column_index]
                            )
                        likelihood_row = np.fromiter(
                            (
                                float(row[column_index])
                                for column_index in likelihood_order
                            ),
                            dtype=np.float64,
                            count=observations,
                        )
                        if not np.isfinite(likelihood_row).all():
                            raise SamplerError(
                                f"CmdStan generated non-finite log_lik at row {row_count + 1}"
                            )
                        likelihood_handle.write(likelihood_row.tobytes(order="C"))
                        if "divergent__" in diagnostic_columns:
                            divergent_count += int(
                                float(row[diagnostic_columns["divergent__"]]) > 0
                            )
                        if "treedepth__" in diagnostic_columns:
                            treedepth_count += int(
                                float(row[diagnostic_columns["treedepth__"]])
                                >= max_treedepth
                            )
                    except (ValueError, OverflowError) as error:
                        raise SamplerError(
                            f"CmdStan CSV contains a non-numeric draw at row {row_count + 1}: {csv_path}"
                        ) from error
                    row_count += 1
                if row_count != iterations:
                    raise SamplerError(
                        f"CmdStan CSV has {row_count} retained draws; expected {iterations}: {csv_path}"
                    )

    log_lik = DiskLogLikelihood(
        path=likelihood_path,
        chains=chains,
        draws=iterations,
        observations=observations,
    )
    for name, values in scalar_draws.items():
        if not np.isfinite(values).all():
            raise SamplerError(f"CmdStan generated non-finite draws for {name!r}")

    if divergent_count:
        warning_messages.append(
            f"CmdStan reported {divergent_count} divergent transitions"
        )
    if treedepth_count:
        warning_messages.append(
            f"CmdStan reported {treedepth_count} transitions at max_treedepth={max_treedepth}"
        )

    posterior_columns = [
        name for name in _POSTERIOR_VARIABLE_ORDER if name in scalar_draws
    ]
    if candidate == "neut":
        posterior_columns = [name for name in posterior_columns if name != "frac_neut"]
    exported_names = {
        "frac_non_neut[1]": "frac_ko",
        "frac_non_neut[2]": "frac_del",
        "frac_non_neut[3]": "frac_pos",
    }
    records: dict[str, Any] = {
        "iteration": np.tile(np.arange(1, iterations + 1, dtype=np.int64), chains),
        "chain": np.repeat(np.arange(1, chains + 1, dtype=np.int64), iterations),
        "draw": np.arange(1, chains * iterations + 1, dtype=np.int64),
    }
    for name in posterior_columns:
        target_name = exported_names.get(name, name)
        records[target_name] = scalar_draws[name].reshape(-1)
    posterior = pd.DataFrame(records)

    # Keep all captured scalar variables for summaries, including the neutral
    # fraction omitted only from the neutral candidate's exported posterior.
    summary_draws = {
        name: values
        for name, values in scalar_draws.items()
        if name in _SCALAR_VARIABLES
    }
    return posterior, summary_draws, log_lik, warning_messages


class CmdStanRunner:
    """Compile/cache and sample one candidate without materializing log_lik."""

    def __init__(self, cache_dir: str | os.PathLike[str] | None = None) -> None:
        configured_cache = cache_dir or os.environ.get("MAVEDB_GMM_CACHE_DIR")
        self.cache_dir = (
            Path(configured_cache) if configured_cache else _DEFAULT_CACHE_DIR
        )

    def get_cmdstan_model(self, candidate: str) -> Any:
        """Return the digest-keyed compiled model, serializing first compilation."""

        if candidate not in _CANDIDATE_STAN_FILES:
            raise ValueError(f"unknown GMM candidate: {candidate!r}")
        import cmdstanpy

        stan_file = STAN_DIR / _CANDIDATE_STAN_FILES[candidate]
        if not stan_file.is_file():
            raise FileNotFoundError(f"GMM Stan source was not found: {stan_file}")
        digest = model_cache_key(candidate)
        model_cache = self.cache_dir / f"{candidate}-{digest}"
        model_cache.mkdir(parents=True, exist_ok=True)
        lock_path = self.cache_dir / f"{candidate}-{digest}.lock"

        with _exclusive_file_lock(lock_path):
            for source in _included_stan_files(stan_file):
                source_copy = model_cache / source.relative_to(STAN_DIR.resolve())
                source_copy.parent.mkdir(parents=True, exist_ok=True)
                if (
                    not source_copy.is_file()
                    or source_copy.read_bytes() != source.read_bytes()
                ):
                    shutil.copy2(source, source_copy)
            cached_stan_file = model_cache / stan_file.name
            executable = model_cache / stan_file.stem
            compiled = executable.is_file() and os.access(executable, os.X_OK)
            if compiled:
                return cmdstanpy.CmdStanModel(
                    stan_file=str(cached_stan_file),
                    exe_file=str(executable),
                    stanc_options={},
                    cpp_options={},
                )
            return cmdstanpy.CmdStanModel(
                stan_file=str(cached_stan_file),
                stanc_options={},
                cpp_options={},
            )

    def __call__(
        self,
        scores: np.ndarray,
        candidate: str,
        settings: Any,
        seed: int | None = None,
    ) -> CandidateFit:
        if candidate not in _CANDIDATE_STAN_FILES:
            raise ValueError(f"unknown GMM candidate: {candidate!r}")
        values = np.asarray(scores, dtype=np.float64).reshape(-1)
        if values.size == 0:
            raise SamplerError("cannot fit an empty score vector")
        if not np.isfinite(values).all():
            raise SamplerError("score vector contains non-finite values")
        median = float(np.median(values))
        mad = float(np.median(np.abs(values - median)) * 1.4826)
        if not np.isfinite(mad) or mad <= 0:
            raise SamplerError("sampled score MAD is zero or non-finite")

        chains = int(settings.chains)
        parallel_chains = int(settings.parallel_chains)
        iter_warmup = int(settings.iter_warmup)
        iter_sampling = int(settings.iter_sampling)
        if chains < 1 or parallel_chains < 1 or parallel_chains > chains:
            raise ValueError("parallel_chains must be between one and chains")
        if iter_warmup < 0 or iter_sampling < 1:
            raise ValueError(
                "iter_warmup must be nonnegative and iter_sampling positive"
            )
        refresh = int(settings.refresh)
        if refresh < 0:
            raise ValueError("refresh must be nonnegative")
        if int(settings.memory_budget_mb) <= 0:
            raise ValueError("memory_budget_mb must be positive")

        # Compile before making a per-fit scratch directory: dependency and
        # global compilation failures should fail the batch, not look like a
        # recoverable candidate failure.
        model = self.get_cmdstan_model(candidate)
        scratch = Path(tempfile.mkdtemp(prefix=f"mavedb-gmm-{candidate}-"))
        selected_seed = seed if seed is not None else getattr(settings, "seed", None)
        sampler_warnings: list[str] = []
        log_lik: DiskLogLikelihood | None = None
        succeeded = False

        try:
            try:
                fit = model.sample(
                    data={"N": int(values.size), "y": values},
                    seed=selected_seed,
                    chains=chains,
                    parallel_chains=parallel_chains,
                    iter_warmup=iter_warmup,
                    iter_sampling=iter_sampling,
                    save_warmup=False,
                    adapt_delta=float(settings.adapt_delta),
                    max_treedepth=int(settings.max_treedepth),
                    refresh=max(refresh, 1),
                    inits=_initial_values(candidate),
                    output_dir=str(scratch),
                    sig_figs=18,
                    show_progress=refresh > 0,
                    show_console=False,
                )
            except RuntimeError as error:
                raise SamplerError(
                    f"CmdStan sampling failed for {candidate}; retained files: {scratch}: {error}",
                    candidate=candidate,
                ) from error

            runset = getattr(fit, "runset", None)
            csv_files = getattr(runset, "csv_files", None)
            if not csv_files:
                raise SamplerError(
                    f"CmdStan returned no CSV files for {candidate}; retained files: {scratch}",
                    candidate=candidate,
                )
            stdout_warnings = _diagnostic_warnings(getattr(runset, "stdout_files", []))
            try:
                posterior, scalar_draws, log_lik, parse_warnings = _read_cmdstan_csvs(
                    csv_files=csv_files,
                    candidate=candidate,
                    chains=chains,
                    iterations=iter_sampling,
                    observations=values.size,
                    scratch=scratch,
                    max_treedepth=int(settings.max_treedepth),
                )
            except SamplerError as error:
                raise SamplerError(
                    str(error),
                    candidate=candidate,
                    warnings=(*stdout_warnings, *error.warnings),
                ) from error
            sampler_warnings.extend(parse_warnings)
            sampler_warnings.extend(stdout_warnings)

            try:
                fit_summary = summarize_draws(scalar_draws)
                params = {
                    str(row.variable): float(row.median)
                    for row in fit_summary.itertuples(index=False)
                }
                loo_result = compute_loo(
                    log_lik, memory_budget_mb=int(settings.memory_budget_mb)
                )
            except (
                SamplerError,
                FloatingPointError,
                OverflowError,
                np.linalg.LinAlgError,
            ) as error:
                previous_warnings = getattr(error, "warnings", ())
                raise SamplerError(
                    str(error),
                    candidate=candidate,
                    warnings=(*sampler_warnings, *previous_warnings),
                ) from error
            sampler_warnings.extend(_pareto_warnings(loo_result.pareto_k))
            candidate_fit = CandidateFit(
                candidate=candidate,
                posterior=posterior,
                fit_summary=fit_summary,
                params=params,
                loo=loo_result.summary,
                pointwise=loo_result.pointwise,
                pareto_k=loo_result.pareto_k,
                warnings=sampler_warnings,
            )
            succeeded = True
            return candidate_fit
        finally:
            if log_lik is not None:
                del log_lik
            if succeeded:
                shutil.rmtree(scratch)
            else:
                LOGGER.warning(
                    "Retaining failed CmdStan scratch directory: %s", scratch
                )


def _pareto_warnings(pareto_k: np.ndarray) -> list[str]:
    if pareto_k.size == 0:
        return []
    count = int(np.count_nonzero(pareto_k > 0.7))
    if count == 0:
        return []
    return [
        f"Pareto-k exceeds 0.7 for {count}/{pareto_k.size} observations "
        f"(maximum {float(np.max(pareto_k)):.3g})"
    ]


def compile_model(
    candidate: str,
    cache_dir: str | os.PathLike[str] | None = None,
) -> Any:
    """Public helper for serial precompilation through the same digest cache."""

    return CmdStanRunner(cache_dir=cache_dir).get_cmdstan_model(candidate)
