"""Matplotlib equivalents of the R GMM diagnostic plots."""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from typing import Any

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpecFromSubplotSpec, SubplotSpec
from matplotlib.lines import Line2D

MODE_ORDER = ("ko", "del", "neut", "pos", "noise")
MODE_COLORS = {
    "ko": "#e41a1c",
    "ko/del": "#fb9a99",
    "del": "#ff7f00",
    "del/neut": "#a6cee3",
    "neut": "#377eb8",
    "neut/pos": "#cab2d6",
    "pos": "#984ea3",
    "incoherent": "green",
    "noise": "grey",
}
CONTROL_COLORS = {
    "WT": "#33a02c",
    "Synonymous": "#b2df8a",
    "Nonsense": "#f781bf",
}
UNKNOWN_COLOR = "#9e9e9e"
W_LIMITS = (-0.01, 2.01)
POSTERIOR_EXCLUDED = {
    "score_urn",
    "iteration",
    "chain",
    "draw",
    "batch",
    "model",
    "top_model",
    "theta_noise",
}
PARAMS_EXCLUDED = POSTERIOR_EXCLUDED | {"median", "mad"}
PARAMETER_ORDER = (
    "frac_neut",
    "frac_ko",
    "frac_del",
    "frac_pos",
    "theta_neut",
    "theta_ko",
    "theta_del",
    "theta_pos",
    "sigma_neut",
    "beta",
    "delta",
    "lp__",
)
PARAMETER_LABELS = {
    "frac_neut": "Fraction Neutral",
    "frac_ko": "Subfraction KO",
    "frac_del": "Subfraction Del",
    "frac_pos": "Subfraction Pos",
    "theta_neut": r"$\theta_{\mathrm{neut}}$",
    "theta_ko": r"$\theta_{\mathrm{ko}}$",
    "theta_del": r"$\theta_{\mathrm{del}}$",
    "theta_pos": r"$\theta_{\mathrm{pos}}$",
    "sigma_neut": r"$\sigma_{\mathrm{neut}}$",
    "beta": r"$\beta$",
    "delta": r"$\alpha$",
    "lp__": "Log-Likelihood",
}
RHAT_CMAP = LinearSegmentedColormap.from_list("gmm_rhat", ["darkgrey", "red"])
RHAT_NORM = Normalize(vmin=0.99, vmax=3.0, clip=True)
RHAT_MISSING_COLOR = "#c994c7"
LOO_CMAP = mpl.colormaps["Reds"]


def _frame(data: pd.DataFrame | None) -> pd.DataFrame:
    return data if isinstance(data, pd.DataFrame) else pd.DataFrame()


def _numeric(data: pd.DataFrame, column: str) -> np.ndarray:
    if column not in data:
        return np.array([], dtype=float)
    return pd.to_numeric(data[column], errors="coerce").to_numpy(
        dtype=float, na_value=np.nan
    )


def _finite_column(data: pd.DataFrame, column: str) -> np.ndarray:
    values = _numeric(data, column)
    return values[np.isfinite(values)]


def _mode_labels(data: pd.DataFrame, indices: np.ndarray | None = None) -> np.ndarray:
    if indices is None:
        indices = np.arange(len(data))
    if "mode_density" not in data:
        return np.full(len(indices), "unclassified", dtype=object)
    labels: list[str] = []
    column = data["mode_density"]
    for index in indices:
        value = column.iloc[int(index)]
        if value is None or pd.isna(value):
            labels.append("unclassified")
        else:
            labels.append(str(value))
    return np.asarray(labels, dtype=object)


def _ordered_modes(values: Sequence[str]) -> list[str]:
    present = set(values)
    ordered = [mode for mode in MODE_ORDER if mode in present]
    ordered.extend(
        mode for mode in MODE_COLORS if mode not in MODE_ORDER and mode in present
    )
    ordered.extend(sorted(present - set(ordered)))
    return ordered


def _mode_color(mode: str) -> str:
    return MODE_COLORS.get(mode, UNKNOWN_COLOR)


def _score_limits(scores: pd.DataFrame) -> tuple[float, float]:
    values = _finite_column(scores, "s")
    if values.size == 0:
        return (-1.0, 1.0)
    low = float(np.min(values))
    high = float(np.max(values))
    width = high - low
    padding = max(abs(low) * 0.05, 0.05) if width == 0 else width * 0.05
    return (low - padding, high + padding)


def _histogram_edges(variable: str, limits: tuple[float, float]) -> np.ndarray:
    if variable == "w":
        return np.linspace(-0.005, 2.005, 100)
    return np.linspace(limits[0], limits[1], 100)


def _scaled_histogram(
    ax: Axes,
    scores: pd.DataFrame,
    variable: str,
    edges: np.ndarray,
    orientation: str,
) -> float:
    """Draw stacked counts scaled by one panel-wide maximum, as in R scaled()."""
    values = _numeric(scores, variable)
    finite = np.isfinite(values)
    valid_indices = np.flatnonzero(finite)
    if valid_indices.size == 0:
        return 0.0
    labels = _mode_labels(scores, valid_indices)
    modes = _ordered_modes(labels.tolist())
    counts_by_mode = {
        mode: np.histogram(values[valid_indices[labels == mode]], bins=edges)[0].astype(
            float
        )
        for mode in modes
    }
    global_max = max(
        (float(counts.max()) for counts in counts_by_mode.values()), default=0.0
    )
    if global_max <= 0:
        return 0.0

    scaled_counts = {
        mode: counts / global_max for mode, counts in counts_by_mode.items()
    }
    bottom = np.zeros(len(edges) - 1, dtype=float)
    widths = np.diff(edges)
    for mode in modes:
        heights = scaled_counts[mode]
        for bin_index in np.flatnonzero(heights > 0):
            if orientation == "vertical":
                ax.bar(
                    edges[bin_index],
                    heights[bin_index],
                    width=widths[bin_index],
                    bottom=bottom[bin_index],
                    align="edge",
                    color=_mode_color(mode),
                    alpha=0.5,
                    linewidth=0,
                    zorder=1,
                )
            else:
                ax.barh(
                    edges[bin_index],
                    heights[bin_index],
                    height=widths[bin_index],
                    left=bottom[bin_index],
                    align="edge",
                    color=_mode_color(mode),
                    alpha=0.5,
                    linewidth=0,
                    zorder=1,
                )
        bottom += heights
    return float(bottom.max(initial=0.0))


def _draw_component_density(
    ax: Axes,
    density: pd.DataFrame,
    name: str,
    orientation: str,
    x_limits: tuple[float, float] | None = None,
) -> None:
    if not {"name", "mode", "x", "density"}.issubset(density.columns):
        return
    rows = density.loc[density["name"].astype(str) == name]
    x_values = pd.to_numeric(rows["x"], errors="coerce").to_numpy(
        dtype=float, na_value=np.nan
    )
    y_values = pd.to_numeric(rows["density"], errors="coerce").to_numpy(
        dtype=float, na_value=np.nan
    )
    finite = np.isfinite(x_values) & np.isfinite(y_values)
    if x_limits is not None:
        finite &= (x_values >= x_limits[0]) & (x_values <= x_limits[1])
    if not finite.any():
        return
    rows = rows.iloc[np.flatnonzero(finite)]
    x_values = x_values[finite]
    y_values = y_values[finite]
    global_max = float(np.max(y_values))
    if not np.isfinite(global_max) or global_max <= 0:
        return
    y_values = y_values / global_max
    modes = rows["mode"].astype(str).to_numpy()
    for mode in _ordered_modes(modes.tolist()):
        selected = modes == mode
        order = np.argsort(x_values[selected], kind="stable")
        x = x_values[selected][order]
        y = y_values[selected][order]
        if orientation == "vertical":
            line = ax.plot(x, y, color=_mode_color(mode), linewidth=1.1, zorder=2)[0]
        else:
            line = ax.plot(y, x, color=_mode_color(mode), linewidth=1.1, zorder=2)[0]
        line.set_gid(f"density:{name}:{mode}")


def _control_type(hgvs_pro: Any) -> str | None:
    if not isinstance(hgvs_pro, str):
        return None
    if hgvs_pro in {"p.=", "_wt"}:
        return "WT"
    if hgvs_pro in {"p.(=)", "_sy"}:
        return "Synonymous"
    if re.search(r"^p\.([A-Za-z]{3})\d*\1", hgvs_pro, flags=re.IGNORECASE):
        return "Synonymous"
    # HGVS uses either Ter or * for a stop; accept positions with any digit count.
    if re.search(r"^p\.[A-Za-z]{3}\d+(?:Ter|\*)", hgvs_pro, flags=re.IGNORECASE):
        return "Nonsense"
    return None


def _draw_control_rugs(ax_w: Axes, ax_s: Axes, scores: pd.DataFrame) -> None:
    if "hgvs_pro" not in scores:
        return
    w_values = _numeric(scores, "w")
    s_values = _numeric(scores, "s")
    for index, hgvs_pro in enumerate(scores["hgvs_pro"].tolist()):
        control = _control_type(hgvs_pro)
        if control is None:
            continue
        color = CONTROL_COLORS[control]
        if index < w_values.size and np.isfinite(w_values[index]):
            line = ax_w.plot(
                [w_values[index], w_values[index]],
                [-0.05, 0],
                color=color,
                linewidth=1.0,
                clip_on=False,
                zorder=4,
            )[0]
            line.set_gid(f"control:{control}")
        if index < s_values.size and np.isfinite(s_values[index]):
            line = ax_s.plot(
                [-0.05, 0],
                [s_values[index], s_values[index]],
                color=color,
                linewidth=1.0,
                clip_on=False,
                zorder=4,
            )[0]
            line.set_gid(f"control:{control}")


def _draw_mode_legend(ax: Axes) -> None:
    handles = [
        Line2D([0], [0], color=_mode_color(mode), linewidth=2, label=mode)
        for mode in MODE_ORDER
    ]
    handles.extend(
        Line2D([0], [0], color=CONTROL_COLORS[label], linewidth=2, label=label)
        for label in ("Nonsense", "Synonymous", "WT")
    )
    ax.axis("off")
    ax.legend(
        handles=handles,
        loc="center left",
        ncol=2,
        frameon=False,
        fontsize=7,
        handlelength=1.2,
        columnspacing=0.65,
        labelspacing=0.5,
    )


def draw_transform(
    fig: Figure,
    subspec: SubplotSpec,
    scores: pd.DataFrame,
    density: pd.DataFrame,
) -> None:
    """Draw transformed scores with w/s marginals and control rugs."""
    scores = _frame(scores)
    density = _frame(density)
    score_limits = _score_limits(scores)
    grid = GridSpecFromSubplotSpec(
        2,
        2,
        subplot_spec=subspec,
        width_ratios=(4, 1),
        height_ratios=(1, 4),
        wspace=0.04,
        hspace=0.04,
    )
    ax_w = fig.add_subplot(grid[0, 0])
    ax_legend = fig.add_subplot(grid[0, 1])
    ax_scatter = fig.add_subplot(grid[1, 0])
    ax_s = fig.add_subplot(grid[1, 1])

    w_edges = _histogram_edges("w", score_limits)
    s_edges = _histogram_edges("s", score_limits)
    w_stack_max = _scaled_histogram(ax_w, scores, "w", w_edges, "vertical")
    s_stack_max = _scaled_histogram(ax_s, scores, "s", s_edges, "horizontal")
    _draw_component_density(ax_w, density, "w", "vertical")
    _draw_component_density(ax_s, density, "s", "horizontal", score_limits)

    w_values = _numeric(scores, "w")
    s_values = _numeric(scores, "s")
    finite = np.isfinite(w_values) & np.isfinite(s_values)
    if finite.any():
        labels = _mode_labels(scores, np.flatnonzero(finite))
        w_finite = w_values[finite]
        s_finite = s_values[finite]
        for mode in _ordered_modes(labels.tolist()):
            selected = labels == mode
            points = ax_scatter.scatter(
                w_finite[selected],
                s_finite[selected],
                s=10,
                color=_mode_color(mode),
                marker="o",
                linewidths=0,
                rasterized=True,
                zorder=3,
            )
            points.set_gid(f"scores:{mode}")

    _draw_control_rugs(ax_w, ax_s, scores)
    _draw_mode_legend(ax_legend)

    ax_w.set_xlim(*W_LIMITS)
    ax_w.set_ylim(-0.075, max(1.0, w_stack_max) * 1.06)
    ax_w.set_ylabel("Frequency")
    ax_w.tick_params(axis="x", labelbottom=False)
    ax_w.grid(axis="x", color="grey", linestyle=":", linewidth=0.6, zorder=0)
    ax_w.tick_params(axis="x", bottom=False)
    ax_w.tick_params(axis="y", labelsize=8)

    ax_scatter.set_xlim(*W_LIMITS)
    ax_scatter.set_ylim(*score_limits)
    ax_scatter.set_xlabel("Latent Fitness (w)")
    ax_scatter.set_ylabel("Experiment score (s)")
    ax_scatter.grid(axis="y", color="grey", linestyle=":", linewidth=0.6, zorder=0)
    ax_scatter.tick_params(labelsize=8)

    ax_s.set_ylim(*score_limits)
    ax_s.set_xlim(-0.075, max(1.0, s_stack_max) * 1.06)
    ax_s.set_xlabel("Frequency")
    ax_s.tick_params(axis="y", labelleft=False)
    ax_s.grid(axis="y", color="grey", linestyle=":", linewidth=0.6, zorder=0)
    ax_s.tick_params(axis="y", left=False, right=False)
    ax_s.tick_params(axis="x", labelsize=8)


def _is_true(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return value == 1
    return isinstance(value, str) and value.strip().lower() in {"true", "t", "1"}


def _selected_parameter_rows(params: pd.DataFrame) -> pd.DataFrame:
    if "top_model" not in params:
        return params
    return params.loc[params["top_model"].map(_is_true)]


def _parameter_estimates(params: pd.DataFrame) -> dict[str, list[float]]:
    params = _selected_parameter_rows(params)
    estimates: dict[str, list[float]] = {}
    for variable in params.columns:
        if variable in PARAMS_EXCLUDED:
            continue
        values = pd.to_numeric(params[variable], errors="coerce").to_numpy(
            dtype=float, na_value=np.nan
        )
        finite = values[np.isfinite(values)]
        if finite.size:
            estimates[variable] = finite.tolist()
    return estimates


def _rhat_by_variable(fit_summary: pd.DataFrame) -> dict[str, float]:
    if not {"variable", "rhat"}.issubset(fit_summary.columns):
        return {}
    rows = fit_summary
    if "top_model" in rows:
        top = rows["top_model"].map(_is_true)
        rows = rows.loc[top]
    result: dict[str, float] = {}
    for variable, group in rows.groupby("variable", sort=False, dropna=True):
        values = pd.to_numeric(group["rhat"], errors="coerce").to_numpy(
            dtype=float, na_value=np.nan
        )
        values = values[np.isfinite(values)]
        if values.size:
            result[str(variable)] = float(np.median(values))
    return result


def _rhat_color(value: float | None) -> str:
    if value is None or not np.isfinite(value) or value < 0.99 or value > 3.0:
        return RHAT_MISSING_COLOR
    return mpl.colors.to_hex(RHAT_CMAP(RHAT_NORM(value)))


def _parameter_variables(posterior: pd.DataFrame) -> list[str]:
    numeric = [
        column
        for column in posterior.columns
        if column not in POSTERIOR_EXCLUDED
        and np.isfinite(
            pd.to_numeric(posterior[column], errors="coerce").to_numpy(
                dtype=float, na_value=np.nan
            )
        ).any()
    ]
    ordered = [variable for variable in PARAMETER_ORDER if variable in numeric]
    ordered.extend(variable for variable in numeric if variable not in ordered)
    return ordered


def _histogram_range(values: np.ndarray) -> tuple[float, float]:
    low = float(np.min(values))
    high = float(np.max(values))
    width = high - low
    padding = max(width * 0.01, abs(low) * 0.01, 0.5 if width == 0 else 0.0)
    return low - padding, high + padding


def draw_posterior(
    fig: Figure,
    subspec: SubplotSpec,
    posterior: pd.DataFrame,
    fit_summary: pd.DataFrame,
    params: pd.DataFrame,
) -> None:
    """Draw posterior parameter histograms with Rhat colors and estimates."""
    posterior = _frame(posterior)
    fit_summary = _frame(fit_summary)
    params = _frame(params)
    variables = _parameter_variables(posterior)
    if not variables:
        ax = fig.add_subplot(subspec)
        ax.axis("off")
        ax.text(0.5, 0.5, "No posterior diagnostics", ha="center", va="center")
        return

    ncols = min(4, len(variables))
    nrows = math.ceil(len(variables) / ncols)
    grid = GridSpecFromSubplotSpec(
        nrows + 1,
        ncols,
        subplot_spec=subspec,
        height_ratios=[1.0] * nrows + [0.12],
        hspace=0.65,
        wspace=0.35,
    )
    rhat = _rhat_by_variable(fit_summary)
    estimates = _parameter_estimates(params)
    for index, variable in enumerate(variables):
        ax = fig.add_subplot(grid[index // ncols, index % ncols])
        values = pd.to_numeric(posterior[variable], errors="coerce").to_numpy(
            dtype=float, na_value=np.nan
        )
        values = values[np.isfinite(values)]
        if not values.size:
            ax.axis("off")
            continue
        low, high = _histogram_range(values)
        ax.hist(
            values,
            bins=np.linspace(low, high, 101),
            color=_rhat_color(rhat.get(variable)),
            edgecolor="none",
            linewidth=0,
        )
        ax.set_xlim(low, high)
        ax.set_xlabel(PARAMETER_LABELS.get(variable, variable.replace("_", " ")))
        if index % ncols == 0:
            ax.set_ylabel("Count")
        else:
            ax.set_ylabel("")
        ax.tick_params(labelsize=7)
        ax.xaxis.label.set_size(8)
        setattr(ax, "mave_variable", variable)
        setattr(ax, "mave_rhat", rhat.get(variable))
        for estimate in estimates.get(variable, []):
            line = ax.axvline(estimate, color="black", linestyle=":", linewidth=1.0)
            line.set_gid("parameter-estimate")

    for index in range(len(variables), nrows * ncols):
        ax = fig.add_subplot(grid[index // ncols, index % ncols])
        ax.axis("off")

    cax = fig.add_subplot(grid[-1, :])
    sm = ScalarMappable(norm=RHAT_NORM, cmap=RHAT_CMAP)
    colorbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    colorbar.set_label("Rhat", fontsize=8)
    colorbar.set_ticks([0.99, 1.5, 2.0, 2.5, 3.0])
    colorbar.ax.tick_params(labelsize=7)


def draw_loo(ax: Axes, loo: pd.DataFrame) -> None:
    """Draw horizontal ELPD-LOO bars and standard-error intervals."""
    ax.clear()
    loo = _frame(loo)
    required = {"model", "elpd_loo"}
    if not required.issubset(loo.columns):
        ax.set_xlabel("ELPD LOO")
        ax.text(
            0.5,
            0.5,
            "No LOO diagnostics",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return

    values = pd.to_numeric(loo["elpd_loo"], errors="coerce").to_numpy(
        dtype=float, na_value=np.nan
    )
    finite = np.isfinite(values)
    rows = loo.loc[finite].copy()
    if rows.empty:
        ax.set_xlabel("ELPD LOO")
        ax.text(
            0.5,
            0.5,
            "No LOO diagnostics",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return

    model_names = rows["model"].astype(str).tolist()
    present = set(model_names)
    order = [model for model in MODE_ORDER if model in present]
    order.extend(model for model in model_names if model not in order)
    order = list(dict.fromkeys(order))
    row_by_model = {str(row["model"]): row for _, row in rows.iterrows()}
    y_positions = np.arange(len(order), dtype=float)
    p_values = np.array(
        [
            pd.to_numeric(
                pd.Series([row_by_model[model].get("p_loo", np.nan)]), errors="coerce"
            ).iloc[0]
            for model in order
        ],
        dtype=float,
    )
    valid_p = np.isfinite(p_values)
    colors = [UNKNOWN_COLOR] * len(order)
    if valid_p.any():
        p_min = float(p_values[valid_p].min())
        p_max = float(p_values[valid_p].max())
        if p_min == p_max:
            pad = max(abs(p_min) * 0.01, 0.5)
            p_min -= pad
            p_max += pad
        norm = Normalize(vmin=p_min, vmax=p_max)
        colors = [
            LOO_CMAP(norm(value)) if np.isfinite(value) else UNKNOWN_COLOR
            for value in p_values
        ]

    center = np.array([float(row_by_model[model]["elpd_loo"]) for model in order])
    errors = np.array(
        [
            pd.to_numeric(
                pd.Series([row_by_model[model].get("se_elpd_loo", np.nan)]),
                errors="coerce",
            ).iloc[0]
            for model in order
        ],
        dtype=float,
    )
    errors = np.where(np.isfinite(errors) & (errors >= 0), errors, 0.0)
    ax.barh(y_positions, center, height=0.6, color=colors, zorder=2)
    ax.errorbar(
        center,
        y_positions,
        xerr=errors,
        fmt="none",
        ecolor="black",
        elinewidth=0.8,
        capsize=2,
        zorder=3,
    )
    ax.set_yticks(y_positions, labels=order)
    ax.invert_yaxis()
    ax.yaxis.tick_right()
    ax.set_ylabel("")
    ax.set_xlabel("ELPD LOO", fontsize=9, labelpad=2)
    ax.grid(axis="x", color="grey", linestyle=":", linewidth=0.6, zorder=0)
    ax.tick_params(axis="y", length=0, labelsize=8)
    ax.tick_params(axis="x", labelsize=8)

    if valid_p.any():
        sm = ScalarMappable(norm=norm, cmap=LOO_CMAP)
        colorbar = ax.figure.colorbar(
            sm, ax=ax, orientation="horizontal", pad=0.32, fraction=0.08
        )
        colorbar.set_label(r"$p_{LOO}$", fontsize=8)
        colorbar.ax.tick_params(labelsize=7)
