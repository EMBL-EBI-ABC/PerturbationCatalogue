"""Matplotlib plots for validating MAVE GMM scores against known variants.

The score histogram approximates ggplot's pseudo-log transform with Matplotlib's
asinh-like symlog transform and bins uniformly in that transformed space. The
density curves use a normalized Gaussian KDE approximation.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PLOT_NAMES = (
    "distributions",
    "control_distribution",
    "control_roc_curves",
    "control_auc_distribution",
    "clinvar_distribution",
    "clinvar_roc_curves",
    "clinvar_auc_distribution",
    "gwas_distribution",
    "gwas_roc_curves",
    "gwas_auc_distribution",
)

SCORE_ORDER = ("s", "ns", "w", "aw", "p", "ap")
CONTROL_SCORE_ORDER = ("s", "ns", "w", "p")
SCORE_LABELS = {
    "s": "Raw Score (s)",
    "ns": "Median Scaled Score",
    "w": "Standardised Score (w)",
    "aw": "Abs. Std. Score (|w - 1|)",
    "p": "Neutral SD Scaled w",
    "ap": "Abs. SD scaled w",
}
SCORE_COLORS = {
    "ko": "#e41a1c",
    "ko/del": "#fb9a99",
    "del": "#ff7f00",
    "del/neut": "#a6cee3",
    "neut": "#377eb8",
    "neut/pos": "#cab2d6",
    "pos": "#984ea3",
    "noise": "#999999",
    "incoherent": "#4daf4a",
}
CONTROL_COLORS = {
    "WT": "#33a02c",
    "Synonymous": "#b2df8a",
    "Nonsense": "#e41a1c",
}
CLINVAR_COLORS = {
    "benign": "#1f78b4",
    "likely benign": "#80b1d3",
    "uncertain significance": "#984ea3",
    "likely pathogenic": "#fb9a99",
    "pathogenic": "#e31a1c",
    "conflicting classifications of pathogenicity": "#a65628",
}
UNKNOWN_COLOR = "#888888"
_MODE_ORDER = ("ko", "del", "neut", "pos", "noise")


def _frame(tables: dict[str, pd.DataFrame], name: str) -> pd.DataFrame:
    value = tables.get(name)
    return value.copy() if isinstance(value, pd.DataFrame) else pd.DataFrame()


def _numeric(data: pd.DataFrame, column: str) -> np.ndarray:
    if column not in data:
        return np.full(len(data), np.nan, dtype=float)
    return pd.to_numeric(data[column], errors="coerce").to_numpy(
        dtype=float, na_value=np.nan
    )


def _text(data: pd.DataFrame, column: str, default: str = "") -> np.ndarray:
    if column not in data:
        return np.full(len(data), default, dtype=object)
    return data[column].astype("string").fillna(default).to_numpy(dtype=object)


def _truthy(value: Any) -> bool:
    if value is None:
        return False
    missing = pd.isna(value)
    if isinstance(missing, (bool, np.bool_)) and missing:
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(value) if np.isfinite(value) else False
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def _mode_values(data: pd.DataFrame) -> np.ndarray:
    values = _text(data, "mode_density", "unclassified")
    return np.asarray(
        [value if value else "unclassified" for value in values], dtype=object
    )


def _ordered(values: Sequence[Any], preferred: Sequence[str] = ()) -> list[str]:
    present = {
        str(value)
        for value in values
        if value is not None and not pd.isna(value) and str(value)
    }
    result = [value for value in preferred if value in present]
    if preferred is _MODE_ORDER or tuple(preferred) == _MODE_ORDER:
        result.extend(value for value in SCORE_COLORS if value in present - set(result))
    result.extend(sorted(present - set(result)))
    return result


def _mode_color(label: str) -> str:
    return SCORE_COLORS.get(label, UNKNOWN_COLOR)


def _empty_axis(ax: Axes, message: str) -> None:
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=8,
        color="#555555",
    )


def _draw_score_histogram(ax: Axes, scores: pd.DataFrame, column: str) -> None:
    values = _numeric(scores, column)
    modes = _mode_values(scores)
    keep = np.isfinite(values)
    values = values[keep]
    modes = modes[keep]
    if values.size == 0:
        _empty_axis(ax, "No finite scores")
        return

    categories = _ordered(modes.tolist(), _MODE_ORDER)
    grouped = [values[modes == category] for category in categories]
    colors = [_mode_color(category) for category in categories]
    if column == "s":
        # Use equal-width bins after the symlog transform, as the R plot does
        # after its pseudo-log scale transform.
        transformed_edges = np.linspace(
            np.arcsinh(values.min()), np.arcsinh(values.max()), 51
        )
        if transformed_edges[0] == transformed_edges[-1]:
            transformed_edges = np.linspace(
                transformed_edges[0] - 0.5, transformed_edges[-1] + 0.5, 51
            )
        bins = np.sinh(transformed_edges)
    else:
        bins = 50
    ax.hist(
        grouped,
        bins=bins,
        stacked=True,
        color=colors,
        edgecolor="white",
        linewidth=0.25,
        label=categories,
    )
    ax.set_ylabel("Count")
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.set_axisbelow(True)
    if categories:
        ax.legend(
            handles=[
                Patch(facecolor=_mode_color(category), label=category)
                for category in categories
            ],
            loc="upper right",
            frameon=False,
            fontsize=7,
        )


def _distributions(scores: pd.DataFrame) -> Figure:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), squeeze=False)
    try:
        raw_ax, fitted_ax = axes[0]
        _draw_score_histogram(raw_ax, scores, "s")
        raw_ax.set_xscale("symlog", linthresh=1.0)
        raw_ax.set_xlabel("Raw Scores (s)")
        raw_ax.set_title("Raw scores")

        _draw_score_histogram(fitted_ax, scores, "w")
        fitted_ax.set_xlabel("Fitted Scores (w)")
        fitted_ax.set_title("Fitted scores")
        fig.suptitle("MAVE score distributions")
        fig.tight_layout()
        return fig
    except Exception:
        plt.close(fig)
        raise


def _background_histogram(ax: Axes, scores: pd.DataFrame) -> list[Patch]:
    values = _numeric(scores, "w")
    modes = _mode_values(scores)
    keep = np.isfinite(values) & (values >= 0.0) & (values <= 2.0)
    edges = np.arange(0.0, 2.0001, 0.05)
    centers = (edges[:-1] + edges[1:]) / 2.0
    categories = _ordered(modes[keep].tolist(), _MODE_ORDER)
    counts = {
        category: np.histogram(values[keep & (modes == category)], bins=edges)[0]
        for category in categories
    }
    global_max = max((int(group.max()) for group in counts.values()), default=0)
    if global_max <= 0:
        _empty_axis(ax, "No background scores in 0–2")
        return []

    bottom = np.zeros(centers.size, dtype=float)
    handles: list[Patch] = []
    for category in categories:
        heights = counts[category].astype(float) / global_max
        ax.bar(
            centers,
            heights,
            bottom=bottom,
            width=0.05,
            align="center",
            color=_mode_color(category),
            alpha=0.25,
            edgecolor="none",
            zorder=1,
        )
        bottom += heights
        handles.append(
            Patch(
                facecolor=_mode_color(category),
                alpha=0.35,
                label=category,
            )
        )
    return handles


def _draw_scaled_density(
    ax: Axes, values: np.ndarray, *, color: str, label: str
) -> tuple[Line2D | None, str | None]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return None, None

    unique = np.unique(values)
    if unique.size < 2:
        rug_y = np.full(values.size, 0.012)
        (line,) = ax.plot(
            values,
            rug_y,
            linestyle="none",
            marker="|",
            markersize=9,
            markeredgewidth=1.2,
            color=color,
            label=label,
            zorder=4,
        )
        return line, f"{label}: singular sample shown as rug (n={values.size})"

    standard_deviation = float(np.std(values, ddof=1))
    q25, q75 = np.quantile(values, [0.25, 0.75])
    robust_scale = float((q75 - q25) / 1.34)
    scale = (
        min(standard_deviation, robust_scale)
        if robust_scale > 0
        else standard_deviation
    )
    bandwidth = max(0.9 * scale * values.size ** (-0.2), np.finfo(float).eps)
    grid = np.linspace(
        float(values.min() - 3 * bandwidth),
        float(values.max() + 3 * bandwidth),
        256,
    )
    density = np.zeros(grid.size, dtype=float)
    for start in range(0, values.size, 512):
        chunk = values[start : start + 512]
        standardized = (grid[:, None] - chunk[None, :]) / bandwidth
        density += np.exp(-0.5 * standardized**2).sum(axis=1)
    peak = float(density.max())
    if not np.isfinite(peak) or peak <= 0:
        return None, f"{label}: density could not be estimated (n={values.size})"
    (line,) = ax.plot(
        grid,
        density / peak,
        color=color,
        linewidth=1.35,
        label=label,
        zorder=3,
    )
    return line, None


def _background_distribution(
    scores: pd.DataFrame,
    values_table: pd.DataFrame,
    *,
    title: str,
    value_column: str,
    category_column: str | None = None,
    category_order: Sequence[str] = (),
    category_colors: dict[str, str] | None = None,
    hit_only: bool = False,
) -> Figure:
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    try:
        handles = _background_histogram(ax, scores)
        notes: list[str] = []
        value_values = _numeric(values_table, value_column)
        eligible = np.isfinite(value_values)
        if hit_only:
            hit_values = (
                values_table["gwas_hit"].map(_truthy).to_numpy()
                if "gwas_hit" in values_table
                else np.zeros(len(values_table), dtype=bool)
            )
            eligible &= hit_values
            if not eligible.any():
                notes.append("No GWAS hits with finite fitted scores")
        categories = (
            _ordered(_text(values_table, category_column).tolist(), category_order)
            if category_column is not None
            else (["GWAS Hit"] if eligible.any() else [])
        )
        if category_column is not None and not categories:
            notes.append(f"No {title.split()[0]} categories with finite scores")

        for category in categories:
            selected = eligible.copy()
            if category_column is not None:
                selected &= _text(values_table, category_column) == category
            if not selected.any():
                continue
            values = value_values[selected]
            if category_colors is not None:
                color = category_colors.get(
                    category, category_colors.get(category.lower(), UNKNOWN_COLOR)
                )
            elif hit_only:
                color = "black"
            else:
                color = UNKNOWN_COLOR
            line, note = _draw_scaled_density(ax, values, color=color, label=category)
            if line is not None:
                handles.append(line)
            if note:
                notes.append(note)

        if hit_only and eligible.any():
            line = ax.axvline(
                1.0,
                color="#555555",
                linestyle="--",
                linewidth=0.8,
                label="Neutral score (w = 1)",
                zorder=2,
            )
            handles.append(line)
        for index, note in enumerate(notes):
            ax.text(
                0.01,
                0.99 - index * 0.06,
                note,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=7,
                color="#555555",
            )
        current_xlim = ax.get_xlim()
        ax.set_xlim(min(0.0, current_xlim[0]), max(2.0, current_xlim[1]))
        ax.set_xlabel("Fitted Scores (w)")
        ax.set_ylabel("Scaled Density")
        ax.set_title(title)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)
        ax.set_axisbelow(True)
        if handles:
            ax.legend(
                handles=handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.16),
                ncol=min(5, len(handles)),
                frameon=False,
                fontsize=7,
            )
        fig.text(
            0.5,
            0.01,
            "Histogram max scaled to 1; each density max scaled to 1.",
            ha="center",
            fontsize=7,
            color="#555555",
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        return fig
    except Exception:
        plt.close(fig)
        raise


def _control_distribution(scores: pd.DataFrame, controls: pd.DataFrame) -> Figure:
    plotted_frames = []
    for flag, label in (
        ("wt_mut", "WT"),
        ("synonymous_mut", "Synonymous"),
        ("nonsense_mut", "Nonsense"),
    ):
        if flag in controls and "w" in controls:
            selected = controls[flag].map(_truthy)
            plotted_frames.append(controls.loc[selected, ["w"]].assign(category=label))
    plotted = (
        pd.concat(plotted_frames, ignore_index=True)
        if plotted_frames
        else pd.DataFrame(columns=("w", "category"))
    )
    return _background_distribution(
        scores,
        plotted,
        title="Internal control score distributions",
        value_column="w",
        category_column="category",
        category_order=("WT", "Synonymous", "Nonsense"),
        category_colors=CONTROL_COLORS,
    )


def _clinvar_distribution(scores: pd.DataFrame, matches: pd.DataFrame) -> Figure:
    matches = matches.copy()
    if "clnsig" in matches:
        matches["_plot_category"] = _text(matches, "clnsig", "Unknown")
    else:
        matches["_plot_category"] = ""
    return _background_distribution(
        scores,
        matches,
        title="ClinVar fitted score distributions",
        value_column="w",
        category_column="_plot_category",
        category_order=(
            "benign",
            "likely benign",
            "uncertain significance",
            "likely pathogenic",
            "pathogenic",
        ),
        category_colors=CLINVAR_COLORS,
    )


def _gwas_distribution(scores: pd.DataFrame, matches: pd.DataFrame) -> Figure:
    return _background_distribution(
        scores,
        matches,
        title="GWAS hit fitted score distribution",
        value_column="w",
        hit_only=True,
    )


def _benchmark_roc(tables: dict[str, pd.DataFrame], benchmark: str) -> pd.DataFrame:
    source = {"controls": "ctrl_roc", "clinvar": "clinvar_roc", "gwas": "gwas_roc"}[
        benchmark
    ]
    return _frame(tables, source)


def _draw_curve(
    ax: Axes,
    rows: pd.DataFrame,
    *,
    metric: str,
    benchmark: str,
) -> bool:
    x_column, y_column, area_column = (
        ("fpr", "tpr", "auc") if metric == "roc" else ("tpr", "precision", "pr_auc")
    )
    x = _numeric(rows, x_column)
    y = _numeric(rows, y_column)
    finite = np.isfinite(x) & np.isfinite(y)
    rows = rows.loc[finite].copy()
    x = x[finite]
    y = y[finite]
    if rows.empty:
        _empty_axis(ax, f"No {metric.upper()} data")
        return False

    urns = _text(rows, "score_urn")
    order = [value for value in _ordered(urns.tolist()) if value.lower() != "all"]
    if np.any(np.char.lower(urns.astype(str)) == "all"):
        order.append("all")
    plotted = False
    for urn in order:
        selected = (
            np.char.lower(urns.astype(str)) == "all" if urn == "all" else urns == urn
        )
        group = rows.loc[selected]
        group_x = _numeric(group, x_column)
        group_y = _numeric(group, y_column)
        valid = np.isfinite(group_x) & np.isfinite(group_y)
        group_x = group_x[valid]
        group_y = group_y[valid]
        if not group_x.size:
            continue
        curve_order = np.argsort(group_x, kind="stable")
        pooled = urn.lower() == "all"
        ax.plot(
            group_x[curve_order],
            group_y[curve_order],
            color="black" if pooled else "lightgrey",
            linestyle="-" if pooled else ":",
            linewidth=1.4 if pooled else 0.8,
            alpha=1.0 if pooled else 0.9,
            zorder=3 if pooled else 2,
        )
        plotted = True
        if pooled and area_column in group:
            area_values = _numeric(group, area_column)
            area_values = area_values[np.isfinite(area_values)]
            if area_values.size:
                metric_label = "AUC" if metric == "roc" else "PR AUC"
                ax.text(
                    0.97,
                    0.04,
                    f"{metric_label} = {area_values[0]:.3g}",
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=7,
                    color="black",
                )
    if not plotted:
        _empty_axis(ax, f"No finite {metric.upper()} curves")
    else:
        if metric == "roc":
            ax.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=0.7)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(color="#dddddd", linewidth=0.5)
    if metric == "roc":
        ax.set_xlabel("1 − Specificity / FPR", fontsize=7)
        ax.set_ylabel("Sensitivity / TPR", fontsize=7)
        ax.set_aspect("equal", adjustable="box")
    else:
        ax.set_xlabel("Recall / Sensitivity / TPR", fontsize=7)
        ax.set_ylabel("Precision / PPV", fontsize=7)
    ax.tick_params(labelsize=7)
    return plotted


def _roc_figure(roc: pd.DataFrame, benchmark: str) -> Figure:
    score_types = CONTROL_SCORE_ORDER if benchmark == "controls" else SCORE_ORDER
    columns = len(score_types)
    width = max(10.0, 2.8 * columns)
    fig, axes = plt.subplots(
        2,
        columns,
        figsize=(width, 6.8),
        squeeze=False,
        sharex=False,
        sharey=False,
    )
    try:
        plotted_any = False
        for column, score_type in enumerate(score_types):
            rows = roc.loc[_text(roc, "score_type") == score_type].copy()
            roc_ax, pr_ax = axes[:, column]
            plotted_any |= _draw_curve(roc_ax, rows, metric="roc", benchmark=benchmark)
            plotted_any |= _draw_curve(pr_ax, rows, metric="pr", benchmark=benchmark)
            roc_ax.set_title(SCORE_LABELS[score_type], fontsize=8)
        handles = [
            Line2D([0], [0], color="lightgrey", linestyle=":", label="Per Study"),
            Line2D([0], [0], color="black", linestyle="-", label="All Studies"),
        ]
        if plotted_any:
            fig.legend(
                handles=handles,
                loc="lower center",
                ncol=2,
                frameon=False,
                fontsize=8,
            )
        fig.suptitle(f"{benchmark.title()} ROC and precision–recall curves")
        fig.tight_layout(rect=(0, 0.06, 1, 0.94), w_pad=0.8, h_pad=1.1)
        return fig
    except Exception:
        plt.close(fig)
        raise


def _auc_summary(
    tables: dict[str, pd.DataFrame], benchmark: str, roc: pd.DataFrame
) -> pd.DataFrame:
    summary = _frame(tables, "auc_summary")
    if {"benchmark", "score_type", "score_urn"}.issubset(summary.columns):
        summary = summary.loc[summary["benchmark"].astype("string") == benchmark]
        if not summary.empty:
            return summary.copy()
    required = {"score_type", "score_urn"}
    if not required.issubset(roc.columns):
        return pd.DataFrame(columns=("score_type", "score_urn", "auc", "pr_auc"))
    columns = [
        column
        for column in ("score_type", "score_urn", "auc", "pr_auc")
        if column in roc
    ]
    return roc[columns].drop_duplicates().copy()


def _auc_distribution(
    tables: dict[str, pd.DataFrame], benchmark: str, roc: pd.DataFrame
) -> Figure:
    summary = _auc_summary(tables, benchmark, roc)
    score_types = CONTROL_SCORE_ORDER if benchmark == "controls" else SCORE_ORDER
    metrics = (("auc", "Sensitivity–Specificity"), ("pr_auc", "Precision–Recall"))
    fig, axes = plt.subplots(
        len(score_types),
        len(metrics),
        figsize=(8.5, max(6.0, 1.45 * len(score_types))),
        squeeze=False,
        sharex=True,
    )
    try:
        score_values = _text(summary, "score_type")
        urn_values = _text(summary, "score_urn")
        for row_index, score_type in enumerate(score_types):
            for column_index, (metric, metric_title) in enumerate(metrics):
                ax = axes[row_index, column_index]
                if row_index == 0:
                    ax.set_title(metric_title, fontsize=8)
                if column_index == 0:
                    ax.set_ylabel(SCORE_LABELS[score_type], fontsize=7)
                if metric not in summary:
                    _empty_axis(ax, "No AUC data")
                    ax.set_xlim(0, 1)
                    continue
                values = _numeric(summary, metric)
                selected = (score_values == score_type) & np.isfinite(values)
                selected_urns = urn_values[selected]
                selected_values = values[selected]
                pooled_mask = np.char.lower(selected_urns.astype(str)) == "all"
                study_values = selected_values[~pooled_mask]
                pooled_values = selected_values[pooled_mask]
                if study_values.size:
                    ax.hist(
                        study_values,
                        bins=np.linspace(0, 1, 51),
                        color="#999999",
                        edgecolor="white",
                        linewidth=0.25,
                    )
                else:
                    ax.text(
                        0.5,
                        0.62,
                        "No per-study AUC data",
                        transform=ax.transAxes,
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="#555555",
                    )
                if pooled_values.size:
                    ax.axvline(
                        pooled_values[0],
                        color="black",
                        linestyle=":",
                        linewidth=1.2,
                    )
                else:
                    ax.text(
                        0.5,
                        0.43,
                        "No pooled AUC data",
                        transform=ax.transAxes,
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="#555555",
                    )
                ax.set_xlim(0, 1)
                ax.grid(axis="y", color="#dddddd", linewidth=0.5)
                ax.set_axisbelow(True)
                ax.tick_params(labelsize=7)
                if row_index == len(score_types) - 1:
                    ax.set_xlabel("AUC", fontsize=7)
        fig.suptitle(f"{benchmark.title()} AUC distributions")
        fig.tight_layout(rect=(0, 0, 1, 0.96), h_pad=0.8, w_pad=0.8)
        return fig
    except Exception:
        plt.close(fig)
        raise


def build_plots(tables: dict[str, pd.DataFrame]) -> dict[str, Figure]:
    """Build the ten validation figures from their input tables."""
    scores = _frame(tables, "scores")
    controls = _frame(tables, "controls")
    clinvar = _frame(tables, "clinvar_matches")
    gwas = _frame(tables, "gwas_matches")
    created: dict[str, Figure] = {}
    try:
        created["distributions"] = _distributions(scores)
        created["control_distribution"] = _control_distribution(scores, controls)
        ctrl_roc = _benchmark_roc(tables, "controls")
        created["control_roc_curves"] = _roc_figure(ctrl_roc, "controls")
        created["control_auc_distribution"] = _auc_distribution(
            tables, "controls", ctrl_roc
        )
        created["clinvar_distribution"] = _clinvar_distribution(scores, clinvar)
        clinvar_roc = _benchmark_roc(tables, "clinvar")
        created["clinvar_roc_curves"] = _roc_figure(clinvar_roc, "clinvar")
        created["clinvar_auc_distribution"] = _auc_distribution(
            tables, "clinvar", clinvar_roc
        )
        created["gwas_distribution"] = _gwas_distribution(scores, gwas)
        gwas_roc = _benchmark_roc(tables, "gwas")
        created["gwas_roc_curves"] = _roc_figure(gwas_roc, "gwas")
        created["gwas_auc_distribution"] = _auc_distribution(tables, "gwas", gwas_roc)
        return created
    except Exception:
        for figure in created.values():
            plt.close(figure)
        raise


def save_plots(
    plots: dict[str, Figure],
    output_dir: Path,
    formats: Sequence[str] = ("png", "pdf"),
    dpi: int = 150,
) -> list[Path]:
    """Save each named figure in every requested format and return the paths."""
    if dpi <= 0:
        raise ValueError("dpi must be a positive integer")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for name, figure in plots.items():
        for requested_format in formats:
            file_format = str(requested_format).strip().lower().lstrip(".")
            if not file_format:
                raise ValueError("formats cannot contain an empty file extension")
            path = output_dir / f"{name}.{file_format}"
            figure.savefig(path, format=file_format, dpi=dpi, bbox_inches="tight")
            paths.append(path)
    return paths
