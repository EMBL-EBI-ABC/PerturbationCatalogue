#!/usr/bin/env python3
"""Plot quality-control summaries for MaveDB experiments.

The data preparation code lives in :mod:`_mave_qc`; this module owns the
matplotlib figures and the command-line entry point.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy.stats import gaussian_kde

_HERE = Path(__file__).resolve().parent
_MAVEDB_DIR = _HERE.parent
_DEFAULT_SCORES = _MAVEDB_DIR / "output" / "processed" / "scores.tsv"
_DEFAULT_METADATA = _MAVEDB_DIR / "output" / "processed" / "metadata.tsv"
_DEFAULT_SINGLES = _MAVEDB_DIR / "output" / "processed" / "maves" / "singles.tsv"
_DEFAULT_OUTPUT = _MAVEDB_DIR / "output" / "figures" / "maves" / "qc"

_MUTATION_ORDER = ("Single", "Multiple", "Both", "Mixed", "Neither??")
_MUTATION_COLORS = {
    "Single": "#4C78A8",
    "Multiple": "#F58518",
    "Both": "#54A24B",
    "Mixed": "#B279A2",
    "Neither??": "#9D9D9D",
    "Unknown": "#BAB0AC",
}
_MATCH_ORDER = ("Same Experiment", "Same Study", "Same Gene")
_MATCH_COLORS = {
    "Same Experiment": "#4C78A8",
    "Same Study": "#F58518",
    "Same Gene": "#54A24B",
}
_HUMAN_ORDER = ("Human", "Other")
_METRIC_LABELS = {
    "completeness": "Substitution completeness",
    "completeness_missense": "Missense completeness",
    "stop_proportion": "Nonsense proportion",
    "del_proportion": "Deletion proportion",
    "synonymous_proportion": "Synonymous proportion",
    "wt_mutants": "WT variants",
    "coverage_50": ">50% missense coverage",
    "coverage_75": ">75% missense coverage",
    "coverage_100": "100% missense coverage",
    "coverage_nonsense": "Nonsense coverage",
    "coverage_synonymous": "Synonymous coverage",
}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#E5E7EB",
        "grid.linewidth": 0.7,
        "axes.axisbelow": True,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def _frame(data: Mapping[str, pd.DataFrame], key: str) -> pd.DataFrame:
    value = data.get(key)
    return value if isinstance(value, pd.DataFrame) else pd.DataFrame()


def _values(frame: pd.DataFrame, column: str, default: object = np.nan) -> pd.Series:
    if column in frame.columns:
        return frame[column]
    return pd.Series(default, index=frame.index)


def _numeric(frame: pd.DataFrame, column: str) -> np.ndarray:
    values = pd.to_numeric(_values(frame, column), errors="coerce")
    return values.to_numpy(dtype=float, na_value=np.nan)


def _finite(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def _nonempty_axis(ax: Axes, message: str = "No data available") -> None:
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        color="#6B7280",
        transform=ax.transAxes,
    )
    ax.set_axis_off()


def _histogram(
    ax: Axes,
    values: np.ndarray,
    *,
    bins: int = 30,
    log_x: bool = False,
    color: str = "#4C78A8",
    label: str | None = None,
    alpha: float = 0.85,
) -> bool:
    values = _finite(values)
    if not len(values):
        _nonempty_axis(ax)
        return False

    bin_spec: int | np.ndarray = bins
    if log_x and np.all(values > 0):
        low, high = float(values.min()), float(values.max())
        if high <= low:
            low, high = low / 1.2, high * 1.2
        bin_spec = np.geomspace(low, high, bins + 1)

    ax.hist(values, bins=bin_spec, color=color, alpha=alpha, label=label)
    if log_x:
        if np.all(values > 0):
            ax.set_xscale("log")
        elif np.all(values >= 0):
            ax.set_xscale("symlog", linthresh=1)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    return True


def _grouped_histogram(
    ax: Axes,
    frame: pd.DataFrame,
    value_column: str,
    group_column: str,
    groups: Sequence[str],
    colors: Mapping[str, str],
    *,
    bins: int = 30,
    log_x: bool = False,
) -> bool:
    arrays: list[np.ndarray] = []
    labels: list[str] = []
    for group in groups:
        mask = (
            _values(frame, group_column, "Unknown")
            .fillna("Unknown")
            .astype(str)
            .eq(group)
        )
        values = _finite(_numeric(frame.loc[mask], value_column))
        if len(values):
            arrays.append(values)
            labels.append(group)
    if not arrays:
        _nonempty_axis(ax)
        return False

    all_values = np.concatenate(arrays)
    bin_spec: int | np.ndarray = bins
    if log_x and np.all(all_values > 0):
        low, high = float(all_values.min()), float(all_values.max())
        if high <= low:
            low, high = low / 1.2, high * 1.2
        bin_spec = np.geomspace(low, high, bins + 1)
    ax.hist(
        arrays,
        bins=bin_spec,
        stacked=True,
        color=[colors.get(group, colors.get("Unknown", "#9D9D9D")) for group in labels],
        label=labels,
        alpha=0.88,
        edgecolor="white",
        linewidth=0.25,
    )
    if log_x:
        if np.all(all_values > 0):
            ax.set_xscale("log")
        elif np.all(all_values >= 0):
            ax.set_xscale("symlog", linthresh=1)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    return True


def _density_scatter(
    ax: Axes,
    x: np.ndarray,
    y: np.ndarray,
    *,
    x_log: bool = False,
    cmap: str = "YlOrRd",
    point_color: str = "#4C78A8",
    size: float = 10,
) -> bool:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if x_log:
        valid &= x > 0
    x, y = x[valid], y[valid]
    if not len(x):
        _nonempty_axis(ax)
        return False

    can_hexbin = len(x) >= 150 and np.unique(x).size > 1 and np.unique(y).size > 1
    if can_hexbin:
        bins = ax.hexbin(
            x,
            y,
            gridsize=38,
            mincnt=1,
            bins="log",
            xscale="log" if x_log else "linear",
            cmap=cmap,
            linewidths=0,
        )
        colorbar = ax.figure.colorbar(bins, ax=ax, pad=0.025, fraction=0.045)
        colorbar.set_label("Score sets per bin", fontsize=7)
        colorbar.ax.tick_params(labelsize=7)
    else:
        ax.scatter(x, y, s=size, alpha=0.6, color=point_color, edgecolors="none")
        if x_log and np.unique(x).size > 1:
            ax.set_xscale("log")
    return True


def _metric_grid(
    frame: pd.DataFrame,
    metric_order: Sequence[str],
    *,
    scatter: bool,
    x_label: str,
    y_label: str,
    title: str,
) -> Figure:
    rows = len(metric_order)
    fig, axes = plt.subplots(
        rows,
        2,
        figsize=(12, max(7.5, rows * 2.25)),
        squeeze=False,
        sharex=False,
        sharey=False,
    )
    human = _values(frame, "human", "Other").fillna("Other").astype(str)
    metrics = _values(frame, "metric", "").fillna("").astype(str)
    for row, metric in enumerate(metric_order):
        for col, label in enumerate(_HUMAN_ORDER):
            ax = axes[row, col]
            subset = frame.loc[metrics.eq(metric) & human.eq(label)]
            ax.set_title(f"{_METRIC_LABELS.get(metric, metric)} · {label}", loc="left")
            if scatter:
                _density_scatter(
                    ax,
                    _numeric(subset, "n"),
                    _numeric(subset, "value"),
                    x_log=True,
                )
                if ax.axison:
                    ax.set_xlabel(x_label)
                    ax.set_ylabel(y_label)
                    ax.set_xscale("log")
            else:
                _histogram(ax, _numeric(subset, "value"), bins=35, color="#4C78A8")
                if ax.axison:
                    ax.set_xlabel(x_label)
                    ax.set_ylabel(y_label)
    fig.suptitle(title, x=0.06, y=0.995, ha="left", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=0.8, w_pad=1.2)
    return fig


def _draw_density_curve(ax: Axes, values: np.ndarray, color: str, label: str) -> None:
    values = _finite(values)
    if not len(values):
        return
    if len(values) > 40_000:
        selected = np.linspace(0, len(values) - 1, 40_000, dtype=int)
        density_values = values[selected]
    else:
        density_values = values

    if len(density_values) > 1 and np.ptp(density_values) > 1e-10:
        try:
            density = gaussian_kde(density_values)
            grid = np.linspace(-1, 1, 301)
            ax.plot(grid, density(grid), color=color, linewidth=1.6, label=label)
            return
        except (np.linalg.LinAlgError, ValueError):
            pass

    ax.hist(
        values,
        bins=np.linspace(-1, 1, 46),
        density=True,
        histtype="step",
        linewidth=1.5,
        color=color,
        label=label,
    )


def _correlation_plot(frame: pd.DataFrame, *, human_only: bool) -> Figure:
    plot_data = frame.copy()
    if human_only:
        plot_data = plot_data.loc[
            _values(plot_data, "organism", "").astype(str).eq("Homo sapiens")
        ]
    spearman = _numeric(plot_data, "spearman")
    pearson = _numeric(plot_data, "pearson")
    match = _values(plot_data, "match_type", "Unknown").fillna("Unknown").astype(str)
    valid = np.isfinite(spearman) & np.isfinite(pearson)
    spearman, pearson, match = (
        spearman[valid],
        pearson[valid],
        match.iloc[np.flatnonzero(valid)],
    )

    fig = plt.figure(figsize=(11, 8))
    grid = fig.add_gridspec(2, 2, height_ratios=(1, 1.15), hspace=0.35, wspace=0.28)
    ax_scatter = fig.add_subplot(grid[0, :])
    ax_spear = fig.add_subplot(grid[1, 0])
    ax_pear = fig.add_subplot(grid[1, 1])
    if len(spearman):
        max_points = 100_000
        if len(spearman) > max_points:
            selected = np.linspace(0, len(spearman) - 1, max_points, dtype=int)
            spearman_plot, pearson_plot = spearman[selected], pearson[selected]
            match_plot = match.iloc[selected]
        else:
            spearman_plot, pearson_plot, match_plot = spearman, pearson, match
        handles = []
        for category in _MATCH_ORDER:
            selected = match_plot.eq(category).to_numpy()
            if not selected.any():
                continue
            ax_scatter.scatter(
                spearman_plot[selected],
                pearson_plot[selected],
                s=8,
                alpha=0.45,
                color=_MATCH_COLORS[category],
                edgecolors="none",
                label=category,
                rasterized=True,
            )
            handles.append(Patch(color=_MATCH_COLORS[category], label=category))
        unknown = ~match_plot.isin(_MATCH_ORDER).to_numpy()
        if unknown.any():
            ax_scatter.scatter(
                spearman_plot[unknown],
                pearson_plot[unknown],
                s=8,
                alpha=0.45,
                color="#9D9D9D",
                edgecolors="none",
                label="Other",
                rasterized=True,
            )
            handles.append(Patch(color="#9D9D9D", label="Other"))
        ax_scatter.plot(
            [-1, 1], [-1, 1], color="#9CA3AF", linewidth=0.8, linestyle="--"
        )
        ax_scatter.set(xlim=(-1.03, 1.03), ylim=(-1.03, 1.03))
        ax_scatter.set_xlabel("Spearman correlation")
        ax_scatter.set_ylabel("Pearson correlation")
        if handles:
            ax_scatter.legend(
                handles=handles, title="Comparison", frameon=False, ncol=3
            )
        for category in _MATCH_ORDER:
            values = spearman[match.eq(category).to_numpy()]
            _draw_density_curve(ax_spear, values, _MATCH_COLORS[category], category)
            values = pearson[match.eq(category).to_numpy()]
            _draw_density_curve(ax_pear, values, _MATCH_COLORS[category], category)
        for ax, label in ((ax_spear, "Spearman"), (ax_pear, "Pearson")):
            ax.set_xlim(-1, 1)
            ax.set_xlabel(f"{label} correlation")
            ax.set_ylabel("Density")
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax_pear.legend(title="Comparison", frameon=False)
    else:
        _nonempty_axis(ax_scatter)
        _nonempty_axis(ax_spear)
        _nonempty_axis(ax_pear)

    suffix = " · Human score sets" if human_only else ""
    fig.suptitle(
        f"Within-target score-set correlation{suffix}",
        x=0.06,
        y=0.99,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.92)
    return fig


def _plot_datasets(frame: pd.DataFrame) -> Figure:
    names = (
        "organism",
        "target_name",
        "experiment_set_urn",
        "experiment_urn",
        "score_urn",
        "meta_score_urn",
    )
    labels = (
        "Organisms",
        "Targets",
        "Experiment sets",
        "Experiments",
        "Score sets",
        "Meta-analyses",
    )
    fig, ax = plt.subplots(figsize=(10, 5.8))
    if frame.empty or not {"name", "mutation_type", "n"}.issubset(frame.columns):
        _nonempty_axis(ax)
    else:
        table = frame.copy()
        table["n"] = pd.to_numeric(table["n"], errors="coerce").fillna(0)
        table = table.pivot_table(
            index="name",
            columns="mutation_type",
            values="n",
            aggfunc="sum",
            fill_value=0,
        ).reindex(index=names, fill_value=0)
        bottom = np.zeros(len(names), dtype=float)
        order = [name for name in _MUTATION_ORDER if name in table.columns]
        for mutation_type in order:
            values = table[mutation_type].to_numpy(dtype=float)
            ax.bar(
                np.arange(len(names)),
                values,
                bottom=bottom,
                width=0.68,
                color=_MUTATION_COLORS[mutation_type],
                label=mutation_type,
            )
            bottom += values
        ax.set_xticks(np.arange(len(names)), labels, rotation=18, ha="right")
        ax.set_ylabel("Distinct catalogue entries")
        ax.legend(title="Mutation type", frameon=False, ncol=3, loc="upper left")
    ax.set_title("Catalogue entities by mutation type", loc="left", fontweight="bold")
    fig.tight_layout()
    return fig


def _organism_label(value: object) -> str:
    if pd.isna(value):
        return "Unknown"
    text = str(value)
    if text.startswith("["):
        try:
            decoded = json.loads(text)
        except json.JSONDecodeError:
            return text
        if isinstance(decoded, list):
            return "Unknown / multiple targets"
    return text


def _plot_organisms(frame: pd.DataFrame) -> Figure:
    fig_height = max(5, 0.24 * max(len(frame), 1) + 2.5)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    if frame.empty or not {"organism", "n"}.issubset(frame.columns):
        _nonempty_axis(ax)
    else:
        plot_data = frame.copy()
        plot_data["n"] = pd.to_numeric(plot_data["n"], errors="coerce")
        plot_data = plot_data.dropna(subset=["n"]).sort_values("n")
        if plot_data.empty:
            _nonempty_axis(ax)
        else:
            kingdoms = (
                _values(plot_data, "kingdom", "Other").fillna("Other").astype(str)
            )
            categories = list(dict.fromkeys(kingdoms))
            palette = plt.get_cmap("tab10")
            color_map = {
                name: palette(index % 10) for index, name in enumerate(categories)
            }
            colors = [color_map[name] for name in kingdoms]
            y = np.arange(len(plot_data))
            bars = ax.barh(y, plot_data["n"], color=colors, height=0.78)
            organism_labels = plot_data["organism"].map(_organism_label)
            ax.set_yticks(y, organism_labels)
            ax.set_xlabel("Score sets (log scale)")
            if np.all(plot_data["n"] > 0):
                ax.set_xscale("log")
            ax.grid(axis="y", visible=False)
            handles = [Patch(color=color_map[name], label=name) for name in categories]
            ax.legend(
                handles=handles, title="Kingdom", frameon=False, loc="lower right"
            )
            ax.set_title(
                "Organisms represented in MaveDB", loc="left", fontweight="bold"
            )
            for bar, value in zip(bars, plot_data["n"]):
                ax.annotate(
                    f"{value:g}",
                    (value, bar.get_y() + bar.get_height() / 2),
                    xytext=(3, 0),
                    textcoords="offset points",
                    va="center",
                    fontsize=7,
                )
    fig.tight_layout()
    return fig


def _plot_set_summary(frame: pd.DataFrame) -> Figure:
    specs = (
        ("experiment_set", "Experiments per experiment set", False),
        ("experiment", "Score sets per experiment", False),
        ("score", "Variants per score set", True),
    )
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), squeeze=False)
    for ax, (kind, title, log_x) in zip(axes[0], specs):
        subset = frame.loc[_values(frame, "type", "").astype(str).eq(kind)]
        _grouped_histogram(
            ax,
            subset,
            "n",
            "mut_types",
            _MUTATION_ORDER,
            _MUTATION_COLORS,
            bins=35,
            log_x=log_x,
        )
        ax.set_title(title, loc="left")
        if ax.axison:
            ax.set_xlabel("Count")
            ax.set_ylabel("Number of sets")
    handles = [Patch(color=_MUTATION_COLORS[key], label=key) for key in _MUTATION_ORDER]
    fig.legend(
        handles=handles,
        title="Mutation types",
        loc="lower center",
        ncol=5,
        frameon=False,
    )
    fig.suptitle(
        "Experiment and score-set sizes",
        x=0.06,
        y=0.99,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0.1, 1, 0.94), w_pad=1.2)
    return fig


def _plot_data_per_gene(frame: pd.DataFrame) -> Figure:
    columns = (
        ("experiment_sets", "Experiment sets per target", False),
        ("experiments", "Experiments per target", False),
        ("score_sets", "Score sets per target", False),
        ("meta_analyses", "Meta-analyses per target", False),
        ("measurements", "Measurements per target", True),
    )
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), squeeze=False)
    for ax, (column, title, log_x) in zip(axes.flat, columns):
        _histogram(ax, _numeric(frame, column), bins=30, log_x=log_x)
        ax.set_title(title, loc="left")
        if ax.axison:
            ax.set_xlabel("Count")
            ax.set_ylabel("Targets")
    axes.flat[-1].set_axis_off()
    fig.suptitle(
        "Data available per target",
        x=0.06,
        y=0.99,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=1.5, w_pad=1.2)
    return fig


def _plot_data_per_variant(frame: pd.DataFrame) -> Figure:
    columns = (
        ("experiment_sets", "Experiment sets per variant"),
        ("experiments", "Experiments per variant"),
        ("score_sets", "Score sets per variant"),
        ("measurements", "Measurements per variant"),
    )
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), squeeze=False)
    flags = _values(frame, "wt_mut", False)
    flag_labels = flags.map(
        lambda value: (
            "WT" if value is True or str(value).lower() in {"true", "1"} else "Mutant"
        )
    )
    plot_frame = frame.copy()
    plot_frame["_variant_class"] = flag_labels.to_numpy()
    color_map = {"WT": "#4C78A8", "Mutant": "#F58518"}
    for row, variant_class in enumerate(("WT", "Mutant")):
        subset = plot_frame.loc[plot_frame["_variant_class"].eq(variant_class)]
        for col, (column, title) in enumerate(columns):
            ax = axes[row, col]
            _histogram(
                ax,
                _numeric(subset, column),
                bins=28,
                color=color_map[variant_class],
            )
            ax.set_title(f"{title} · {variant_class}", loc="left")
            if ax.axison:
                ax.set_xlabel("Count")
                ax.set_ylabel("Variants")
    fig.suptitle(
        "Study coverage per protein variant",
        x=0.06,
        y=0.99,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=1.5, w_pad=1.2)
    return fig


def _plot_multiple_mutation(frame: pd.DataFrame) -> Figure:
    fig, ax = plt.subplots(figsize=(8, 6))
    _density_scatter(
        ax,
        _numeric(frame, "mean_muts"),
        _numeric(frame, "mean_backgrounds"),
    )
    ax.set_xlabel("Mean variants per measurement")
    ax.set_ylabel("Mean backgrounds per variant")
    ax.set_title("Multiple-mutation score-set summaries", loc="left", fontweight="bold")
    fig.tight_layout()
    return fig


def _plot_evolution(frame: pd.DataFrame) -> Figure:
    fig, ax = plt.subplots(figsize=(12, 6))
    if frame.empty or "estimate" not in frame:
        _nonempty_axis(ax)
    else:
        plot_data = frame.copy()
        plot_data["estimate"] = pd.to_numeric(plot_data["estimate"], errors="coerce")
        plot_data = plot_data.dropna(subset=["estimate"]).sort_values("estimate")
        if plot_data.empty:
            _nonempty_axis(ax)
        else:
            organism = _values(plot_data, "organism", "").fillna("").astype(str)
            kingdom = _values(plot_data, "kingdom", "Other").fillna("Other").astype(str)
            groups = np.where(organism.eq("Homo sapiens"), "Human", kingdom)
            categories = list(dict.fromkeys(groups))
            palette = plt.get_cmap("Dark2")
            colors = {
                group: palette(index % palette.N)
                for index, group in enumerate(categories)
            }
            x = np.arange(len(plot_data))
            ax.bar(
                x,
                plot_data["estimate"],
                color=[colors[group] for group in groups],
                width=0.88,
            )
            ax.axhline(0, color="#374151", linewidth=0.8)
            ax.set_ylim(-1.05, 1.05)
            ax.set_yticks(np.arange(-1, 1.01, 0.25))
            ax.set_ylabel("Correlation with BLOSUM100")
            ax.set_xlabel("Score sets (ordered by correlation)")
            ax.set_xticks([])
            handles = [Patch(color=colors[group], label=group) for group in categories]
            ax.legend(
                handles=handles,
                title="Organism group",
                frameon=False,
                ncol=min(4, len(handles)),
            )
            ax.set_title(
                "Study score correlation with evolutionary substitution scores",
                loc="left",
                fontweight="bold",
            )
    fig.tight_layout()
    return fig


def build_plots(data: Mapping[str, pd.DataFrame]) -> dict[str, Figure]:
    """Build all thirteen QC figures from prepared summary tables."""
    completeness_order = (
        "completeness",
        "completeness_missense",
        "stop_proportion",
        "del_proportion",
        "synonymous_proportion",
        "wt_mutants",
    )
    coverage_order = (
        "coverage_50",
        "coverage_75",
        "coverage_100",
        "coverage_nonsense",
        "coverage_synonymous",
    )
    return {
        "datasets": _plot_datasets(_frame(data, "datasets")),
        "organisms": _plot_organisms(_frame(data, "organisms")),
        "set_summary": _plot_set_summary(_frame(data, "set_summary")),
        "data_per_gene": _plot_data_per_gene(_frame(data, "data_per_gene")),
        "data_per_variant": _plot_data_per_variant(_frame(data, "data_per_variant")),
        "completeness_scatter": _metric_grid(
            _frame(data, "completeness"),
            completeness_order,
            scatter=True,
            x_label="Measurements per score set",
            y_label="Proportion / count",
            title="Score-set completeness",
        ),
        "completeness_hist": _metric_grid(
            _frame(data, "completeness"),
            completeness_order,
            scatter=False,
            x_label="Proportion / count",
            y_label="Score sets",
            title="Score-set completeness distributions",
        ),
        "coverage_scatter": _metric_grid(
            _frame(data, "coverage"),
            coverage_order,
            scatter=True,
            x_label="Measurements per score set",
            y_label="Proportion of positions",
            title="Score-set coverage",
        ),
        "coverage_hist": _metric_grid(
            _frame(data, "coverage"),
            coverage_order,
            scatter=False,
            x_label="Proportion of positions",
            y_label="Score sets",
            title="Score-set coverage distributions",
        ),
        "multiple_mutation": _plot_multiple_mutation(_frame(data, "multiple_mutation")),
        "score_set_correlation": _correlation_plot(
            _frame(data, "score_set_correlation"), human_only=False
        ),
        "score_set_correlation_human": _correlation_plot(
            _frame(data, "score_set_correlation"), human_only=True
        ),
        "blosum_study_correlation": _plot_evolution(
            _frame(data, "study_evo_correlation")
        ),
    }


def save_plots(
    plots: Mapping[str, Figure],
    output: Path,
    *,
    formats: Sequence[str] = ("png", "pdf"),
    dpi: int = 220,
) -> list[Path]:
    """Save figures, replacing files with matching figure names and formats."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
    for name, figure in plots.items():
        for file_format in formats:
            destination = output / f"{name}.{file_format}"
            figure.savefig(
                destination,
                format=file_format,
                dpi=dpi,
                bbox_inches="tight",
                metadata={"Creator": "MaveDB QC plotting pipeline"},
            )
            saved.append(destination)
    return saved


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scores", type=Path, default=_DEFAULT_SCORES, help="Complete score table TSV"
    )
    parser.add_argument(
        "--metadata", type=Path, default=_DEFAULT_METADATA, help="MaveDB metadata TSV"
    )
    parser.add_argument(
        "--singles",
        type=Path,
        default=_DEFAULT_SINGLES,
        help="Single-mutant score table TSV",
    )
    parser.add_argument(
        "--output", type=Path, default=_DEFAULT_OUTPUT, help="Directory for QC figures"
    )
    parser.add_argument(
        "--chunk-size", type=int, default=250_000, help="Score rows to read per chunk"
    )
    parser.add_argument(
        "--formats", nargs="+", choices=("png", "pdf"), default=("png", "pdf")
    )
    parser.add_argument("--dpi", type=int, default=220, help="PNG output resolution")
    parser.add_argument(
        "--matrix", type=Path, default=None, help="Optional BLOSUM100 matrix file"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parser().parse_args(argv)
    if args.chunk_size <= 0:
        raise SystemExit("--chunk-size must be a positive integer")
    if args.dpi <= 0:
        raise SystemExit("--dpi must be a positive integer")

    if __package__:
        from ._mave_qc import prepare_qc_data
    else:  # Direct execution: python /path/to/mave_qc.py
        from _mave_qc import prepare_qc_data

    data = prepare_qc_data(
        args.scores,
        args.metadata,
        args.singles,
        chunk_size=args.chunk_size,
        matrix_path=args.matrix,
    )
    plots = build_plots(data)
    saved = save_plots(plots, args.output, formats=args.formats, dpi=args.dpi)
    print(f"Saved {len(saved)} figure files to {args.output}")
    for figure in plots.values():
        plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
