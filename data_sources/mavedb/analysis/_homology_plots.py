"""Matplotlib plots for the MAVE homology analysis."""

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
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

PLOT_NAMES = (
    "overlapping_variants",
    "pearson_distribution",
    "spearman_distribution",
    "pearson_bars",
    "correlation_vs_percent_identity",
    "correlation_vs_percent_identity_filtered",
    "correlation_vs_percent_identity_boxplot",
    "clinvar_roc_curves",
    "clinvar_roc_reduced",
)

RELATIONSHIP_ORDER = (
    "Same Experiment - Identical",
    "Same Study - Identical",
    "Same Study - Paralog (same species)",
    "Different Study - Identical",
    "Different Study - Ortholog (one to one)",
    "Different Study - Paralog (same species)",
    "Different Study - Paralog (diff. species)",
)
GENE_RELATIONSHIP_COLORS = {
    "Identical": "#80b1d3",
    "Ortholog (one to one)": "#b3de69",
    "Paralog (same species)": "#fdb462",
    "Paralog (diff. species)": "#fb9a99",
}
RELATIONSHIP_COLORS = {
    "Same Experiment - Identical": "#1f78b4",
    "Same Study - Identical": "#a6cee3",
    "Same Study - Paralog (same species)": "#e31a1c",
    "Different Study - Identical": "#80b1d3",
    "Different Study - Ortholog (one to one)": "#b3de69",
    "Different Study - Paralog (same species)": "#fdb462",
    "Different Study - Paralog (diff. species)": "#fb9a99",
}
SIGNIFICANCE_COLORS = {
    "q > 0.01": "#e41a1c",
    "q <= 0.01": "#377eb8",
    "q unavailable": "#999999",
}
CLINVAR_SCORE_LABELS = {
    "w": "Standardised Score (w)",
    "aw": "Abs. Std. Score (|w - 1|)",
    "p": "Neutral SD Scaled w",
    "ap": "Abs. SD scaled w",
}
CLINVAR_SCORE_ORDER = ("w", "aw", "p", "ap")
CLINVAR_TYPE_ORDER = (
    "All",
    "Ortholog (one to one)",
    "Ortholog (one to many)",
    "Ortholog (many to many)",
    "Paralog (diff. species)",
    "Paralog (same species)",
    "Alt Allele",
    "Identical",
)
SAME_WT_ORDER = ("Both", "Same WT", "Diff. WT")
CLINVAR_TYPE_COLORS = {
    "All": "#e41a1c",
    "Ortholog (one to one)": "#377eb8",
    "Paralog (diff. species)": "#4daf4a",
    "Paralog (same species)": "#984ea3",
}
SAME_WT_LINESTYLES = {"Both": "-", "Diff. WT": "--", "Same WT": ":"}


def _frame(data: pd.DataFrame | None) -> pd.DataFrame:
    return data if isinstance(data, pd.DataFrame) else pd.DataFrame()


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


def _observed(values: Sequence[Any]) -> list[str]:
    result: list[str] = []
    for value in values:
        if value is None or pd.isna(value):
            continue
        label = str(value)
        if label and label not in result:
            result.append(label)
    return result


def _ordered(values: Sequence[Any], preferred: Sequence[str]) -> list[str]:
    present = set(_observed(values))
    ordered = [label for label in preferred if label in present]
    ordered.extend(sorted(present - set(ordered)))
    return ordered


def _empty_axis(ax: Axes, message: str) -> None:
    ax.set_axis_off()
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=9,
    )


def _relationship_values(data: pd.DataFrame) -> np.ndarray:
    return _text(data, "relationship")


def _gene_relationship_values(data: pd.DataFrame) -> np.ndarray:
    if "gene_relationship" in data:
        return _text(data, "gene_relationship")
    relationship = _relationship_values(data)
    return np.asarray(
        [value.split(" - ", 1)[-1] if value else "" for value in relationship],
        dtype=object,
    )


def _significance_values(data: pd.DataFrame) -> np.ndarray:
    if "pearson_q" not in data:
        return np.full(len(data), "q unavailable", dtype=object)
    raw = data["pearson_q"]
    numeric = pd.to_numeric(raw, errors="coerce").to_numpy(dtype=float, na_value=np.nan)
    labels = np.full(len(data), "q unavailable", dtype=object)
    for index, value in enumerate(raw.tolist()):
        if isinstance(value, (bool, np.bool_)):
            labels[index] = "q <= 0.01" if value else "q > 0.01"
        elif np.isfinite(numeric[index]):
            labels[index] = "q <= 0.01" if numeric[index] <= 0.01 else "q > 0.01"
    return labels


def _draw_empty_figure(title: str, message: str) -> Figure:
    fig, ax = plt.subplots(figsize=(7, 3.5), constrained_layout=True)
    ax.set_title(title)
    _empty_axis(ax, message)
    return fig


def _overlapping_variants(data: pd.DataFrame) -> Figure:
    gene_relationship = _gene_relationship_values(data)
    relationships = _ordered(
        gene_relationship,
        tuple(GENE_RELATIONSHIP_COLORS),
    )
    if not relationships:
        return _draw_empty_figure("Overlapping variants", "No gene relationship data")

    fig, axes = plt.subplots(
        len(relationships),
        1,
        figsize=(8.5, max(2.2, 1.9 * len(relationships))),
        squeeze=False,
        sharex=False,
        constrained_layout=True,
    )
    shared = _numeric(data, "shared_variants")
    for index, relationship in enumerate(relationships):
        ax = axes[index, 0]
        values = shared[(gene_relationship == relationship) & np.isfinite(shared)]
        values = values[values > 0]
        ax.set_title(relationship, loc="left", fontsize=9)
        if values.size:
            lower = float(values.min())
            upper = float(values.max())
            if lower == upper:
                lower = max(lower / 2, np.finfo(float).tiny)
                upper *= 2
            edges = np.geomspace(lower, upper, 51)
            ax.hist(
                values,
                bins=edges,
                color=GENE_RELATIONSHIP_COLORS.get(relationship, "#999999"),
                edgecolor="white",
                linewidth=0.35,
            )
            ax.set_xscale("log")
            ax.axvline(50, color="black", linewidth=0.8)
            ax.grid(axis="y", color="#dddddd", linewidth=0.5)
            ax.set_axisbelow(True)
        else:
            ax.text(
                0.5,
                0.5,
                "No positive shared-variant counts",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=8,
            )
        if index == len(relationships) - 1:
            ax.set_xlabel("Number of shared variants")
        ax.set_ylabel("Study pairs")
    return fig


def _correlation_distribution(data: pd.DataFrame, variable: str) -> Figure:
    values = _numeric(data, variable)
    relationships = _relationship_values(data)
    valid_labels = [
        relationships[index]
        for index in range(len(data))
        if relationships[index] and np.isfinite(values[index])
    ]
    groups = _ordered(valid_labels, RELATIONSHIP_ORDER)
    title = "Pearson correlation" if variable == "pearson" else "Spearman correlation"
    if not groups:
        return _draw_empty_figure(title, "No finite correlation values")

    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.72 * len(groups) + 1.4)))
    y_positions = np.arange(len(groups), dtype=float)
    for position, relationship in zip(y_positions, groups):
        selected = (relationships == relationship) & np.isfinite(values)
        group_values = values[selected]
        gene_relationship = relationship.split(" - ", 1)[-1]
        color = GENE_RELATIONSHIP_COLORS.get(gene_relationship, "#999999")
        if group_values.size > 1 and np.ptp(group_values) > 0:
            violin = ax.violinplot(
                [group_values],
                positions=[position],
                vert=False,
                widths=0.72,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )
            for body in violin["bodies"]:
                body.set_facecolor(color)
                body.set_edgecolor(color)
                body.set_alpha(0.5)
        if group_values.size > 1:
            box = ax.boxplot(
                [group_values],
                positions=[position],
                vert=False,
                widths=0.16,
                patch_artist=True,
                manage_ticks=False,
                showfliers=True,
                flierprops={"markersize": 2.5, "markeredgewidth": 0.5},
            )
            for patch in box["boxes"]:
                patch.set_facecolor(color)
                patch.set_edgecolor(color)
            for median in box["medians"]:
                median.set_color("black")
                median.set_linewidth(0.9)
        elif group_values.size == 1:
            ax.scatter(group_values, [position], color=color, s=18, zorder=4)
        ax.text(
            1.015,
            position,
            f"n = {group_values.size}",
            transform=blended_transform_factory(ax.transAxes, ax.transData),
            ha="left",
            va="center",
            fontsize=8,
            clip_on=False,
        )

    ax.axvline(0, color="lightgrey", linewidth=0.8, zorder=0)
    ax.set_xlim(-1.05, 1.05)
    ax.set_yticks(y_positions, labels=groups)
    ax.set_xlabel(
        "Pearson Correlation" if variable == "pearson" else "Spearman Correlation"
    )
    ax.set_ylabel("")
    ax.set_title(title)
    ax.grid(axis="x", linestyle=":", color="lightgrey", linewidth=0.6)
    ax.tick_params(axis="y", length=0)
    ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.36, right=0.86, top=0.90, bottom=0.12)
    return fig


def _pearson_bars(data: pd.DataFrame) -> Figure:
    values = _numeric(data, "pearson")
    valid = np.isfinite(values)
    if not valid.any():
        return _draw_empty_figure(
            "Pearson correlations", "No finite Pearson correlations"
        )

    frame = data.loc[valid].copy()
    values = values[valid]
    frame["_pearson"] = values
    frame = frame.sort_values("_pearson", ascending=False, kind="stable")
    values = frame["_pearson"].to_numpy(dtype=float)
    relationship = _relationship_values(frame)
    significance = _significance_values(frame)
    colors = [
        to_rgba(
            RELATIONSHIP_COLORS.get(label, "#999999"),
            (
                1.0
                if significant == "q <= 0.01"
                else 0.05 if significant == "q > 0.01" else 0.4
            ),
        )
        for label, significant in zip(relationship, significance)
    ]

    fig, ax = plt.subplots(figsize=(11, 5.5))
    positions = np.arange(len(values))
    ax.bar(positions, values, color=colors, width=0.9, linewidth=0)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlim(-0.7, max(len(values) - 0.3, 0.7))
    ax.set_xticks([])
    ax.set_xlabel("Study pairs, sorted by Pearson correlation")
    ax.set_ylabel("Pearson Correlation")
    ax.set_title("Pearson correlations across MAVE study pairs")

    present_relationships = _ordered(relationship, tuple(RELATIONSHIP_COLORS))
    relationship_handles = [
        Patch(
            facecolor=RELATIONSHIP_COLORS.get(label, "#999999"),
            label=label,
        )
        for label in present_relationships
    ]
    present_significance = _ordered(
        significance,
        ("q > 0.01", "q <= 0.01", "q unavailable"),
    )
    significance_handles = [
        Patch(
            facecolor="grey",
            alpha=1.0 if label == "q <= 0.01" else 0.05 if label == "q > 0.01" else 0.4,
            label=label,
        )
        for label in present_significance
    ]
    handles = relationship_handles + significance_handles
    if handles:
        ax.legend(
            handles=handles,
            title="Relationship / FDR corrected Pearson p-value",
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            frameon=False,
            fontsize=7,
        )
        fig.subplots_adjust(left=0.08, right=0.75, top=0.92, bottom=0.12)
    else:
        fig.subplots_adjust(left=0.08, right=0.96, top=0.92, bottom=0.12)
    return fig


def _regression(
    ax: Axes, x: np.ndarray, y: np.ndarray, color: str, label: str
) -> str | None:
    valid = np.isfinite(x) & np.isfinite(y)
    x_values = x[valid]
    y_values = y[valid]
    if x_values.size < 2 or np.unique(x_values).size < 2:
        return
    slope, intercept = np.polyfit(x_values, y_values, 1)
    x_line = np.linspace(float(x_values.min()), float(x_values.max()), 100)
    line = ax.plot(
        x_line,
        slope * x_line + intercept,
        color=color,
        linewidth=1.1,
        zorder=2,
    )[0]
    line.set_gid(f"regression:{label}")
    predicted = slope * x_values + intercept
    total = float(np.square(y_values - y_values.mean()).sum())
    residual = float(np.square(y_values - predicted).sum())
    r_squared = 1.0 - residual / total if total > 0 else np.nan
    if np.isfinite(r_squared):
        return f"{label}: R² = {r_squared:.2f}"
    return None


def _correlation_identity(data: pd.DataFrame, *, filtered: bool = False) -> Figure:
    title = (
        "Correlation vs. homolog identity (filtered)"
        if filtered
        else "Correlation vs. homolog identity"
    )
    x = _numeric(data, "percent_identity")
    y = _numeric(data, "pearson")
    shared = _numeric(data, "shared_variants")
    keep = np.isfinite(x) & np.isfinite(y) & (x < 100)
    if filtered:
        keep &= ((x < 85) | (x > 90)) & (shared > 50)
    if not keep.any():
        return _draw_empty_figure(
            title, "No homolog correlations meet the plot filters"
        )

    frame = data.loc[keep].copy()
    x = x[keep]
    y = y[keep]
    significance = _significance_values(frame)
    low = _numeric(frame, "pearson_low")
    high = _numeric(frame, "pearson_high")
    fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)

    err_valid = np.isfinite(low) & np.isfinite(high)
    if err_valid.any():
        lower = np.maximum(0, y[err_valid] - low[err_valid])
        upper = np.maximum(0, high[err_valid] - y[err_valid])
        ax.errorbar(
            x[err_valid],
            y[err_valid],
            yerr=np.vstack([lower, upper]),
            fmt="none",
            ecolor="black",
            linewidth=0.55,
            capsize=0,
            zorder=1,
        )

    fit_labels: list[str] = []
    for label in ("q > 0.01", "q <= 0.01", "q unavailable"):
        selected = significance == label
        if not selected.any():
            continue
        color = SIGNIFICANCE_COLORS[label]
        points = ax.scatter(
            x[selected],
            y[selected],
            color=color,
            s=18,
            marker="o",
            linewidths=0,
            zorder=3,
        )
        points.set_gid(
            "points:correlation_vs_percent_identity_filtered"
            if filtered
            else "points:correlation_vs_percent_identity"
        )
        fit_label = _regression(ax, x[selected], y[selected], color, label)
        if fit_label:
            fit_labels.append(fit_label)

    ax.set_xlim(0, 100)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel("Homolog % identity")
    ax.set_ylabel("Pearson Correlation")
    ax.set_title(title)
    ax.grid(color="#dddddd", linewidth=0.5)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color=color,
            markersize=5,
            label=label,
        )
        for label, color in SIGNIFICANCE_COLORS.items()
        if label in significance
    ]
    if handles:
        ax.legend(
            handles=handles,
            title="FDR corrected Pearson p-value",
            loc="lower right",
            frameon=False,
            fontsize=8,
        )
    if fit_labels:
        ax.text(
            0.02,
            0.98,
            "\n".join(fit_labels),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
        )
    if filtered:
        ax.text(
            0.5,
            -0.17,
            "85–90% identity and pairs with 50 or fewer shared variants excluded",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8,
        )
    return fig


def _identity_boxplot(data: pd.DataFrame) -> Figure:
    identity = _numeric(data, "percent_identity")
    correlation = _numeric(data, "pearson")
    keep = np.isfinite(identity) & np.isfinite(correlation) & (identity < 100)
    keep &= identity >= 0
    if not keep.any():
        return _draw_empty_figure(
            "Pearson correlation by homolog identity",
            "No correlations with homolog identity below 100%",
        )

    frame = data.loc[keep].copy()
    identity = identity[keep]
    correlation = correlation[keep]
    significance = _significance_values(frame)
    edges = np.arange(0, 105, 5, dtype=float)
    labels = [f"{int(left)}–{int(right)}" for left, right in zip(edges[:-1], edges[1:])]
    bins = pd.cut(identity, bins=edges, labels=labels, right=True, include_lowest=True)

    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)
    offsets = {"q > 0.01": -0.16, "q <= 0.01": 0.16, "q unavailable": 0.0}
    widths = {"q unavailable": 0.22, "q > 0.01": 0.25, "q <= 0.01": 0.25}
    for significance_label in ("q > 0.01", "q <= 0.01", "q unavailable"):
        for bin_index, bin_label in enumerate(labels):
            selected = (bins == bin_label) & (significance == significance_label)
            values = correlation[selected]
            if values.size == 0:
                continue
            position = bin_index + offsets[significance_label]
            color = SIGNIFICANCE_COLORS[significance_label]
            if values.size == 1:
                points = ax.scatter([position], values, color=color, s=12, zorder=3)
                points.set_gid(f"identity-box:{bin_label}:{significance_label}")
                continue
            box = ax.boxplot(
                [values],
                positions=[position],
                widths=widths[significance_label],
                patch_artist=True,
                manage_ticks=False,
                showfliers=True,
                flierprops={"markersize": 2.0, "markeredgewidth": 0.4},
            )
            for patch in box["boxes"]:
                patch.set_facecolor(color)
                patch.set_edgecolor(color)
                patch.set_alpha(0.8)
                patch.set_gid(f"identity-box:{bin_label}:{significance_label}")
            for median in box["medians"]:
                median.set_color("black")
                median.set_linewidth(0.8)

    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xticks(np.arange(len(labels)), labels=labels, rotation=55, ha="right")
    ax.set_xlabel("Homolog % identity")
    ax.set_ylabel("Pearson Correlation")
    ax.set_title("Pearson correlation by homolog identity")
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    handles = [
        Patch(facecolor=color, alpha=0.8, label=label)
        for label, color in SIGNIFICANCE_COLORS.items()
        if label in significance
    ]
    if handles:
        ax.legend(
            handles=handles,
            title="FDR corrected Pearson p-value",
            loc="upper right",
            frameon=False,
            fontsize=8,
        )
    return fig


def _facet_pairs(data: pd.DataFrame) -> list[tuple[str, str]]:
    type_values = _text(data, "type")
    same_wt_values = _text(data, "same_wt")
    pairs = {
        (str(type_value), str(same_wt))
        for type_value, same_wt in zip(type_values, same_wt_values)
        if type_value and same_wt
    }
    type_rank = {label: index for index, label in enumerate(CLINVAR_TYPE_ORDER)}
    wt_rank = {label: index for index, label in enumerate(SAME_WT_ORDER)}
    return sorted(
        pairs,
        key=lambda pair: (
            type_rank.get(pair[0], len(type_rank)),
            pair[0],
            wt_rank.get(pair[1], len(wt_rank)),
            pair[1],
        ),
    )


def _matching_clinvar_rows(
    data: pd.DataFrame, score_type: str, facet: tuple[str, str]
) -> pd.DataFrame:
    if not {"score_type", "type", "same_wt"}.issubset(data.columns):
        return data.iloc[0:0].copy()
    score_types = _text(data, "score_type")
    types = _text(data, "type")
    same_wt = _text(data, "same_wt")
    return data.loc[
        (score_types == score_type) & (types == facet[0]) & (same_wt == facet[1])
    ].copy()


def _curve_rows(data: pd.DataFrame, column: str) -> np.ndarray:
    values = _numeric(data, column)
    return np.isfinite(values)


def _draw_roc_or_pr(
    ax: Axes,
    rows: pd.DataFrame,
    *,
    score_type: str,
    facet: tuple[str, str],
    metric: str,
) -> None:
    x_column, y_column, auc_column = (
        ("fpr", "tpr", "auc") if metric == "roc" else ("tpr", "precision", "pr_auc")
    )
    valid = _curve_rows(rows, x_column) & _curve_rows(rows, y_column)
    rows = rows.loc[valid].copy()
    urns = _text(rows, "homolog_urn")
    if rows.empty:
        _empty_axis(ax, f"No {metric.upper()} data")
        return

    order = _ordered(urns, ("All",))
    order = [label for label in order if label != "All"]
    if "All" in urns:
        order.append("All")
    for urn in order:
        selected = urns == urn
        group = rows.loc[selected].copy()
        x = _numeric(group, x_column)
        y = _numeric(group, y_column)
        finite = np.isfinite(x) & np.isfinite(y)
        group = group.loc[finite]
        x = x[finite]
        y = y[finite]
        if not x.size:
            continue
        curve_order = np.argsort(x, kind="stable")
        x = x[curve_order]
        y = y[curve_order]
        pooled = urn == "All"
        line = ax.plot(
            x,
            y,
            color="black" if pooled else "lightgrey",
            linestyle="-" if pooled else ":",
            linewidth=1.4 if pooled else 0.75,
            zorder=3 if pooled else 2,
        )[0]
        line.set_gid(f"clinvar-{metric}:{score_type}:{facet[0]}:{facet[1]}:{urn}")
        if pooled and auc_column in group:
            auc_values = _numeric(group, auc_column)
            finite_auc = auc_values[np.isfinite(auc_values)]
            if finite_auc.size:
                label = "AUC" if metric == "roc" else "PR AUC"
                ax.text(
                    0.97,
                    0.05,
                    f"{label} = {finite_auc[0]:.3g}",
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=6.5,
                )

    if metric == "roc":
        diagonal = ax.plot(
            [0, 1], [0, 1], color="#555555", linestyle="--", linewidth=0.65, zorder=1
        )[0]
        diagonal.set_gid("clinvar-roc-reference")
        ax.set_ylabel("Sensitivity / TPR", fontsize=7)
        ax.set_xlabel("1 − Specificity / FPR", fontsize=7)
    else:
        ax.set_ylabel("Precision", fontsize=7)
        ax.set_xlabel("Recall (TPR)", fontsize=7)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(color="#dddddd", linewidth=0.45)
    ax.tick_params(labelsize=6, length=2)


def _clinvar_roc_curves(data: pd.DataFrame) -> Figure:
    score_values = _text(data, "score_type")
    score_types = _ordered(score_values, CLINVAR_SCORE_ORDER)
    facets = _facet_pairs(data)
    if not facets:
        facets = [("All", "Both")]
    if not score_types:
        score_types = ["w"]
    n_columns = len(facets)
    n_rows = 2 * len(score_types)
    width = max(8.0, min(2.35 * n_columns, 32.0))
    height = max(4.5, 1.9 * n_rows)
    fig, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(width, height),
        squeeze=False,
        sharex=False,
        sharey=False,
    )

    for score_index, score_type in enumerate(score_types):
        score_label = CLINVAR_SCORE_LABELS.get(score_type, score_type)
        roc_row = 2 * score_index
        pr_row = roc_row + 1
        for facet_index, facet in enumerate(facets):
            rows = _matching_clinvar_rows(data, score_type, facet)
            roc_ax = axes[roc_row, facet_index]
            pr_ax = axes[pr_row, facet_index]
            _draw_roc_or_pr(
                roc_ax,
                rows,
                score_type=score_type,
                facet=facet,
                metric="roc",
            )
            _draw_roc_or_pr(
                pr_ax,
                rows,
                score_type=score_type,
                facet=facet,
                metric="pr",
            )
            if score_index == 0:
                roc_ax.set_title(f"{facet[0]}\n{facet[1]}", fontsize=7)
            if facet_index == 0:
                roc_ax.set_ylabel(f"{score_label}\nSensitivity / TPR", fontsize=7)
                pr_ax.set_ylabel(f"{score_label}\nPrecision", fontsize=7)

    legend_handles = [
        Line2D(
            [0], [0], color="lightgrey", linestyle=":", linewidth=1, label="Per Study"
        ),
        Line2D([0], [0], color="black", linestyle="-", linewidth=1.4, label="Overall"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=False,
        fontsize=8,
        ncol=2,
    )
    fig.suptitle("ClinVar pathogenicity prediction by homolog relationship", y=0.995)
    fig.subplots_adjust(
        left=0.07,
        right=0.99,
        top=0.94,
        bottom=0.035,
        wspace=0.35,
        hspace=0.42,
    )
    return fig


def _clinvar_roc_reduced(data: pd.DataFrame) -> Figure:
    if not {"score_type", "homolog_urn"}.issubset(data.columns):
        rows = data.iloc[0:0].copy()
    else:
        score_types = _text(data, "score_type")
        urns = _text(data, "homolog_urn")
        rows = data.loc[(score_types == "p") & (urns == "All")].copy()
    valid = _curve_rows(rows, "fpr") & _curve_rows(rows, "tpr")
    rows = rows.loc[valid].copy()
    if rows.empty:
        fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        _empty_axis(ax, "No pooled neutral-SD ClinVar ROC data (score_type = 'p')")
        ax.set_axis_on()
        ax.set_xlabel("1 − Specificity / FPR")
        ax.set_ylabel("Sensitivity / TPR")
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_title("Pooled ClinVar ROC curves")
        return fig

    types = _ordered(_text(rows, "type"), CLINVAR_TYPE_ORDER)
    same_wt_values = _text(rows, "same_wt")
    pairs = _facet_pairs(rows)
    fig, ax = plt.subplots(
        figsize=(8.5, max(5.5, 0.28 * len(pairs) + 3.5)),
        constrained_layout=True,
    )
    for type_label in types:
        for same_wt in SAME_WT_ORDER:
            selected = (_text(rows, "type") == type_label) & (same_wt_values == same_wt)
            group = rows.loc[selected].copy()
            if group.empty:
                continue
            x = _numeric(group, "fpr")
            y = _numeric(group, "tpr")
            order = np.argsort(x, kind="stable")
            line = ax.plot(
                x[order],
                y[order],
                color=CLINVAR_TYPE_COLORS.get(type_label, "#999999"),
                linestyle=SAME_WT_LINESTYLES.get(same_wt, "-"),
                linewidth=1.2,
                label=f"{type_label}, {same_wt}",
            )[0]
            line.set_gid(f"clinvar-reduced:{type_label}:{same_wt}")

    diagonal = ax.plot([0, 1], [0, 1], color="#555555", linestyle="--", linewidth=0.75)[
        0
    ]
    diagonal.set_gid("clinvar-reduced-reference")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("1 − Specificity / FPR")
    ax.set_ylabel("Sensitivity / TPR")
    ax.set_title("Pooled neutral-SD ClinVar ROC curves")
    ax.grid(color="#dddddd", linewidth=0.55)
    ax.legend(
        title="Relationship, WT match",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=7,
    )
    return fig


def build_plots(
    correlations: pd.DataFrame, clinvar_roc: pd.DataFrame
) -> dict[str, Figure]:
    """Build the nine figures produced by the R MAVE homology analysis."""
    correlations = _frame(correlations)
    clinvar_roc = _frame(clinvar_roc)
    return {
        "overlapping_variants": _overlapping_variants(correlations),
        "pearson_distribution": _correlation_distribution(correlations, "pearson"),
        "spearman_distribution": _correlation_distribution(correlations, "spearman"),
        "pearson_bars": _pearson_bars(correlations),
        "correlation_vs_percent_identity": _correlation_identity(correlations),
        "correlation_vs_percent_identity_filtered": _correlation_identity(
            correlations, filtered=True
        ),
        "correlation_vs_percent_identity_boxplot": _identity_boxplot(correlations),
        "clinvar_roc_curves": _clinvar_roc_curves(clinvar_roc),
        "clinvar_roc_reduced": _clinvar_roc_reduced(clinvar_roc),
    }


def save_plots(
    plots: dict[str, Figure],
    output: Path,
    formats: Sequence[str] = ("png", "pdf"),
    dpi: int = 150,
) -> list[Path]:
    """Save each named figure in every requested format and return the paths."""
    if dpi <= 0:
        raise ValueError("dpi must be a positive integer")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for name, figure in plots.items():
        for requested_format in formats:
            file_format = str(requested_format).strip().lower().lstrip(".")
            if not file_format:
                raise ValueError("formats cannot contain an empty file extension")
            path = output / f"{name}.{file_format}"
            figure.savefig(path, format=file_format, dpi=dpi, bbox_inches="tight")
            paths.append(path)
    return paths
