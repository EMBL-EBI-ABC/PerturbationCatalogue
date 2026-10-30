"""Summary reports for MAVE score sets."""

from __future__ import annotations

import textwrap

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

MUTS_GROUPED = (
    "stop",
    "del",
    "P",
    "G",
    "A",
    "I",
    "L",
    "M",
    "V",
    "F",
    "W",
    "Y",
    "H",
    "N",
    "Q",
    "S",
    "T",
    "C",
    "R",
    "K",
    "D",
    "E",
)
HEATMAP_STRIP_WIDTH = 100
HEATMAP_ROW_HEIGHT_CM = 7
REPORT_WIDTH_CM = 30


def _required_tables(tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    required = (
        "meta",
        "scores",
        "singles",
        "density",
        "posterior",
        "fit_summary",
        "params",
        "loo",
    )
    missing = [key for key in required if key not in tables]
    if missing:
        raise ValueError(f"Missing required report tables: {', '.join(missing)}")
    invalid = [key for key in required if not isinstance(tables[key], pd.DataFrame)]
    if invalid:
        raise TypeError(
            f"Report tables must be pandas DataFrames: {', '.join(invalid)}"
        )
    return {key: tables[key] for key in required}


def _filter_urn(table: pd.DataFrame, target_score_urn: str) -> pd.DataFrame:
    if "score_urn" not in table.columns:
        if table.empty:
            return table.copy()
        raise ValueError("Non-empty report tables must contain a 'score_urn' column")
    return table.loc[table["score_urn"].eq(target_score_urn)].copy()


def _truthy_model(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False)
    if pd.api.types.is_numeric_dtype(values.dtype):
        return values.eq(1).fillna(False)
    return (
        values.astype("string").str.strip().str.lower().isin({"true", "t", "1", "yes"})
    )


def _top_model(table: pd.DataFrame) -> pd.DataFrame:
    if "top_model" not in table.columns:
        return table.copy()
    return table.loc[_truthy_model(table["top_model"])].copy()


def _finite_values(table: pd.DataFrame, column: str) -> pd.Series:
    if column not in table.columns:
        return pd.Series(dtype=float)
    values = pd.to_numeric(table[column], errors="coerce")
    return values.loc[np.isfinite(values)]


def _first_value(table: pd.DataFrame, column: str, default: str) -> str:
    if column not in table.columns:
        return default
    for value in table[column].tolist():
        if pd.isna(value):
            continue
        if isinstance(value, (float, np.floating)) and float(value).is_integer():
            return str(int(value))
        return str(value)
    return default


def _study_description(
    meta: pd.DataFrame, params: pd.DataFrame, fit_success: bool
) -> str:
    fields = [
        f"Target: {_first_value(meta, 'target_name', 'Unknown')}",
        f"Uniprot: {_first_value(meta, 'identifier_uniprot', 'Unknown')}",
        f"Organism: {_first_value(meta, 'organism', 'Unknown')}",
        f"N: {_first_value(meta, 'n', 'Missing')}",
        f"Positions: {_first_value(meta, 'n_positions', 'Missing')}",
        f"Multi-mutants: {_first_value(meta, 'multiple_mutations', 'Missing')}",
    ]
    if fit_success:
        median = _first_value(params, "median", "Missing")
        mad = _first_value(params, "mad", "Missing")
        fields.append(f"GMM Fit: Median = {median}, MAD = {mad}")
    else:
        fields.append("GMM Fit: Failed")
    return " | ".join(fields)


def _valid_heatmap_data(singles: pd.DataFrame) -> pd.DataFrame:
    required = {"position", "wt", "mut"}
    if not required.issubset(singles.columns):
        return pd.DataFrame(columns=[*required, "s", "w"])
    data = singles.copy()
    data["position"] = pd.to_numeric(data["position"], errors="coerce")
    data = data.dropna(subset=["position", "wt", "mut"])
    data["position"] = data["position"].astype(int)
    data["wt"] = data["wt"].astype(str)
    data["mut"] = data["mut"].astype(str)
    data = data.loc[data["mut"].isin(MUTS_GROUPED)]
    return data


def _heatmap_strip_count(singles: pd.DataFrame) -> int:
    data = _valid_heatmap_data(singles)
    if data.empty:
        return 0
    span = int(data["position"].max() - data["position"].min())
    return span // HEATMAP_STRIP_WIDTH + 1


def _two_slope_norm(values: np.ndarray, center: float) -> mcolors.TwoSlopeNorm:
    finite = values[np.isfinite(values)]
    if finite.size:
        radius = float(np.max(np.abs(finite - center)))
    else:
        radius = 1.0
    if radius == 0:
        radius = 1.0
    return mcolors.TwoSlopeNorm(
        vmin=center - radius,
        vcenter=center,
        vmax=center + radius,
    )


def _heatmap_colors(singles: pd.DataFrame, standardised: bool):
    column = "w" if standardised else "s"
    values = pd.to_numeric(singles.get(column, pd.Series(dtype=float)), errors="coerce")
    values_array = values.to_numpy(dtype=float, na_value=np.nan)
    if standardised:
        cmap = mcolors.LinearSegmentedColormap.from_list(
            "mave_standardised", ["#b2182b", "#ffffff", "#2166ac"]
        )
        norm = _two_slope_norm(values_array, 1.0)
        title = "Standardised Scores"
    else:
        cmap = mcolors.LinearSegmentedColormap.from_list(
            "mave_experimental", ["#762a83", "#ffffff", "#1b7837"]
        )
        finite = values_array[np.isfinite(values_array)]
        midpoint = float(np.median(finite)) if finite.size else 0.0
        norm = _two_slope_norm(values_array, midpoint)
        title = "Experimental Scores"
    cmap.set_bad("#e5e5e5")
    return column, cmap, norm, title


def _draw_empty(ax: Axes, text: str) -> None:
    ax.set_axis_off()
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes)


def _draw_heatmap(
    fig: Figure,
    subspec,
    singles: pd.DataFrame,
    standardised: bool,
) -> None:
    data = _valid_heatmap_data(singles)
    rows = _heatmap_strip_count(data)
    if rows == 0:
        ax = fig.add_subplot(subspec)
        _draw_empty(ax, "No single substitution scores")
        return

    value_column, cmap, norm, title = _heatmap_colors(data, standardised)
    if value_column not in data.columns:
        data[value_column] = np.nan
    values = pd.to_numeric(data[value_column], errors="coerce")
    data[value_column] = values.where(np.isfinite(values), np.nan)
    first_position = int(data["position"].min())
    strip_spec = subspec.subgridspec(
        rows + 1,
        1,
        height_ratios=[1.0] * rows + [0.12],
        hspace=0.18,
    )

    for strip in range(rows):
        ax = fig.add_subplot(strip_spec[strip])
        strip_start = first_position + strip * HEATMAP_STRIP_WIDTH
        strip_end = strip_start + HEATMAP_STRIP_WIDTH
        strip_data = data.loc[
            data["position"].ge(strip_start) & data["position"].lt(strip_end)
        ]
        for row_index, mutation in enumerate(MUTS_GROUPED):
            mutation_data = strip_data.loc[strip_data["mut"].eq(mutation)]
            for _, variant in mutation_data.iterrows():
                score = variant[value_column]
                face = cmap(norm(float(score))) if pd.notna(score) else "#e5e5e5"
                ax.add_patch(
                    Rectangle(
                        (float(variant["position"]) - 0.5, row_index - 0.5),
                        1,
                        1,
                        facecolor=face,
                        edgecolor="none",
                    )
                )
        # Match the source report's WT reference: outline the WT cell only where
        # the corresponding score is present.
        wt_rows = strip_data.loc[strip_data[value_column].notna()]
        for _, variant in wt_rows.drop_duplicates(["position", "wt"]).iterrows():
            if variant["wt"] not in MUTS_GROUPED:
                continue
            wt_index = MUTS_GROUPED.index(variant["wt"])
            ax.add_patch(
                Rectangle(
                    (float(variant["position"]) - 0.5, wt_index - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor="black",
                    linewidth=0.55,
                )
            )

        strip_xmax = (
            strip_end - 0.5
            if rows > 1
            else min(strip_end, int(data["position"].max()) + 1) - 0.5
        )
        ax.set_xlim(strip_start - 0.5, strip_xmax)
        ax.set_ylim(-0.5, len(MUTS_GROUPED) - 0.5)
        ax.set_yticks(range(len(MUTS_GROUPED)))
        ax.set_yticklabels(MUTS_GROUPED, fontsize=7)
        tick_positions = (
            pd.Series(range(strip_start, strip_end))
            if rows > 1
            else strip_data["position"]
        )
        ax.set_xticks(_position_ticks(tick_positions))
        ax.tick_params(axis="x", labelsize=7)
        ax.tick_params(axis="y", length=0)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if strip == 0:
            ax.set_title(title, fontsize=10, pad=3)
        ax.set_xlabel("")

    # Reserve a thin row within the allocated panel for its color key.
    cax = fig.add_subplot(strip_spec[rows])
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(scalar, cax=cax, orientation="horizontal")
    colorbar.ax.tick_params(labelsize=6, length=2)
    colorbar.set_label("w" if standardised else "s", fontsize=7, labelpad=1)


def _position_ticks(positions: pd.Series) -> list[int]:
    if positions.empty:
        return []
    low, high = int(positions.min()), int(positions.max())
    count = min(6, high - low + 1)
    return sorted({int(round(value)) for value in np.linspace(low, high, count)})


def _correlation_matrix(
    scores: pd.DataFrame,
    meta: pd.DataFrame,
    target_score_urn: str,
    variable: str,
) -> tuple[pd.DataFrame | None, str | None]:
    target_meta = meta.loc[
        meta.get("score_urn", pd.Series(dtype=object)).eq(target_score_urn)
    ]
    if target_meta.empty:
        return None, "No other studies on this gene"

    target_names = set(target_meta.get("target_name", pd.Series(dtype=object)).dropna())
    target_uniprots = set(
        target_meta.get("identifier_uniprot", pd.Series(dtype=object)).dropna()
    )
    related = pd.Series(False, index=meta.index)
    if "target_name" in meta.columns and target_names:
        related |= meta["target_name"].isin(target_names)
    if "identifier_uniprot" in meta.columns and target_uniprots:
        related |= meta["identifier_uniprot"].isin(target_uniprots)
    urns = (
        set(meta.loc[related, "score_urn"].dropna()) if "score_urn" in meta else set()
    )
    related_scores = scores.loc[
        scores.get("score_urn", pd.Series(dtype=object)).isin(urns)
    ]
    if related_scores.get("score_urn", pd.Series(dtype=object)).nunique() < 2:
        return None, "No other studies on this gene"
    if not {"score_urn", "hgvs_pro", variable}.issubset(related_scores.columns):
        return None, "No studies share variants"

    observations = related_scores[["score_urn", "hgvs_pro", variable]].copy()
    observations[variable] = pd.to_numeric(observations[variable], errors="coerce")
    observations[variable] = observations[variable].where(
        np.isfinite(observations[variable]), np.nan
    )
    observations = observations.dropna(subset=["hgvs_pro"])
    observations["score_urn"] = (
        observations["score_urn"]
        .astype(str)
        .str.replace("urn:mavedb:", "", n=1, regex=False)
    )
    averages = observations.groupby(
        ["score_urn", "hgvs_pro"], as_index=False, sort=False
    )[variable].mean()
    wide = averages.pivot(index="hgvs_pro", columns="score_urn", values=variable)
    correlations = wide.corr(method="pearson", min_periods=2)
    if correlations.isna().to_numpy().all():
        return None, "No studies share variants"
    return correlations, None


def _correlation_distances(correlations: pd.DataFrame) -> np.ndarray:
    values = correlations.to_numpy(dtype=float)
    count = values.shape[0]
    distances = np.zeros((count, count), dtype=float)
    dimensions = values.shape[1]
    for left in range(count):
        for right in range(left + 1, count):
            complete = np.isfinite(values[left]) & np.isfinite(values[right])
            if not complete.any():
                distance = 5.0
            else:
                differences = values[left, complete] - values[right, complete]
                distance = float(
                    np.sqrt(np.square(differences).sum() * dimensions / complete.sum())
                )
            distances[left, right] = distances[right, left] = distance
    return distances


def _cluster_order(correlations: pd.DataFrame) -> np.ndarray:
    count = correlations.shape[0]
    if count < 2:
        return np.arange(count)
    distances = _correlation_distances(correlations)
    tree = linkage(squareform(distances, checks=False), method="complete")
    return leaves_list(tree)


def _draw_correlation(
    fig: Figure,
    subspec,
    scores: pd.DataFrame,
    meta: pd.DataFrame,
    target_score_urn: str,
    variable: str,
) -> None:
    ax = fig.add_subplot(subspec)
    correlations, empty_text = _correlation_matrix(
        scores, meta, target_score_urn, variable
    )
    if correlations is None:
        _draw_empty(ax, empty_text or "No studies share variants")
        return

    order = _cluster_order(correlations)
    labels = correlations.columns.to_list()
    ordered = correlations.to_numpy(dtype=float)[np.ix_(order, order)]
    ordered_labels = [labels[index] for index in order]
    cmap = plt.get_cmap("RdBu").copy()
    cmap.set_bad("#c994c7")
    image = ax.imshow(
        ordered,
        vmin=-1,
        vmax=1,
        cmap=cmap,
        interpolation="nearest",
        aspect="equal",
    )
    ticks = np.arange(len(ordered_labels))
    if len(ordered_labels) <= 15:
        axis_labels = ordered_labels
    else:
        target_label = target_score_urn.replace("urn:mavedb:", "", 1)
        axis_labels = [
            label if label == target_label else "" for label in ordered_labels
        ]
    ax.set_xticks(ticks, labels=axis_labels, rotation=90, fontsize=7)
    ax.set_yticks(ticks, labels=axis_labels, fontsize=7)
    ax.tick_params(length=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    if len(ordered_labels) <= 15:
        for row in range(len(ordered_labels)):
            for column in range(len(ordered_labels)):
                value = ordered[row, column]
                label = f"{value:.2g}" if np.isfinite(value) else "NA"
                ax.text(column, row, label, ha="center", va="center", fontsize=7)
    colorbar_axes = ax.inset_axes([1.04, 0, 0.035, 1])
    colorbar = fig.colorbar(image, cax=colorbar_axes, orientation="vertical")
    colorbar.set_label("R", fontsize=8)
    colorbar.ax.tick_params(labelsize=7, length=2)


def _draw_failed_distribution(
    fig: Figure, subspec, scores: pd.DataFrame, singles: pd.DataFrame
) -> None:
    ax = fig.add_subplot(subspec)
    values = _finite_values(scores, "s").to_numpy(dtype=float)
    if not values.size:
        _draw_empty(ax, "No finite experimental scores")
        ax.set_axis_on()
        ax.set_xlabel("Experimental Score (s)")
        ax.set_ylabel("Count")
        return
    single_values = _finite_values(singles, "s").to_numpy(dtype=float)
    midpoint = (
        float(np.median(single_values))
        if single_values.size
        else float(np.median(values))
    )
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "mave_failed_distribution", ["#762a83", "#ffffff", "#1b7837"]
    )
    norm = _two_slope_norm(values, midpoint)
    counts, edges = np.histogram(values, bins=100)
    centers = (edges[:-1] + edges[1:]) / 2
    for count, left, right, center in zip(counts, edges[:-1], edges[1:], centers):
        ax.bar(left, count, width=right - left, align="edge", color=cmap(norm(center)))
    ax.set_xlabel("Experimental Score (s)")
    ax.set_ylabel("Count")
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.set_axisbelow(True)


def make_summary_plot(target_score_urn: str, tables: dict[str, pd.DataFrame]) -> Figure:
    """Compose a single-page, 30 cm-wide MAVE summary report."""
    source = _required_tables(tables)
    filtered = {
        key: _filter_urn(table, target_score_urn) for key, table in source.items()
    }
    fit_summary = _top_model(filtered["fit_summary"])
    params = _top_model(filtered["params"])
    fit_success = bool(_finite_values(filtered["scores"], "w").size)
    description = _study_description(filtered["meta"], params, fit_success)

    strip_count = _heatmap_strip_count(filtered["singles"])
    heatmap_height = HEATMAP_ROW_HEIGHT_CM * strip_count if strip_count else 2.0
    if fit_success:
        loo_height = max(3.0, 2.0 * len(filtered["loo"]))
        heights = [
            heatmap_height,
            heatmap_height,
            15.0,
            15.0,
            loo_height,
            10.0,
        ]
    else:
        heights = [heatmap_height, 5.0, 10.0]
    total_height_cm = sum(heights)
    subtitle = f"{description} | GMM Fit Failed" if not fit_success else description
    subtitle_lines = textwrap.fill(subtitle, width=145)
    subtitle_line_count = subtitle_lines.count("\n") + 1
    height_points = total_height_cm / 2.54 * 72
    title_top = 1.0 - 4.0 / height_points
    subtitle_top = title_top - (14.0 + 4.0) / height_points
    plot_top = subtitle_top - (subtitle_line_count * 11.0 + 22.0) / height_points

    fig = plt.figure(
        figsize=(REPORT_WIDTH_CM / 2.54, total_height_cm / 2.54),
        facecolor="white",
    )
    fig.text(0.5, title_top, target_score_urn, ha="center", va="top", fontsize=14)
    fig.text(0.5, subtitle_top, subtitle_lines, ha="center", va="top", fontsize=9)
    grid = fig.add_gridspec(
        len(heights),
        1,
        height_ratios=heights,
        top=plot_top,
        bottom=0.025,
        left=0.09,
        right=0.95,
        hspace=0.34,
    )

    panel_index = 0
    _draw_heatmap(fig, grid[panel_index], filtered["singles"], standardised=False)
    panel_index += 1
    if fit_success:
        _draw_heatmap(fig, grid[panel_index], filtered["singles"], standardised=True)
        panel_index += 1

    if fit_success:
        from .diagnostic_plots import draw_loo, draw_posterior, draw_transform

        draw_transform(fig, grid[panel_index], filtered["scores"], filtered["density"])
        panel_index += 1
        draw_posterior(
            fig,
            grid[panel_index],
            filtered["posterior"],
            fit_summary,
            params,
        )
        panel_index += 1
        loo_ax = fig.add_subplot(grid[panel_index])
        draw_loo(loo_ax, filtered["loo"])
        panel_index += 1
        _draw_correlation(
            fig,
            grid[panel_index],
            source["scores"],
            source["meta"],
            target_score_urn,
            "w",
        )
    else:
        _draw_failed_distribution(
            fig,
            grid[panel_index],
            filtered["scores"],
            filtered["singles"],
        )
        panel_index += 1
        _draw_correlation(
            fig,
            grid[panel_index],
            source["scores"],
            source["meta"],
            target_score_urn,
            "s",
        )
    return fig
