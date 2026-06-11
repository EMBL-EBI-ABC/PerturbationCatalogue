"""Reusable target data table component with CSS grid layout."""

from __future__ import annotations

import hashlib
from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional

import pandas as pd
import plotly.express as px
from dash import html, dcc
import dash_bootstrap_components as dbc

from utils import COLORS, format_number, format_target_label

GridControlFactory = Optional[Callable[[str, Dict[str, Any]], Any]]

GRID_STYLE = {
    "display": "grid",
    "gridTemplateColumns": "5fr 7fr",
    "columnGap": "1rem",
    "rowGap": "0.15rem",
    "gridAutoRows": "min-content",
}

HEADER_TITLES = ["Dataset", "Effect"]

DATASET_METADATA_FIELDS = [
    ("dataset_tissues", "Tissue"),
    ("dataset_cell_types", "Cell type"),
    ("dataset_cell_lines", "Cell line"),
    ("dataset_library_perturbation_types", "Library perturbation"),
    ("dataset_diseases", "Disease"),
    ("dataset_sexes", "Sex"),
    ("dataset_developmental_stages", "Developmental stage"),
    ("dataset_score_interpretation", "Score interpretation"),
    ("dataset_readout_technology_labels", "Readout technology"),
]
DATASET_FIELD_FALLBACKS = {
    "dataset_id": ["id"],
    "dataset_tissues": ["tissue_labels", "tissues", "tissue"],
    "dataset_cell_types": ["cell_type_labels", "cell_types", "cell_type"],
    "dataset_cell_lines": ["cell_line_labels", "cell_lines", "cell_line"],
    "dataset_library_perturbation_types": [
        "library_perturbation_type_labels",
        "library_perturbation_types",
        "library_perturbation_type",
        "library_type",
    ],
    "dataset_diseases": ["disease_labels", "diseases", "disease"],
    "dataset_sexes": ["sex_labels", "sexes", "sex"],
    "dataset_developmental_stages": [
        "developmental_stage_labels",
        "developmental_stages",
        "developmental_stage",
    ],
    "dataset_score_interpretation": ["score_interpretation"],
    "dataset_readout_technology_labels": [
        "readout_technology_labels",
        "readout_technology",
    ],
}

GREEN = "#2acc06"
RED = "#ff4824"
MINUS = "−"


def _format_target(perturbation: Dict[str, Any]) -> str:
    """Format a perturbation target, falling back to target ID."""
    return format_target_label(
        perturbation,
        symbol_key="target_symbol",
        ensg_key="target_ensg",
        id_key="target_id",
    )


# Color mapping for metadata field badges
METADATA_FIELD_COLORS = {
    "Tissue": "primary",
    "Cell type": "info",
    "Cell line": "secondary",
    "Library perturbation": "warning",
    "Disease": "danger",
    "Sex": "success",
    "Developmental stage": "dark",
    "Score interpretation": "secondary",
    "Readout technology": "info",
}


def TargetDataTable(
    data: Optional[List[Dict[str, Any]]],
    modality: str,
    table_id: Optional[str] = None,
    rows_per_dataset_limit: Optional[int] = None,
    empty_message: str = "No results found.",
    error_message: Optional[str] = None,
    dataset_control_factory: GridControlFactory = None,
    effect_gene_source: str = "effect",
    section_id: Optional[str] = None,
    download_url_base: Optional[str] = None,
    perturbed_gene_name: Optional[str] = None,
):
    """Render the reusable data table."""
    datasets = data or []

    grid_children: List[Any] = []
    grid_children.extend(_build_header_cells(effect_gene_source, section_id))

    if error_message:
        grid_children.append(_grid_message(error_message, tone="error"))

    if not datasets:
        grid_children.append(_grid_message(empty_message))
    else:
        for entry in datasets:
            grid_children.extend(
                _build_dataset_rows(
                    entry,
                    modality,
                    rows_per_dataset_limit,
                    dataset_control_factory,
                    effect_gene_source,
                    section_id,
                    download_url_base,
                    perturbed_gene_name,
                )
            )

    div_kwargs = {
        "children": grid_children,
        "className": "target-data-table-grid",
        "style": GRID_STYLE,
    }
    if table_id is not None:
        div_kwargs["id"] = table_id

    return html.Div(**div_kwargs)


def _build_header_cells(
    effect_gene_source: str = "effect",
    section_id: Optional[str] = None,
) -> List[html.Div]:
    """Create the table header with compact titles."""
    # For both Perturb-Seq sections, always show "Effect" header
    if section_id in ("perturb_seq_perturbed", "perturb_seq_affected"):
        header_titles = ["Dataset", "Effect"]
    else:
        header_titles = [
            "Dataset",
            "Perturbation" if effect_gene_source == "perturbation" else "Effect",
        ]

    headers = []
    for title in header_titles:
        headers.append(
            html.Div(html.H3(title, className="fw-semibold mb-1"), className="pb-1")
        )

    return headers


def _build_dataset_rows(
    entry: Dict[str, Any],
    modality: str,
    rows_per_dataset_limit: Optional[int],
    dataset_control_factory: GridControlFactory,
    effect_gene_source: str,
    section_id: Optional[str],
    download_url_base: Optional[str] = None,
    perturbed_gene_name: Optional[str] = None,
) -> List[Any]:
    dataset_meta = entry.get("dataset") or {}
    dataset_id = _resolve_meta_value(dataset_meta, "dataset_id") or "Dataset"
    results = entry.get("results") or []

    # For MAVE, Perturb-Seq, and CRISPR we show single components (heatmap/table), so row_span should be 1
    is_perturb_seq_table = modality == "perturb-seq" and section_id in (
        "perturb_seq_perturbed",
        "perturb_seq_affected",
    )
    is_crispr_table = modality == "crispr-screen"
    uses_single_component = (
        modality == "mave" or is_perturb_seq_table or is_crispr_table
    )
    row_span = 1 if uses_single_component else max(len(results), 1)

    children: List[Any] = [
        _render_dataset_cell(dataset_meta, row_span, modality),
    ]

    # Build download URL for this dataset
    download_url = None
    if download_url_base and dataset_id:
        # Append dataset_id filter to the base URL
        separator = "&" if "?" in download_url_base else "?"
        download_url = f"{download_url_base}{separator}dataset_id={dataset_id}"

    # Get dataset cell_type for fallback when effect cell_type is N/A
    ds_cell_types = _resolve_meta_value(dataset_meta, "dataset_cell_types")
    ds_cell_type = (
        ds_cell_types[0]
        if isinstance(ds_cell_types, list) and ds_cell_types
        else ds_cell_types
    )

    # Build GSEA button data for this dataset (perturb_seq_perturbed only)
    gsea_button_data = None
    if perturbed_gene_name and dataset_id and section_id == "perturb_seq_perturbed":
        gsea_button_data = {
            "dataset_id": dataset_id,
            "perturbed_gene_name": perturbed_gene_name,
            "dataset_cell_types": ds_cell_type,
        }

    if results:
        # For MAVE modality, render a single heatmap instead of individual result cells
        if modality == "mave":
            children.append(_mave_heatmap_effect(results, download_url))
        # For Perturb-Seq sections, render as a table
        elif section_id in ("perturb_seq_perturbed", "perturb_seq_affected"):
            children.append(
                _perturb_seq_table(
                    results, section_id, download_url, gsea_button_data, ds_cell_type
                )
            )
        # For CRISPR, render as a table
        elif is_crispr_table:
            children.append(_crispr_table(results, download_url))
        else:
            # Fallback for other modalities
            for result in results:
                children.append(
                    _render_result_cell(
                        result, modality, effect_gene_source, section_id
                    )
                )
    else:
        children.append(
            _grid_message(
                "No rows available for this dataset.",
                columns="2 / 3",
            )
        )

    truncated_flag = entry.get("truncated")
    if rows_per_dataset_limit and truncated_flag:
        limit_value = entry.get("rows_per_dataset_limit", rows_per_dataset_limit)
        children.append(_truncation_notice(limit_value))

    if dataset_control_factory:
        control = dataset_control_factory(dataset_id, entry)
        if control:
            children.append(
                html.Div(
                    control,
                    style={
                        "gridColumn": "2 / 3",
                        "marginTop": "0.5rem",
                        "justifySelf": "end",
                    },
                )
            )

    children.append(_dataset_separator())
    return children


def _render_dataset_cell(
    dataset_meta: Dict[str, Any], span_rows: int, modality: str = ""
):
    dataset_id = _resolve_meta_value(dataset_meta, "dataset_id")
    formatted_id = _format_dataset_id(dataset_id)
    url_dataset_id = _dataset_id_to_url_format(dataset_id)

    # Create the dataset title with [more info] link
    title_elements = [
        html.Span(formatted_id, className="h4 fw-semibold text-break"),
    ]

    # Perturb-seq provenance badge, shown just before the [more info] link.
    # The modality search maps ES fields to their api_name, so this is the api_name
    # (the /dataset and /search endpoints instead return the raw es_field).
    provenance_badge = reprocessed_badge(
        dataset_meta.get("dataset_perturb_seq_reprocessed"),
        class_name="ms-2 align-self-center",
    )
    if provenance_badge is not None:
        title_elements.append(provenance_badge)

    if url_dataset_id:
        title_elements.append(
            html.A(
                "[more info]",
                href=f"/perturbation-catalogue/dataset/{url_dataset_id}",
                className="text-decoration-none ms-2 small align-self-center",
                style={"color": COLORS["primary"]},
            )
        )

    # Add MaveDB link for MAVE datasets
    if modality == "mave" and dataset_id:
        title_elements.append(
            html.A(
                "[MaveDB info]",
                href=f"https://mavedb.org/score-sets/{dataset_id}",
                className="text-decoration-none ms-2 small align-self-center",
                style={"color": COLORS["primary"]},
                target="_blank",
            )
        )

    title_content = html.Div(
        title_elements,
        className="d-flex align-items-baseline flex-wrap mb-2",
    )

    metadata_lines = [
        _dataset_meta_line(label, _resolve_meta_value(dataset_meta, field))
        for field, label in DATASET_METADATA_FIELDS
        if _resolve_meta_value(dataset_meta, field)
    ]

    metadata_section = (
        html.Div(metadata_lines, className="d-flex flex-column gap-1 small")
        if metadata_lines
        else None
    )

    return html.Div(
        [
            title_content,
            metadata_section,
        ],
        className="dataset-column px-2 py-2 border rounded-3 bg-white",
        style={"gridRow": f"span {span_rows}"},
    )


def _render_result_cell(
    result: Dict[str, Any],
    modality: str,
    effect_gene_source: str,
    section_id: Optional[str] = None,
):
    if modality == "perturb-seq":
        return _perturb_seq_effect(
            result.get("perturbation") or {},
            result.get("effect") or {},
            effect_gene_source,
            section_id,
        )
    return _score_effect(
        result.get("perturbation") or {},
        result.get("effect") or {},
        modality,
    )


def _perturb_seq_effect(
    perturbation: Dict[str, Any],
    effect: Dict[str, Any],
    effect_gene_source: str,
    section_id: Optional[str] = None,
) -> html.Div:
    """Render a single Perturb-Seq result row (legacy card format for non-table sections)."""
    perturbation_gene_name = _format_target(perturbation)
    effect_gene_name = effect.get("gene_name") or "N/A"

    log2fc_value = effect.get("log2fc")
    padj_value = _format_numeric(effect.get("padj"))
    base_mean_value = _format_numeric(effect.get("base_mean"))

    log2fc_display = _format_numeric(log2fc_value)
    log2fc_tile_value: Any = log2fc_display

    if isinstance(log2fc_value, (int, float)):
        if log2fc_value > 0:
            log2fc_tile_value = _arrow_value("▲", log2fc_display, GREEN)
        elif log2fc_value < 0:
            log2fc_tile_value = _arrow_value("▼", log2fc_display, RED)

    grid_items = [
        _field_tile("log2FC", log2fc_tile_value),
        _field_tile("padj", padj_value),
        _field_tile("base mean", base_mean_value),
    ]

    grid = html.Div(
        grid_items,
        className="effect-grid d-grid",
        style={
            "display": "grid",
            "gridTemplateColumns": "repeat(auto-fit, minmax(140px, 1fr))",
            "gap": "0.5rem",
        },
    )

    # For both Perturb-Seq sections, show both Perturbation and Effect gene
    if section_id in ("perturb_seq_perturbed", "perturb_seq_affected"):
        gene_section = html.Div(
            [
                html.Div(
                    [
                        html.Span("Perturbation", className="fw-light text-muted me-2"),
                        html.Span(
                            perturbation_gene_name,
                            className="h4 fw-bold mb-0 text-break",
                        ),
                    ],
                    className="d-flex flex-column flex-md-row gap-1 mb-2",
                ),
                html.Div(
                    [
                        html.Span("Effect gene", className="fw-light text-muted me-2"),
                        html.Span(
                            effect_gene_name, className="h4 fw-bold mb-0 text-break"
                        ),
                    ],
                    className="d-flex flex-column flex-md-row gap-1 mb-2",
                ),
            ],
            className="mb-2",
        )
    else:
        # For other sections, show only one gene based on effect_gene_source
        if effect_gene_source == "perturbation":
            gene_name = perturbation_gene_name
            gene_label = "Perturbation gene"
        else:
            gene_name = effect_gene_name
            gene_label = "Effect gene"
        gene_section = html.Div(
            [
                html.Span(gene_label, className="fw-light text-muted me-2"),
                html.Span(gene_name, className="h4 fw-bold mb-0 text-break"),
            ],
            className="d-flex flex-column flex-md-row gap-1 mb-2",
        )

    return html.Div(
        [
            gene_section,
            grid,
        ],
        className="effect-column px-2 py-2 border rounded-3 bg-white",
    )


def _perturb_seq_table(
    results: List[Dict[str, Any]],
    section_id: Optional[str] = None,
    download_url: Optional[str] = None,
    gsea_button_data: Optional[Dict[str, str]] = None,
    dataset_cell_types: Optional[str] = None,
    extra_controls: Optional[Any] = None,
) -> html.Div:
    """Render Perturb-Seq results as a traditional table with columns."""
    if not results:
        return html.Div(
            "No results available.",
            className="text-muted fst-italic py-2",
        )

    # Build table header
    header_row = html.Tr(
        [
            html.Th("Perturbation", className="text-start"),
            html.Th("Effect Gene", className="text-start"),
            html.Th("Log2FC", className="text-end"),
            html.Th("Padj", className="text-end"),
            html.Th("Statistical Score", className="text-start"),
            html.Th("Cell Type", className="text-start"),
        ]
    )

    # Build table rows
    table_rows = []
    for result in results:
        perturbation = result.get("perturbation") or {}
        effect = result.get("effect") or {}

        perturbation_gene_name = _format_target(perturbation)
        effect_gene_name = effect.get("gene_name") or "N/A"

        log2fc_value = effect.get("log2fc")
        log2fc_display = _format_numeric(log2fc_value)

        # Apply color styling for log2fc
        if isinstance(log2fc_value, (int, float)):
            if log2fc_value > 0:
                log2fc_cell = html.Td(
                    _arrow_value("▲", log2fc_display, GREEN),
                    className="text-end",
                )
            elif log2fc_value < 0:
                log2fc_cell = html.Td(
                    _arrow_value("▼", log2fc_display, RED),
                    className="text-end",
                )
            else:
                log2fc_cell = html.Td(log2fc_display, className="text-end")
        else:
            log2fc_cell = html.Td(log2fc_display, className="text-end")

        padj_raw = effect.get("padj")
        padj_value = _format_numeric(padj_raw)

        # Apply green color for significant padj values (<= 0.05)
        if isinstance(padj_raw, (int, float)) and padj_raw <= 0.05:
            padj_cell = html.Td(
                html.Span(padj_value, style={"color": GREEN, "fontWeight": "bold"}),
                className="text-end",
            )
        else:
            padj_cell = html.Td(padj_value, className="text-end")

        # Build Statistical Score cell (score_name: score_value)
        score_name = effect.get("score_name")
        score_value = effect.get("score_value")
        if score_name and score_value is not None:
            statistical_score = f"{score_name}: {_format_numeric(score_value)}"
        elif score_name:
            statistical_score = score_name
        elif score_value is not None:
            statistical_score = _format_numeric(score_value)
        else:
            statistical_score = "N/A"

        # Get cell type (use dataset cell_types as fallback if effect cell_type is N/A)
        cell_type = effect.get("cell_type") or dataset_cell_types or "N/A"

        table_rows.append(
            html.Tr(
                [
                    html.Td(
                        perturbation_gene_name,
                        className="text-start fw-semibold",
                    ),
                    html.Td(effect_gene_name, className="text-start fw-semibold"),
                    log2fc_cell,
                    padj_cell,
                    html.Td(statistical_score, className="text-start"),
                    html.Td(cell_type, className="text-start"),
                ]
            )
        )

    table = html.Table(
        [
            html.Thead(header_row, className="table-light"),
            html.Tbody(table_rows),
        ],
        className="table table-sm table-hover mb-0",
        style={"fontSize": "0.9rem"},
    )

    # Build content with optional download button and GSEA button
    content_children = []
    button_row_children = []

    if download_url:
        button_row_children.append(
            html.A(
                dbc.Button(
                    [
                        html.I(className="bi bi-download me-2"),
                        "Download Data",
                    ],
                    color="primary",
                    size="sm",
                    style={
                        "backgroundColor": COLORS["primary"],
                        "borderColor": COLORS["primary"],
                        "borderRadius": "6px",
                    },
                ),
                href=download_url,
                target="_blank",
                className="text-decoration-none me-2",
            )
        )

    # Add GSEA button for perturb_seq_perturbed section only
    if section_id == "perturb_seq_perturbed" and gsea_button_data:
        dataset_id = gsea_button_data.get("dataset_id", "")
        perturbed_gene = gsea_button_data.get("perturbed_gene_name", "")
        gsea_dataset_cell_type = gsea_button_data.get("dataset_cell_types") or ""
        # Generate unique ID for the popover target
        unique_key = f"{dataset_id}_{perturbed_gene}"
        gsea_icon_id = (
            f"gsea-info-icon-{hashlib.md5(unique_key.encode()).hexdigest()[:8]}"
        )
        button_row_children.extend(
            [
                dbc.Button(
                    [
                        html.I(className="bi bi-bar-chart-line me-2"),
                        "GSEA",
                    ],
                    id={
                        "type": "gsea-modal-trigger",
                        "dataset_id": dataset_id,
                        "perturbed_gene": perturbed_gene,
                        "dataset_cell_types": gsea_dataset_cell_type,
                    },
                    color="success",
                    size="sm",
                    className="me-1",
                    style={
                        "borderRadius": "6px",
                    },
                ),
                html.Span(
                    html.I(className="bi bi-question-circle"),
                    id=gsea_icon_id,
                    style={"cursor": "pointer", "color": "#6c757d"},
                ),
                dbc.Popover(
                    [
                        dbc.PopoverHeader("Pathway enrichment (GSEA)"),
                        dbc.PopoverBody(
                            "Shows biological pathways whose genes are collectively up- or down-regulated after a genetic perturbation, based on single-cell Perturb-seq data and MSigDB Hallmark gene sets."
                        ),
                    ],
                    target=gsea_icon_id,
                    trigger="click",
                    placement="bottom",
                ),
            ]
        )

    # Add extra controls if provided
    if extra_controls:
        button_row_children.append(extra_controls)

    if button_row_children:
        content_children.append(
            html.Div(
                button_row_children,
                className="mb-2 d-flex align-items-center",
            )
        )
    content_children.append(table)

    return html.Div(
        content_children,
        className="effect-column px-2 py-2 border rounded-3 bg-white",
        style={"overflowX": "auto"},
    )


def _crispr_table(
    results: List[Dict[str, Any]],
    download_url: Optional[str] = None,
) -> html.Div:
    """Render CRISPR screen results as a traditional table with columns."""
    if not results:
        return html.Div(
            "No results available.",
            className="text-muted fst-italic py-2",
        )

    # Build table header
    header_row = html.Tr(
        [
            html.Th("Perturbation", className="text-start"),
            html.Th("Score Name", className="text-start"),
            html.Th("Score Value", className="text-end"),
            html.Th("Significant", className="text-center"),
            html.Th("Significance Criteria", className="text-start"),
        ]
    )

    # Build table rows
    table_rows = []
    for result in results:
        perturbation = result.get("perturbation") or {}
        effect = result.get("effect") or {}

        perturbation_gene_name = _format_target(perturbation)
        score_name = effect.get("score_name") or "N/A"
        score_value = _format_numeric(effect.get("score_value"))
        significant = effect.get("significant")
        significance_criteria = effect.get("significance_criteria") or "N/A"

        # Apply color styling for significant field
        if significant is not None:
            significant_str = str(significant).lower()
            if significant_str == "true":
                significant_cell = html.Td(
                    html.Span(
                        str(significant), style={"color": GREEN, "fontWeight": "bold"}
                    ),
                    className="text-center",
                )
            else:
                significant_cell = html.Td(
                    html.Span(
                        str(significant), style={"color": RED, "fontWeight": "bold"}
                    ),
                    className="text-center",
                )
        else:
            significant_cell = html.Td("N/A", className="text-center")

        table_rows.append(
            html.Tr(
                [
                    html.Td(
                        perturbation_gene_name,
                        className="text-start fw-semibold",
                    ),
                    html.Td(score_name, className="text-start"),
                    html.Td(score_value, className="text-end"),
                    significant_cell,
                    html.Td(significance_criteria, className="text-start"),
                ]
            )
        )

    table = html.Table(
        [
            html.Thead(header_row, className="table-light"),
            html.Tbody(table_rows),
        ],
        className="table table-sm table-hover mb-0",
        style={"fontSize": "0.9rem"},
    )

    # Build content with optional download button
    content_children = []
    if download_url:
        content_children.append(
            html.Div(
                html.A(
                    dbc.Button(
                        [
                            html.I(className="bi bi-download me-2"),
                            "Download Data",
                        ],
                        color="primary",
                        size="sm",
                        style={
                            "backgroundColor": COLORS["primary"],
                            "borderColor": COLORS["primary"],
                            "borderRadius": "6px",
                        },
                    ),
                    href=download_url,
                    target="_blank",
                    className="text-decoration-none",
                ),
                className="mb-2",
            )
        )
    content_children.append(table)

    return html.Div(
        content_children,
        className="effect-column px-2 py-2 border rounded-3 bg-white",
        style={"overflowX": "auto"},
    )


def _score_effect(
    perturbation: Dict[str, Any], effect: Dict[str, Any], modality: str
) -> html.Div:
    pert_gene = _format_target(perturbation)
    variant = perturbation.get("name")
    score_name = effect.get("score_name")
    score_value = _format_numeric(effect.get("score_value"))
    significant = effect.get("significant")
    significance_criteria = effect.get("significance_criteria")

    headline_children = [
        html.Div(
            [
                html.Span("Perturbation", className="fw-light text-muted me-2"),
                html.Span(pert_gene, className="h4 fw-bold mb-0 text-break"),
            ],
            className="d-flex flex-column flex-md-row gap-1",
        )
    ]
    if variant:
        headline_children.append(
            html.Div(
                [
                    html.Span("Variant", className="fw-light text-muted me-2"),
                    html.Span(variant, className="fw-semibold text-break"),
                ],
                className="d-flex flex-column flex-md-row gap-1",
            )
        )

    grid_items = []
    if score_name:
        grid_items.append(_field_tile("Effect score", score_name))
    if score_value is not None:
        grid_items.append(_field_tile("Value", score_value))

    # Add significant field with conditional coloring
    if significant is not None:
        significant_str = str(significant).lower()
        if significant_str == "true":
            significant_display = html.Span(
                str(significant), style={"color": GREEN, "fontWeight": "bold"}
            )
        else:
            significant_display = html.Span(
                str(significant), style={"color": RED, "fontWeight": "bold"}
            )
        grid_items.append(_field_tile("Significant", significant_display))

    # Add significance_criteria field
    if significance_criteria is not None:
        grid_items.append(
            _field_tile("Significance criteria", str(significance_criteria))
        )

    if not grid_items:
        grid_items.append(
            _field_tile("Value", html.Span("N/A", className="text-muted"))
        )

    grid = html.Div(
        grid_items,
        style={
            "display": "grid",
            "gridTemplateColumns": "repeat(auto-fit, minmax(180px, 1fr))",
            "gap": "0.5rem",
        },
    )

    return html.Div(
        [html.Div(headline_children, className="mb-2 d-flex flex-column gap-1"), grid],
        className="effect-column px-2 py-2 border rounded-3 bg-white",
    )


def _mave_heatmap_effect(
    results: List[Dict[str, Any]],
    download_url: Optional[str] = None,
) -> html.Div:
    """Create a heatmap visualization for MAVE data showing position-based scores."""
    if not results:
        return html.Div(
            "No data available for heatmap.",
            className="effect-column px-2 py-2 border rounded-3 bg-white text-muted fst-italic",
        )

    # Build positions dictionary
    positions = defaultdict(list)

    for item in results:
        perturbation = item.get("perturbation") or {}
        effect = item.get("effect") or {}

        position = perturbation.get("position")
        if position is None:
            continue

        aa_change = perturbation.get("aa_change")
        aa_wt = perturbation.get("aa_wt")
        score_value = effect.get("score_value")

        is_reference = aa_change == "="
        aa = aa_wt if is_reference else aa_change

        if aa is None or score_value is None:
            continue

        positions[position].append(
            {
                "aa": aa,
                "score": score_value,
                "is_ref": is_reference,
            }
        )

    if not positions:
        return html.Div(
            "No valid position data for heatmap.",
            className="effect-column px-2 py-2 border rounded-3 bg-white text-muted fst-italic",
        )

    # Collect all unique amino acids across all positions
    all_aas = set()
    for lst in positions.values():
        for d in lst:
            all_aas.add(d["aa"])

    # Sort amino acids for consistent ordering
    aa_index = sorted(all_aas)

    if not aa_index:
        return html.Div(
            "No valid amino acid data for heatmap.",
            className="effect-column px-2 py-2 border rounded-3 bg-white text-muted fst-italic",
        )

    # Sort positions numerically
    sorted_positions = dict(sorted(positions.items(), key=lambda x: x[0]))

    # Build score data + ref mask dictionaries: position -> aa -> value
    # This allows us to handle different AAs per position
    data_dict = {}
    ref_mask_dict = {}

    for pos, lst in sorted_positions.items():
        pos_str = str(pos)
        # Create dictionaries for this position: aa -> score/ref
        pos_scores = {d["aa"]: d["score"] for d in lst}
        pos_refs = {d["aa"]: d["is_ref"] for d in lst}

        # Build lists for all AAs, using NaN/False for missing ones
        data_dict[pos_str] = [pos_scores.get(aa, None) for aa in aa_index]
        ref_mask_dict[pos_str] = [pos_refs.get(aa, False) for aa in aa_index]

    df = pd.DataFrame(data_dict, index=aa_index)
    ref_df = pd.DataFrame(ref_mask_dict, index=aa_index)

    # Create heatmap
    fig = px.imshow(df, color_continuous_scale="RdYlGn")

    # Disable hover tooltips
    fig.update_traces(hoverinfo="skip", hovertemplate="")

    # Add annotations for reference cells
    annotations = []
    for j, col in enumerate(ref_df.columns):
        for i, row in enumerate(ref_df.index):
            # Use iloc to get scalar value and avoid pandas boolean ambiguity
            is_ref_value = ref_df.iloc[i, j]
            # Explicitly check if value is True
            if pd.notna(is_ref_value) and is_ref_value == True:
                (
                    dict(
                        x=j,
                        y=i,
                        text="WT",
                        showarrow=False,
                        font=dict(color="black", size=10),
                    )
                )

    if annotations:
        fig.update_layout(annotations=annotations)

    # Update layout for better display
    fig.update_layout(
        title=dict(
            text="<b>Functional Score by Variant</b>",
            x=0.5,  # Center the title
            xanchor="center",
            font=dict(size=18),
        ),
        margin=dict(l=40, r=40, t=60, b=40),  # Increased top margin for title
        xaxis_title="Position",
        yaxis_title="Amino Acid",
        hovermode=False,  # Disable hover mode completely
    )

    fig.update_xaxes(tickmode="linear")
    fig.update_yaxes(tickmode="linear")

    # Build content with optional download button
    content_children = []
    if download_url:
        content_children.append(
            html.Div(
                html.A(
                    dbc.Button(
                        [
                            html.I(className="bi bi-download me-2"),
                            "Download Data",
                        ],
                        color="primary",
                        size="sm",
                        style={
                            "backgroundColor": COLORS["primary"],
                            "borderColor": COLORS["primary"],
                            "borderRadius": "6px",
                        },
                    ),
                    href=download_url,
                    target="_blank",
                    className="text-decoration-none",
                ),
                className="mb-2",
            )
        )
    content_children.append(
        dcc.Graph(
            figure=fig,
            config={"displayModeBar": False},
            style={"height": "100%", "width": "100%"},
        )
    )

    return html.Div(
        content_children,
        className="effect-column px-2 py-2 border rounded-3 bg-white",
        style={
            "height": "500px" if download_url else "450px",
            "minHeight": "450px",
            "overflow": "hidden",
            "marginBottom": "0.5rem",
        },
    )


def _dataset_meta_line(label: str, value: Optional[Any]) -> html.Div:
    if value is None or (isinstance(value, list) and not value):
        return html.Div()

    if isinstance(value, list):
        pretty_value = ", ".join([_capitalize_value(str(v)) for v in value])
    else:
        pretty_value = _capitalize_value(str(value))

    # Get badge color for this field, default to "secondary" if not found
    badge_color = METADATA_FIELD_COLORS.get(label, "secondary")
    return html.Div(
        [
            dbc.Badge(
                label,
                color=badge_color,
                className="text-uppercase small",
                style={"fontSize": "0.7rem", "width": "fit-content"},
            ),
            html.Span(pretty_value, className="fw-semibold small text-break"),
        ],
        className="d-flex flex-column flex-md-row align-items-md-center gap-1 gap-md-2",
    )


def _grid_message(
    message: str,
    tone: str = "info",
    span_all: bool = True,
    columns: Optional[str] = None,
) -> html.Div:
    color = {"error": "#dc3545", "info": "#6c757d"}.get(tone, "#6c757d")
    style = {
        "gridColumn": columns or ("1 / -1" if span_all else "auto"),
        "color": color,
    }
    return html.Div(message, className="fst-italic py-2", style=style)


def _truncation_notice(limit_value: int) -> html.Div:
    return html.Div(
        f"Displaying top {limit_value} results",
        className="text-center fst-italic",
        style={
            "gridColumn": "2 / 3",
            "backgroundColor": "#f8f9fa",
            "padding": "0.5rem 0.75rem",
            "borderRadius": "0.5rem",
            "justifySelf": "end",
            "textAlign": "right",
        },
    )


def _dataset_separator() -> html.Div:
    return html.Div(
        "",
        style={
            "gridColumn": "1 / -1",
            "borderBottom": "2px solid #dee2e6",
            "margin": "0.5rem 0",
        },
    )


def _format_dataset_id(dataset_id: Optional[str]) -> str:
    if not dataset_id:
        return "Dataset"
    return dataset_id


def _dataset_id_to_url_format(dataset_id: Optional[str]) -> Optional[str]:
    """Convert dataset ID to URL format (e.g., 'Depmap ACH000558' -> 'depmap_ACH000558')."""
    if not dataset_id:
        return None
    # Split by space or underscore
    parts = dataset_id.replace("_", " ").split()
    if not parts:
        return None
    # Lowercase first part, keep rest as is, join with underscore
    if len(parts) == 1:
        return parts[0].lower()
    return "_".join([parts[0].lower()] + parts[1:])


def _resolve_meta_value(meta: Dict[str, Any], field: str) -> Optional[Any]:
    if field in meta and meta[field] not in (None, ""):
        return meta[field]
    for fallback in DATASET_FIELD_FALLBACKS.get(field, []):
        value = meta.get(fallback)
        if value not in (None, ""):
            return value
    return None


def _format_numeric(value: Any) -> str:
    if value is None:
        return "N/A"
    formatted = format_number(value)
    return formatted.replace("-", MINUS)


def _field_tile(label: str, value: Any) -> html.Div:
    return html.Div(
        [
            html.Div(label, className="fw-light text-muted small text-uppercase"),
            html.Div(value, className="fw-semibold small"),
        ],
        className="d-flex flex-column gap-1 border rounded-3 px-2 py-1 bg-light",
        style={"minHeight": "60px"},
    )


def _arrow_value(symbol: str, value: str, color: str) -> html.Span:
    return html.Span(
        [
            html.Span(symbol, style={"color": color, "marginRight": "0.25rem"}),
            value,
        ],
        style={"color": color},
    )


def _capitalize_value(value: str) -> str:
    return value[:1].upper() + value[1:] if value else value


# Public wrapper functions for use in other modules


def mave_heatmap(
    results: List[Dict[str, Any]],
    download_url: Optional[str] = None,
) -> html.Div:
    """Create a heatmap visualization for MAVE data showing position-based scores.

    Public wrapper for _mave_heatmap_effect.

    Args:
        results: List of result dicts with perturbation and effect data.
        download_url: Optional URL for downloading the data.

    Returns:
        A Dash html.Div containing the heatmap visualization.
    """
    return _mave_heatmap_effect(results, download_url)


def perturb_seq_table(
    results: List[Dict[str, Any]],
    section_id: Optional[str] = None,
    download_url: Optional[str] = None,
    dataset_cell_types: Optional[str] = None,
    extra_controls: Optional[Any] = None,
) -> html.Div:
    """Render Perturb-Seq results as a table.

    Public wrapper for _perturb_seq_table.

    Args:
        results: List of result dicts with perturbation and effect data.
        section_id: Section identifier for styling.
        download_url: Optional URL for downloading the data.
        dataset_cell_types: Fallback cell type from dataset metadata.
        extra_controls: Optional extra controls to render alongside the download button.

    Returns:
        A Dash html.Div containing the table.
    """
    return _perturb_seq_table(
        results,
        section_id=section_id,
        download_url=download_url,
        gsea_button_data=None,
        dataset_cell_types=dataset_cell_types,
        extra_controls=extra_controls,
    )


def crispr_table(
    results: List[Dict[str, Any]],
    download_url: Optional[str] = None,
) -> html.Div:
    """Render CRISPR screen results as a table.

    Public wrapper for _crispr_table.

    Args:
        results: List of result dicts with perturbation and effect data.
        download_url: Optional URL for downloading the data.

    Returns:
        A Dash html.Div containing the table.
    """
    return _crispr_table(results, download_url)
