"""Reusable target data table component with CSS grid layout."""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional

from dash import html

from utils import format_number

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
    ("dataset_tissue", "Tissue"),
    ("dataset_cell_type", "Cell type"),
    ("dataset_cell_line", "Cell line"),
    ("dataset_library_perturbation_type", "Library perturbation"),
    ("dataset_disease", "Disease"),
    ("dataset_sex", "Sex"),
    ("dataset_developmental_stage", "Developmental stage"),
]
DATASET_FIELD_FALLBACKS = {
    "dataset_id": ["id"],
    "dataset_tissue": ["tissue"],
    "dataset_cell_type": ["cell_type"],
    "dataset_cell_line": ["cell_line"],
    "dataset_library_perturbation_type": [
        "library_perturbation_type",
        "library_type",
    ],
    "dataset_disease": ["disease"],
    "dataset_sex": ["sex"],
    "dataset_developmental_stage": ["developmental_stage"],
}

GREEN = "#2acc06"
RED = "#ff4824"
MINUS = "−"


def TargetDataTable(
    data: Optional[List[Dict[str, Any]]],
    modality: str,
    table_id: Optional[str] = None,
    rows_per_dataset_limit: Optional[int] = None,
    empty_message: str = "No results found.",
    error_message: Optional[str] = None,
    dataset_control_factory: GridControlFactory = None,
    effect_gene_source: str = "effect",
):
    """Render the reusable data table."""
    datasets = data or []

    grid_children: List[Any] = []
    grid_children.extend(_build_header_cells())

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


def _build_header_cells() -> List[html.Div]:
    """Create the table header with compact titles."""
    return [
        html.Div(html.H3(title, className="fw-semibold mb-1"), className="pb-1")
        for title in HEADER_TITLES
    ]


def _build_dataset_rows(
    entry: Dict[str, Any],
    modality: str,
    rows_per_dataset_limit: Optional[int],
    dataset_control_factory: GridControlFactory,
    effect_gene_source: str,
) -> List[Any]:
    dataset_meta = entry.get("dataset") or {}
    dataset_id = _resolve_meta_value(dataset_meta, "dataset_id") or "Dataset"
    results = entry.get("results") or []
    row_span = max(len(results), 1)

    children: List[Any] = [
        _render_dataset_cell(dataset_meta, row_span),
    ]

    if results:
        for result in results:
            children.append(_render_result_cell(result, modality, effect_gene_source))
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


def _render_dataset_cell(dataset_meta: Dict[str, Any], span_rows: int):
    formatted_id = _format_dataset_id(_resolve_meta_value(dataset_meta, "dataset_id"))
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
            html.Div(formatted_id, className="h4 fw-semibold mb-2 text-break"),
            metadata_section,
        ],
        className="dataset-column px-2 py-2 border rounded-3 bg-white",
        style={"gridRow": f"span {span_rows}"},
    )


def _render_result_cell(result: Dict[str, Any], modality: str, effect_gene_source: str):
    if modality == "perturb-seq":
        return _perturb_seq_effect(
            result.get("perturbation") or {},
            result.get("effect") or {},
            effect_gene_source,
        )
    return _score_effect(
        result.get("perturbation") or {},
        result.get("effect") or {},
        modality,
    )


def _perturb_seq_effect(
    perturbation: Dict[str, Any], effect: Dict[str, Any], effect_gene_source: str
) -> html.Div:
    if effect_gene_source == "perturbation":
        gene_name = perturbation.get("gene_name") or "N/A"
    else:
        gene_name = effect.get("gene_name") or "N/A"
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

    return html.Div(
        [
            html.Div(
                [
                    html.Span("Effect gene", className="fw-light text-muted me-2"),
                    html.Span(gene_name, className="h4 fw-bold mb-0 text-break"),
                ],
                className="d-flex flex-column flex-md-row gap-1 mb-2",
            ),
            grid,
        ],
        className="effect-column px-2 py-2 border rounded-3 bg-white",
    )


def _score_effect(
    perturbation: Dict[str, Any], effect: Dict[str, Any], modality: str
) -> html.Div:
    pert_gene = perturbation.get("gene_name") or "N/A"
    variant = perturbation.get("name")
    score_name = effect.get("score_name")
    score_value = _format_numeric(effect.get("score_value"))

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


def _dataset_meta_line(label: str, value: Optional[Any]) -> html.Div:
    if value is None:
        return html.Div()
    pretty_value = _capitalize_value(str(value))
    return html.Div(
        [
            html.Span(
                label,
                className="fw-light text-muted text-uppercase small me-2",
            ),
            html.Span(pretty_value, className="fw-semibold small text-break"),
        ]
    )


def _grid_message(
    message: str,
    tone: str = "info",
    span_all: bool = True,
    columns: Optional[str] = None,
) -> html.Div:
    color = {"error": "#dc3545", "info": "#6c757d"}.get(tone, "#6c757d")
    style = {"gridColumn": columns or ("1 / -1" if span_all else "auto"), "color": color}
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
    formatted = dataset_id.replace("_", " ")
    return formatted[:1].upper() + formatted[1:]


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


