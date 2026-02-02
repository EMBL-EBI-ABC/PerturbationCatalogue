"""Home page for Perturbation Search"""

import dash
from dash import dcc, html, Input, Output, State, callback, ALL, callback_context
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from urllib.parse import quote
from utils import (
    COLORS,
    DATA_MODALITIES_COLOURS,
    FACET_FIELDS,
    fetch_search_results,
    fetch_all_search_results,
    get_landing_page_summary,
    results_store,
    format_value,
)

SEARCH_RESULTS_PAGE_SIZE = 15

# Mapping from canonical (target) field names to dataset index field names.
# Used for client-side filtering of dataset results.
TARGET_TO_DATASET_FIELD = {
    "license": "license_labels",
    "data_modalities": "data_modalities",
    "tissues_tested": "tissue_labels",
    "cell_types_tested": "cell_type_labels",
    "cell_lines_tested": "cell_line_labels",
    "sex_tested": "sex_labels",
    "developmental_stages_tested": "developmental_stage_labels",
    "diseases_tested": "disease_labels",
}


# Register this page with Dash Pages
dash.register_page(__name__, path="/")


def _format_count(value):
    """Format numeric counts for display, falling back to '0' when missing."""
    if value in (None, "", "N/A"):
        return "0"
    if isinstance(value, (int, float)):
        return f"{int(value):,}" if isinstance(value, int) else format_value(value)
    return format_value(value)


def _render_search_results(results):
    """Render the search results table for the provided entries."""
    # Clear previous cache and repopulate with current results
    results_store.clear()

    if not results:
        return dbc.Alert(
            [
                html.I(className="bi bi-info-circle me-2"),
                "No results found. Try another target name.",
            ],
            color="info",
            className="shadow-sm",
            style={"borderRadius": "10px"},
        )

    header = html.Thead(
        html.Tr(
            [
                html.Th(
                    [
                        "Target name ",
                        html.Span(
                            html.I(className="bi bi-question-circle"),
                            id="target-info-icon",
                            style={"cursor": "pointer", "color": "#6c757d"},
                        ),
                        dbc.Popover(
                            [
                                dbc.PopoverHeader("Target Gene Aggregation"),
                                dbc.PopoverBody(
                                    "Results are aggregated by target gene name. Each row represents all datasets and modalities available for a specific target, combining data from Perturb-seq, CRISPR screens, and MAVE experiments across different tissues, cell types, and experimental conditions."
                                ),
                            ],
                            target="target-info-icon",
                            trigger="click",
                            placement="bottom",
                        ),
                    ],
                    className="fw-semibold",
                ),
                html.Th("Perturb-seq datasets", className="fw-semibold text-center"),
                html.Th(
                    "CRISPR-screen datasets",
                    className="fw-semibold text-center",
                ),
                html.Th(
                    "MAVE datasets",
                    className="fw-semibold text-center",
                ),
                html.Th(
                    [
                        "Top GSEA Terms ",
                        html.Span(
                            html.I(className="bi bi-question-circle"),
                            id="gsea-info-icon",
                            style={"cursor": "pointer", "color": "#6c757d"},
                        ),
                        dbc.Popover(
                            [
                                dbc.PopoverHeader("Pathway enrichment (GSEA)"),
                                dbc.PopoverBody(
                                    "Shows biological pathways whose genes are collectively up- or down-regulated after a genetic perturbation, based on single-cell Perturb-seq data and MSigDB Hallmark gene sets."
                                ),
                            ],
                            target="gsea-info-icon",
                            trigger="click",
                            placement="bottom",
                        ),
                    ],
                    className="fw-semibold",
                ),
                html.Th("Data Modalities", className="fw-semibold"),
            ],
            style={"backgroundColor": "#f1f3f5"},
        )
    )

    rows = []
    for record in results:
        symbol = record.get("perturbed_target_symbol", "N/A")
        results_store[symbol] = record

        n_sc_perturb_seq = _format_count(record.get("n_perturb_seq"))
        n_sc_perturb_seq_up = _format_count(record.get("n_sig_perturb_pairs_up"))
        n_sc_perturb_seq_down = _format_count(record.get("n_sig_perturb_pairs_down"))
        n_crispr = _format_count(record.get("n_crispr"))
        n_sig_crispr = _format_count(record.get("n_sig_crispr"))
        n_mave = _format_count(record.get("n_mave"))
        top_gsea_terms = record.get("top_gsea_terms") or []
        data_modalities = record.get("data_modalities") or []

        gsea_badges = []
        for term in top_gsea_terms[:5]:  # Limit to 5 terms
            # Truncate long terms for display
            display_term = term if len(term) <= 25 else term[:22] + "..."
            gsea_badges.append(
                html.Span(
                    dbc.Badge(
                        display_term,
                        color="light",
                        className="me-1 mb-1",
                        style={
                            "fontSize": "0.65rem",
                            "border": "1px solid #6c757d",
                            "color": "#495057",
                            "backgroundColor": "#f8f9fa",
                            "whiteSpace": "nowrap",
                            "overflow": "hidden",
                            "textOverflow": "ellipsis",
                            "maxWidth": "150px",
                            "display": "inline-block",
                        },
                    ),
                    title=term,  # Full term shown on hover
                )
            )
        if not gsea_badges:
            gsea_badges = [html.Span("No significant hits were found", className="text-muted", style={"fontSize": "0.8rem"})]

        modalities_badges = []
        for modality in data_modalities:
            badge_color = DATA_MODALITIES_COLOURS.get(modality, COLORS["primary"])
            modalities_badges.append(
                dbc.Badge(
                    modality,
                    color="light",
                    className="me-1 mb-1",
                    style={
                        "fontSize": "0.75rem",
                        "border": f"1px solid {badge_color}",
                        "color": badge_color,
                        "backgroundColor": f"{badge_color}1A",
                    },
                )
            )
        if not modalities_badges:
            modalities_badges = [html.Span("N/A", className="text-muted")]

        rows.append(
            html.Tr(
                [
                    html.Td(
                        dcc.Link(
                            symbol,
                            href=f"/perturbation-catalogue/target/{quote(symbol, safe='')}",
                            className="text-decoration-none fw-semibold",
                            style={"color": COLORS["primary"]},
                        )
                    ),
                    html.Td(
                        [
                            html.Div(
                                [
                                    html.Span(
                                        f"Datasets: {n_sc_perturb_seq} ",
                                        style={
                                            "marginRight": "5px",
                                        },
                                    ),
                                    html.Span(
                                        [
                                            html.Span("("),
                                            html.Span(
                                                "Sig. genes up ↑: ",
                                                style={
                                                    "color": "#1f9d55",
                                                    "marginRight": "0.25rem",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(
                                                n_sc_perturb_seq_up,
                                                style={
                                                    "color": "#1f9d55",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(";"),
                                        ],
                                        className="d-inline-flex align-items-center me-3",
                                    ),
                                    html.Span(
                                        [
                                            html.Span(
                                                "Sig. genes down ↓: ",
                                                style={
                                                    "color": "#c53030",
                                                    "marginRight": "0.25rem",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(
                                                n_sc_perturb_seq_down,
                                                style={
                                                    "color": "#c53030",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(")"),
                                        ],
                                        className="d-inline-flex align-items-center",
                                    ),
                                ],
                                className="d-flex justify-content-center flex-wrap",
                            )
                        ],
                        className="text-center",
                    ),
                    html.Td(
                        html.Div(
                            [
                                html.Span(
                                    f"Datasets: {n_crispr}",
                                    style={
                                        "marginRight": "5px",
                                    },
                                ),
                                html.Span(
                                    f"Sig. hits: {n_sig_crispr}",
                                    style={
                                        "color": "#1f9d55",
                                        "marginRight": "0.25rem",
                                        "fontWeight": "600",
                                    },
                                ),
                            ],
                            className="d-flex justify-content-center flex-wrap",
                        ),
                        className="text-center",
                    ),
                    html.Td(n_mave, className="text-center"),
                    html.Td(
                        html.Div(gsea_badges, className="d-flex flex-wrap"),
                        style={
                            "maxWidth": "180px",
                            "minWidth": "120px",
                            "overflow": "hidden",
                        },
                    ),
                    html.Td(modalities_badges),
                ]
            )
        )

    table = dbc.Table(
        [header, html.Tbody(rows)],
        bordered=False,
        hover=True,
        responsive=True,
        striped=False,
        className="align-middle shadow-sm",
        style={"borderRadius": "12px", "overflow": "hidden"},
    )

    return table


def _recalculate_facets(filtered_results, original_facets, search_mode="targets"):
    """Recalculate facet counts based on currently filtered results."""
    is_dataset_mode = search_mode == "datasets"
    recalculated = {}
    for field in FACET_FIELDS:
        value_counts = {}
        # In dataset mode, result fields use different names (e.g. cell_type_labels vs cell_types_tested)
        result_field = TARGET_TO_DATASET_FIELD.get(field, field) if is_dataset_mode else field
        for result in filtered_results:
            field_value = result.get(result_field)
            if field_value is None:
                continue
            if isinstance(field_value, list):
                for v in field_value:
                    v_str = str(v).strip()
                    if v_str:
                        value_counts[v_str] = value_counts.get(v_str, 0) + 1
            else:
                v_str = str(field_value).strip()
                if v_str:
                    value_counts[v_str] = value_counts.get(v_str, 0) + 1

        # Preserve original ordering, update counts
        original_values = original_facets.get(field, [])
        facet_list = []
        seen = set()
        for item in original_values:
            val = item.get("value")
            if val is not None:
                val_str = str(val).strip()
                count = value_counts.get(val_str, 0)
                facet_list.append({"value": val_str, "count": count})
                seen.add(val_str)
        for val_str, count in value_counts.items():
            if val_str not in seen:
                facet_list.append({"value": val_str, "count": count})

        recalculated[field] = facet_list

    return recalculated


def _render_dataset_search_results(results):
    """Render the search results table for dataset search results."""
    if not results:
        return dbc.Alert(
            [
                html.I(className="bi bi-info-circle me-2"),
                "No datasets found. Try another search term.",
            ],
            color="info",
            className="shadow-sm",
            style={"borderRadius": "10px"},
        )

    header = html.Thead(
        html.Tr(
            [
                html.Th("Dataset ID", className="fw-semibold"),
                html.Th("Study Title", className="fw-semibold"),
                html.Th("Study Year", className="fw-semibold text-center"),
                html.Th("Data Modality", className="fw-semibold"),
            ],
            style={"backgroundColor": "#f1f3f5"},
        )
    )

    rows = []
    for record in results:
        dataset_id = record.get("dataset_id", "N/A")
        study_title = record.get("study_title") or record.get("experiment_title") or "N/A"
        # Truncate long titles
        display_title = study_title if len(study_title) <= 60 else study_title[:57] + "..."
        modality = record.get("data_modalities") or []
        year = record.get("study_year") or "N/A"

        modality_badges = []
        for mod in modality:
            badge_color = DATA_MODALITIES_COLOURS.get(mod, COLORS["primary"])
            modality_badges.append(
                dbc.Badge(
                    mod,
                    color="light",
                    className="me-1 mb-1",
                    style={
                        "fontSize": "0.75rem",
                        "border": f"1px solid {badge_color}",
                        "color": badge_color,
                        "backgroundColor": f"{badge_color}1A",
                    },
                )
            )
        if not modality_badges:
            modality_badges = [html.Span("N/A", className="text-muted")]

        rows.append(
            html.Tr(
                [
                    html.Td(
                        dcc.Link(
                            dataset_id,
                            href=f"/perturbation-catalogue/dataset/{dataset_id}",
                            className="text-decoration-none fw-semibold",
                            style={"color": COLORS["primary"]},
                        )
                    ),
                    html.Td(
                        html.Span(display_title, title=study_title),
                        style={"maxWidth": "300px"},
                    ),
                    html.Td(str(year), className="text-center"),
                    html.Td(html.Div(modality_badges, className="d-flex flex-wrap")),
                ]
            )
        )

    table = dbc.Table(
        [header, html.Tbody(rows)],
        bordered=False,
        hover=True,
        responsive=True,
        striped=False,
        className="align-middle shadow-sm",
        style={"borderRadius": "12px", "overflow": "hidden"},
    )

    return table


def _filter_placeholder(message="Search to enable filters."):
    """Placeholder message when filters are unavailable."""
    return dbc.Alert(
        [html.I(className="bi bi-funnel me-2"), message],
        color="light",
        className="shadow-sm",
        style={"borderRadius": "10px"},
    )


def _format_summary_number(value):
    if value is None:
        return "N/A"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:,.2f}"
    return str(value)


SUMMARY_BAR_COLORS = [
    "#007B53",
    "#193F90",
    "#A6093D",
    "#563D82",
    "#3B6FB6",
    "#F0A202",
    "#0A5032",
    "#5B7FC7",
    "#8B6FA8",
    "#6B9BD4",
]


def _build_summary_list_section(title, items, max_items=6):
    items = items or []
    total_count = len(items)

    if total_count == 0:
        display_items = []
    elif total_count < 7:
        display_items = items  # show all when using pie chart
    else:
        display_limit = max(max_items, 10)
        display_items = items[:display_limit]

    if not display_items:
        content = html.Div("No data available", className="text-muted small")
    else:
        labels = [entry.get("value", "N/A") for entry in display_items]
        values = [entry.get("n_datasets", 0) for entry in display_items]

        # Determine chart type
        if total_count < 7:
            fig = go.Figure(
                data=[
                    go.Pie(
                        labels=labels,
                        values=values,
                        hole=0.3,
                        textinfo="label+percent",
                        hoverinfo="skip",
                    )
                ]
            )
            fig.update_layout(
                margin=dict(l=10, r=10, t=10, b=10),
                legend=dict(
                    orientation="h",
                    yanchor="bottom",
                    y=-0.25,
                    x=0.5,
                    xanchor="center",
                    font=dict(size=10),
                ),
                height=360,
            )
        else:
            bar_colors = [
                SUMMARY_BAR_COLORS[i % len(SUMMARY_BAR_COLORS)]
                for i in range(len(display_items))
            ]
            fig = go.Figure(
                data=[
                    go.Bar(
                        x=values,
                        y=labels,
                        orientation="h",
                        text=[f"{v:,}" for v in values],
                        textposition="auto",
                        hoverinfo="skip",
                        marker=dict(color=bar_colors),
                    )
                ]
            )
            fig.update_layout(
                margin=dict(l=120, r=10, t=10, b=10),
                height=350,
                yaxis=dict(autorange="reversed"),
            )

        # Add indicator if we truncated items
        chart = dcc.Graph(
            figure=fig, config={"displayModeBar": False}, className="summary-chart"
        )
        content = html.Div(chart)

    return dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(title, className="fw-semibold bg-light"),
                dbc.CardBody(content, className="summary-card-body"),
            ],
            className="h-100 shadow-sm summary-card",
        ),
        xs=12,
        lg=4,
        className="mb-4",
    )


def _build_summary_component(summary_data, error=None):
    """Build the landing page summary component."""
    if error:
        return dbc.Alert(
            [
                html.H4("Error", className="alert-heading"),
                html.P("Summary data is currently unavailable due to an error:"),
                html.Hr(),
                html.P(error, className="mb-0"),
            ],
            color="danger",
            className="shadow-sm",
            style={"borderRadius": "10px"},
        )

    if not summary_data:
        return dbc.Alert(
            [
                html.I(className="bi bi-info-circle me-2"),
                "Summary data is currently unavailable.",
            ],
            color="light",
            className="shadow-sm",
            style={"borderRadius": "10px"},
        )

    # Build Datasets card with year range
    datasets_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.Div(
                            [
                                html.H6(
                                    "Datasets",
                                    className="text-uppercase text-muted mb-0",
                                ),
                                html.I(
                                    className="bi bi-database",
                                    style={"fontSize": "1.5rem", "color": "#6c757d"},
                                ),
                            ],
                            className="d-flex justify-content-between align-items-center mb-1",
                        ),
                        html.Div(
                            [
                                html.H3(
                                    _format_summary_number(
                                        summary_data.get("n_datasets")
                                    ),
                                    className="mb-0 fw-bold d-inline-block me-2",
                                ),
                                html.Span(
                                    f"(spanning {_format_summary_number(summary_data.get('min_year'))} - {_format_summary_number(summary_data.get('max_year'))})",
                                    className="text-muted",
                                    style={"fontSize": "0.9rem"},
                                ),
                            ],
                            className="d-flex align-items-baseline",
                        ),
                    ]
                )
            ],
            className="shadow-sm summary-stat-card",
        ),
        xs=6,
        md=3,
        lg=3,
        className="mb-4",
    )

    # Build Targets card
    targets_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.Div(
                            [
                                html.H6(
                                    "Targets",
                                    className="text-uppercase text-muted mb-0",
                                ),
                                html.I(
                                    className="bi bi-bullseye",
                                    style={"fontSize": "1.5rem", "color": "#6c757d"},
                                ),
                            ],
                            className="d-flex justify-content-between align-items-center mb-1",
                        ),
                        html.H3(
                            _format_summary_number(summary_data.get("n_targets")),
                            className="mb-0 fw-bold",
                        ),
                    ]
                )
            ],
            className="shadow-sm summary-stat-card",
        ),
        xs=6,
        md=3,
        lg=3,
        className="mb-4",
    )

    # Build Unique tissues card
    tissues_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.Div(
                            [
                                html.H6(
                                    "Unique tissues",
                                    className="text-uppercase text-muted mb-0",
                                ),
                                html.I(
                                    className="bi bi-universal-access-circle",
                                    style={"fontSize": "1.5rem", "color": "#6c757d"},
                                ),
                            ],
                            className="d-flex justify-content-between align-items-center mb-1",
                        ),
                        html.H3(
                            _format_summary_number(summary_data.get("n_tissues")),
                            className="mb-0 fw-bold",
                        ),
                    ]
                )
            ],
            className="shadow-sm summary-stat-card",
        ),
        xs=6,
        md=3,
        lg=3,
        className="mb-4",
    )

    # Build Unique cell types card
    cell_types_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.Div(
                            [
                                html.H6(
                                    "Unique cell types",
                                    className="text-uppercase text-muted mb-0",
                                ),
                                html.I(
                                    className="bi bi-puzzle",
                                    style={"fontSize": "1.5rem", "color": "#6c757d"},
                                ),
                            ],
                            className="d-flex justify-content-between align-items-center mb-1",
                        ),
                        html.H3(
                            _format_summary_number(summary_data.get("n_cell_types")),
                            className="mb-0 fw-bold",
                        ),
                    ]
                )
            ],
            className="shadow-sm summary-stat-card",
        ),
        xs=6,
        md=3,
        lg=3,
        className="mb-4",
    )

    # Build Unique cell lines card
    cell_lines_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.Div(
                            [
                                html.H6(
                                    "Unique cell lines",
                                    className="text-uppercase text-muted mb-0",
                                ),
                                html.I(
                                    className="bi bi-puzzle-fill",
                                    style={"fontSize": "1.5rem", "color": "#6c757d"},
                                ),
                            ],
                            className="d-flex justify-content-between align-items-center mb-1",
                        ),
                        html.H3(
                            _format_summary_number(summary_data.get("n_cell_lines")),
                            className="mb-0 fw-bold",
                        ),
                    ]
                )
            ],
            className="shadow-sm summary-stat-card",
        ),
        xs=6,
        md=3,
        lg=3,
        className="mb-4",
    )

    # Build Unique diseases card
    diseases_card = dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.Div(
                            [
                                html.H6(
                                    "Unique diseases",
                                    className="text-uppercase text-muted mb-0",
                                ),
                                html.I(
                                    className="bi bi-virus2",
                                    style={"fontSize": "1.5rem", "color": "#6c757d"},
                                ),
                            ],
                            className="d-flex justify-content-between align-items-center mb-1",
                        ),
                        html.H3(
                            _format_summary_number(summary_data.get("n_diseases")),
                            className="mb-0 fw-bold",
                        ),
                    ]
                )
            ],
            className="shadow-sm summary-stat-card",
        ),
        xs=6,
        md=3,
        lg=3,
        className="mb-4",
    )

    stat_cards = dbc.Row(
        [
            datasets_card,
            targets_card,
            tissues_card,
            cell_types_card,
            cell_lines_card,
            diseases_card,
        ]
    )

    list_sections = [
        ("Top Modalities", summary_data.get("top_modalities")),
        ("Top Tissues", summary_data.get("top_tissues")),
        ("Top Cell Types", summary_data.get("top_cell_types")),
        ("Top Cell Lines", summary_data.get("top_cell_lines")),
        ("Top Perturbation Types", summary_data.get("top_perturbation_types")),
        ("Top Diseases", summary_data.get("top_diseases")),
        ("Top Sexes", summary_data.get("top_sexes")),
        ("Top Developmental Stages", summary_data.get("top_dev_stages")),
    ]

    list_rows = []
    columns = []
    for idx, (title, items) in enumerate(list_sections, start=1):
        columns.append(_build_summary_list_section(title, items))
        if idx % 3 == 0 or idx == len(list_sections):
            list_rows.append(dbc.Row(columns, className="gy-4"))
            columns = []

    return html.Div(
        [
            html.Div(
                [
                    html.H2("Perturbation Catalogue at a glance", className="mb-3"),
                    html.P(
                        "Explore high-level statistics across the catalogue before diving into specific targets.",
                        className="text-muted",
                    ),
                ],
                className="mb-4",
            ),
            stat_cards,
            html.Div(list_rows, className="summary-lists mt-2"),
        ],
        className="summary-container",
    )


def _build_filter_controls(facets, selected_filters=None, facet_fields=None, search_mode="targets"):
    """Build filter controls for facet fields."""
    if not facets:
        return _filter_placeholder()

    if selected_filters is None:
        selected_filters = {}

    if facet_fields is None:
        facet_fields = FACET_FIELDS

    # Icon mapping for facet fields (covers both target and dataset field names)
    field_icons = {
        "license": "bi-award-fill",
        "data_modalities": "bi-database",
        "tissues_tested": "bi-universal-access-circle",
        "cell_types_tested": "bi-puzzle",
        "cell_lines_tested": "bi-puzzle-fill",
        "diseases_tested": "bi-virus2",
        "sex_tested": "bi-gender-ambiguous",
        "developmental_stages_tested": "bi-graph-up-arrow",
        # Dataset-mode fields
        "license_labels": "bi-award-fill",
        "library_perturbation_type_labels": "bi-database",
        "tissue_labels": "bi-universal-access-circle",
        "cell_type_labels": "bi-puzzle",
        "cell_line_labels": "bi-puzzle-fill",
        "disease_labels": "bi-virus2",
        "sex_labels": "bi-gender-ambiguous",
        "developmental_stage_labels": "bi-graph-up-arrow",
    }

    # Explanations for facet fields (shown in popover)
    field_explanations = {
        "license": {
            "header": "License Filter",
            "body": "Filter targets by the license types of their underlying datasets. Selecting a license shows all targets that have at least one dataset released under that license. Note: A single target may have datasets with different licenses.",
        },
        "data_modalities": {
            "header": "Data Modalities Filter",
            "body": "Filter targets by experimental approach (Perturb-seq, CRISPR screens, MAVE). Selecting a modality shows all targets that have data from that experimental type. A target may appear in multiple modalities.",
        },
        "tissues_tested": {
            "header": "Tissues Filter",
            "body": "Filter targets by the tissues used in experiments. Selecting a tissue shows all targets that have been studied in at least one dataset using that tissue. The aggregated results for a target may include data from multiple tissues.",
        },
        "cell_types_tested": {
            "header": "Cell Types Filter",
            "body": "Filter targets by cell types used in experiments. Selecting a cell type shows all targets studied in at least one dataset using that cell type. A single target's aggregated data may span multiple cell types.",
        },
        "cell_lines_tested": {
            "header": "Cell Lines Filter",
            "body": "Filter targets by cell lines used in experiments. Selecting a cell line shows all targets studied in at least one dataset using that cell line. A target's aggregated results may include data from multiple cell lines.",
        },
        "diseases_tested": {
            "header": "Diseases Filter",
            "body": "Filter targets by disease context of experiments. Selecting a disease shows all targets studied in at least one dataset related to that disease. A target may have been studied across multiple disease contexts.",
        },
        "sex_tested": {
            "header": "Sex Filter",
            "body": "Filter targets by the biological sex of samples used in experiments. Selecting a sex shows all targets studied in at least one dataset using samples of that sex. A target's aggregated data may include both sexes.",
        },
        "developmental_stages_tested": {
            "header": "Developmental Stages Filter",
            "body": "Filter targets by developmental stage of samples. Selecting a stage shows all targets studied in at least one dataset at that developmental stage. A target may have data across multiple developmental stages.",
        },
        # Dataset-mode field explanations
        "license_labels": {
            "header": "License Filter",
            "body": "Filter datasets by the license types under which they are released. Selecting a license shows all datasets released under that license.",
        },
        "library_perturbation_type_labels": {
            "header": "Data Modalities Filter",
            "body": "Filter datasets by experimental approach (Perturb-seq, CRISPR screens, MAVE). Selecting a modality shows all datasets that have data from that experimental type.",
        },
        "tissue_labels": {
            "header": "Tissues Filter",
            "body": "Filter datasets by the tissues used in experiments. Selecting a tissue shows all datasets that used that tissue.",
        },
        "cell_type_labels": {
            "header": "Cell Types Filter",
            "body": "Filter datasets by cell types used in experiments. Selecting a cell type shows all datasets that used that cell type.",
        },
        "cell_line_labels": {
            "header": "Cell Lines Filter",
            "body": "Filter datasets by cell lines used in experiments. Selecting a cell line shows all datasets that used that cell line.",
        },
        "disease_labels": {
            "header": "Diseases Filter",
            "body": "Filter datasets by disease context of experiments. Selecting a disease shows all datasets related to that disease.",
        },
        "sex_labels": {
            "header": "Sex Filter",
            "body": "Filter datasets by the biological sex of samples used in experiments. Selecting a sex shows all datasets using samples of that sex.",
        },
        "developmental_stage_labels": {
            "header": "Developmental Stages Filter",
            "body": "Filter datasets by developmental stage of samples. Selecting a stage shows all datasets at that developmental stage.",
        },
    }

    controls = []
    for field in facet_fields:
        values = facets.get(field, [])
        if not values:
            continue

        # Custom display names for specific fields
        display_name_map = {
            # Target-mode fields
            "license": "License",
            "data_modalities": "Data Modalities",
            "tissues_tested": "Tissues",
            "cell_types_tested": "Cell Types",
            "cell_lines_tested": "Cell Lines",
            "sex_tested": "Sex",
            "developmental_stages_tested": "Developmental Stages",
            "diseases_tested": "Diseases",
            # Dataset-mode fields
            "license_labels": "License",
            "library_perturbation_type_labels": "Data Modalities",
            "tissue_labels": "Tissues",
            "cell_type_labels": "Cell Types",
            "cell_line_labels": "Cell Lines",
            "disease_labels": "Diseases",
            "sex_labels": "Sex",
            "developmental_stage_labels": "Developmental Stages",
        }
        if field in display_name_map:
            display_name = display_name_map[field]
        else:
            display_name = field.replace("_", " ").title()
        options = []
        option_value_map = {}
        field_selected = [str(v).strip().lower() for v in selected_filters.get(field, []) if v is not None]
        for item in values:
            raw_value = item.get("value")
            count = item.get("count", 0)
            if raw_value is None:
                continue
            value = str(raw_value).strip()
            if not value:
                continue
            # Keep values with count 0 if they are currently selected, so user can deselect
            if count <= 0 and value.lower() not in field_selected:
                continue
            options.append({"label": f"{value} ({count})", "value": value})
            option_value_map[value.lower()] = value

        if not options:
            continue

        selected_values = []
        for item in selected_filters.get(field, []):
            if item is None:
                continue
            cleaned_value = str(item).strip()
            if not cleaned_value:
                continue
            mapped_value = option_value_map.get(cleaned_value.lower())
            if mapped_value and mapped_value not in selected_values:
                selected_values.append(mapped_value)

        if len(options) <= 10:
            control = dcc.Checklist(
                id={"type": "facet-filter", "field": field},
                options=options,
                value=selected_values,
                inputStyle={"marginRight": "0.5rem"},
                labelStyle={
                    "display": "block",
                    "marginBottom": "0.35rem",
                    "fontSize": "0.9rem",
                },
            )
        else:
            control = dcc.Dropdown(
                id={"type": "facet-filter", "field": field},
                options=options,
                value=selected_values,
                multi=True,
                placeholder=f"Filter by {display_name}",
                className="facet-dropdown",
                style={"fontSize": "0.9rem", "zIndex": 2000, "position": "relative"},
            )

        # Build header with icon and explanation popover
        icon_class = field_icons.get(field)
        explanation = field_explanations.get(field)
        info_icon_id = f"facet-info-{field}"

        # Build the header content matching GSEA tooltip structure exactly
        header_children = []
        if icon_class:
            header_children.append(html.I(className=f"bi {icon_class} me-2"))

        header_children.append(f"{display_name} ")

        # Add info icon and popover if explanation exists (target mode only)
        if explanation and search_mode != "datasets":
            header_children.append(
                html.Span(
                    html.I(className="bi bi-question-circle me-2"),
                    id=info_icon_id,
                    style={"cursor": "pointer", "marginLeft": "2px"},
                )
            )
            header_children.append(
                dbc.Popover(
                    [
                        dbc.PopoverHeader(explanation["header"]),
                        dbc.PopoverBody(explanation["body"]),
                    ],
                    target=info_icon_id,
                    trigger="click",
                    placement="right",
                )
            )


        controls.append(
            dbc.Card(
                [
                    dbc.CardHeader(
                        header_children,
                        className="fw-semibold",
                        style={"backgroundColor": "#f8f9fa"},
                    ),
                    dbc.CardBody(control, style={"padding": "0.75rem"}),
                ],
                className="mb-3 shadow-sm facet-filter-card",
                style={
                    "borderRadius": "12px",
                    "overflow": "visible",
                    "position": "relative",
                    "zIndex": 1,
                },
            )
        )

    if not controls:
        return _filter_placeholder("No filters available for these results.")

    return html.Div(controls, className="facet-controls")


# Layout for home page
layout = html.Div(
    [
        html.Div(
            [
                dbc.Container(
                    dbc.Row(
                        dbc.Col(
                            dbc.Card(
                                [
                                    html.Div(
                                        [
                                            html.H1(
                                                "Perturbation Catalogue",
                                                className="banner-title",
                                                style={
                                                    "fontSize": "2.5rem",
                                                    "color": "#212529",
                                                },
                                            ),
                                            html.Div(
                                                [
                                                    dbc.InputGroup(
                                                        [
                                                            dbc.Select(
                                                                id="search-mode-dropdown",
                                                                options=[
                                                                    {"label": "Targets", "value": "targets"},
                                                                    {"label": "Datasets", "value": "datasets"},
                                                                ],
                                                                value="targets",
                                                                className="form-select-lg",
                                                                style={
                                                                    "borderRadius": "8px 0 0 8px",
                                                                    "maxWidth": "140px",
                                                                    "borderRight": "none",
                                                                    "backgroundColor": "#f8f9fa",
                                                                    "fontWeight": "500",
                                                                },
                                                            ),
                                                            dbc.Input(
                                                                id="search-input",
                                                                placeholder="Search by metadata fields...",
                                                                type="text",
                                                                value="",
                                                                debounce=True,
                                                                className="form-control-lg",
                                                                style={
                                                                    "borderRadius": "0"
                                                                },
                                                            ),
                                                            dbc.Button(
                                                                [
                                                                    html.I(
                                                                        className="bi bi-search me-2"
                                                                    ),
                                                                    "Search",
                                                                ],
                                                                id="search-button",
                                                                color="primary",
                                                                n_clicks=0,
                                                                className="btn-lg",
                                                                style={
                                                                    "backgroundColor": COLORS[
                                                                        "primary"
                                                                    ],
                                                                    "borderColor": COLORS[
                                                                        "primary"
                                                                    ],
                                                                    "borderRadius": "0 8px 8px 0",
                                                                    "paddingLeft": "2rem",
                                                                    "paddingRight": "2rem",
                                                                },
                                                            ),
                                                        ],
                                                        className="search-input-group banner-search",
                                                    ),
                                                    html.Div(
                                                        [
                                                            html.Span(
                                                                "Try searching for: ",
                                                                className="text-muted me-2",
                                                                style={
                                                                    "fontSize": "0.9rem"
                                                                },
                                                            ),
                                                            dbc.Button(
                                                                "SUMO1",
                                                                id="search-example-1",
                                                                n_clicks=0,
                                                                color="link",
                                                                className="p-1 me-2",
                                                                style={
                                                                    "color": COLORS[
                                                                        "primary"
                                                                    ],
                                                                    "textDecoration": "none",
                                                                    "fontSize": "0.9rem",
                                                                    "border": f"1px solid {COLORS['primary']}",
                                                                    "borderRadius": "4px",
                                                                    "background": "transparent",
                                                                    "padding": "0.25rem 0.5rem",
                                                                    "verticalAlign": "baseline",
                                                                },
                                                            ),
                                                            html.Span(
                                                                ", ",
                                                                className="text-muted me-1",
                                                            ),
                                                            dbc.Button(
                                                                "retina",
                                                                id="search-example-2",
                                                                n_clicks=0,
                                                                color="link",
                                                                className="p-1 me-2",
                                                                style={
                                                                    "color": COLORS[
                                                                        "primary"
                                                                    ],
                                                                    "textDecoration": "none",
                                                                    "fontSize": "0.9rem",
                                                                    "border": f"1px solid {COLORS['primary']}",
                                                                    "borderRadius": "4px",
                                                                    "background": "transparent",
                                                                    "padding": "0.25rem 0.5rem",
                                                                    "verticalAlign": "baseline",
                                                                },
                                                            ),
                                                            html.Span(
                                                                ", ",
                                                                className="text-muted me-1",
                                                            ),
                                                            dbc.Button(
                                                                "acute myeloid leukemia",
                                                                id="search-example-3",
                                                                n_clicks=0,
                                                                color="link",
                                                                className="p-1",
                                                                style={
                                                                    "color": COLORS[
                                                                        "primary"
                                                                    ],
                                                                    "textDecoration": "none",
                                                                    "fontSize": "0.9rem",
                                                                    "border": f"1px solid {COLORS['primary']}",
                                                                    "borderRadius": "4px",
                                                                    "background": "transparent",
                                                                    "padding": "0.25rem 0.5rem",
                                                                    "verticalAlign": "baseline",
                                                                },
                                                            ),
                                                        ],
                                                        className="mt-2",
                                                        style={"textAlign": "left"},
                                                    ),
                                                ]
                                            ),
                                            html.P(
                                                "Perturbation Catalogue is a curated database that brings together data from various genetic perturbation experiments, including Perturb-Seq, CRISPR and MAVE screens, making it easier for researchers to study how modifying genes or proteins affects biological function across various biological and molecular contexts.",
                                                className="banner-description",
                                                style={"color": "#495057"},
                                            ),
                                        ],
                                        className="banner-content",
                                    )
                                ],
                                className="banner-card",
                                style={
                                    "backgroundColor": "#ffffff",
                                    "borderRadius": "0",
                                    "boxShadow": "0 4px 20px rgba(0, 0, 0, 0.15)",
                                    "padding": "1rem",
                                },
                            ),
                        )
                    ),
                    className="banner-container",
                )
            ],
            className="homepage-banner",
        ),
        dbc.Container(
            [
                dbc.Row(
                    dbc.Col(
                        html.Div(id="homepage-summary", className="mt-4"),
                        xs=12,
                    )
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(
                                id="facet-filters",
                                className="mt-4",
                                style={
                                    "position": "relative",
                                    "zIndex": 50,
                                    "overflow": "visible",
                                    "display": "none",
                                },
                            ),
                            xs=12,
                            sm=12,
                            md=12,
                            lg=3,
                            className="mb-4",
                            style={
                                "overflow": "visible",
                                "position": "relative",
                                "zIndex": 50,
                            },
                        ),
                        dbc.Col(
                            dcc.Loading(
                                id="search-results-loading",
                                type="circle",
                                color=COLORS["primary"],
                                children=html.Div(
                                    id="search-results",
                                    className="mt-4",
                                    children=[
                                        html.Div(
                                            id="download-metadata-container",
                                            className="mb-3",
                                            style={"display": "none"},
                                            children=[
                                                dbc.Button(
                                                    [
                                                        html.I(
                                                            className="bi bi-download me-2"
                                                        ),
                                                        "Download Metadata",
                                                    ],
                                                    id="download-metadata-btn",
                                                    color="primary",
                                                    size="sm",
                                                    style={
                                                        "backgroundColor": COLORS[
                                                            "primary"
                                                        ],
                                                        "borderColor": COLORS[
                                                            "primary"
                                                        ],
                                                        "borderRadius": "6px",
                                                    },
                                                ),
                                                dcc.Download(id="download-metadata"),
                                            ],
                                        ),
                                        html.Div(id="search-results-table"),
                                        html.Div(
                                            id="search-results-pagination",
                                            className="mt-3 pagination-bar",
                                            style={"display": "none"},
                                            children=[
                                                dbc.Button(
                                                    [
                                                        html.I(
                                                            className="bi bi-arrow-left me-1"
                                                        ),
                                                        "Previous",
                                                    ],
                                                    id="search-page-prev",
                                                    color="secondary",
                                                    size="sm",
                                                    disabled=True,
                                                    style={"borderRadius": "6px"},
                                                    className="me-3",
                                                ),
                                                html.Span(
                                                    "Page 1 of 1",
                                                    id="search-page-info",
                                                    className="fw-semibold me-3",
                                                ),
                                                dbc.Button(
                                                    [
                                                        "Next",
                                                        html.I(
                                                            className="bi bi-arrow-right ms-1"
                                                        ),
                                                    ],
                                                    id="search-page-next",
                                                    color="secondary",
                                                    size="sm",
                                                    disabled=True,
                                                    style={"borderRadius": "6px"},
                                                ),
                                            ],
                                        ),
                                    ],
                                ),
                                target_components={
                                    "search-results-table": "children",
                                    "search-results-pagination": "style",
                                },
                                delay_show=200,
                                delay_hide=100,
                            ),
                            xs=12,
                            sm=12,
                            md=12,
                            lg=9,
                            style={"position": "relative", "zIndex": 10},
                        ),
                    ],
                    className="py-4 g-4",
                ),
                dcc.Store(
                    id="search-results-store",
                    data={"results": [], "facets": {}, "query": "", "total": 0},
                ),
                dcc.Store(id="search-results-page", data=1),
            ],
            className="content-container",
        ),
    ]
)


@callback(
    Output("search-input", "value"),
    Output("search-button", "n_clicks"),
    Input("search-example-1", "n_clicks"),
    Input("search-example-2", "n_clicks"),
    Input("search-example-3", "n_clicks"),
    State("search-button", "n_clicks"),
    prevent_initial_call=True,
)
def handle_search_examples(
    example1_clicks, example2_clicks, example3_clicks, current_button_clicks
):
    """Handle clicks on search example links."""
    ctx = callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update

    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    current_clicks = current_button_clicks or 0

    if trigger_id == "search-example-1":
        return "SUMO1", current_clicks + 1
    elif trigger_id == "search-example-2":
        return "retina", current_clicks + 1
    elif trigger_id == "search-example-3":
        return "acute myeloid leukemia", current_clicks + 1

    return dash.no_update, dash.no_update


@callback(
    Output("search-results-store", "data"),
    Output("homepage-summary", "children", allow_duplicate=True),
    Output("facet-filters", "style"),
    Output("search-results-table", "children", allow_duplicate=True),
    Input("search-input", "value"),
    Input("search-button", "n_clicks"),
    Input("search-mode-dropdown", "value"),
    prevent_initial_call="initial_duplicate",
)
def update_search_results(query, _, search_mode):
    """Fetch fuzzy search results and store them."""
    search_term = (query or "").strip() if query else ""
    search_mode = search_mode or "targets"

    summary_content = dash.no_update
    filters_style = dash.no_update

    if not search_term:
        summary_data, error = get_landing_page_summary()
        summary_content = _build_summary_component(summary_data, error)
        filters_style = {
            "display": "none",
            "position": "relative",
            "zIndex": 50,
            "overflow": "visible",
        }
        return (
            {
                "results": [],
                "facets": {},
                "query": "",
                "total": 0,
                "page": 1,
                "size": SEARCH_RESULTS_PAGE_SIZE,
                "search_mode": search_mode,
            },
            summary_content,
            filters_style,
            [],  # Clear search results table
        )

    data = fetch_search_results(
        query=search_term if search_term else None,
        page=1,
        size=SEARCH_RESULTS_PAGE_SIZE,
        search_mode=search_mode,
    )
    results = data.get("results", [])
    facets = data.get("facets", {})

    store_payload = {
        "results": results,
        "facets": facets,
        "query": search_term,
        "total": data.get("total", len(results)),
        "page": data.get("page", 1),
        "size": data.get("size", SEARCH_RESULTS_PAGE_SIZE),
        "total_pages": data.get("total_pages"),
        "search_mode": search_mode,
    }

    summary_content = ""
    filters_style = {
        "position": "relative",
        "zIndex": 50,
        "overflow": "visible",
        "display": "block",
    }

    # Return dash.no_update for search-results-table so the second callback handles rendering
    return store_payload, summary_content, filters_style, dash.no_update


@callback(
    Output("search-results-table", "children", allow_duplicate=True),
    Output("search-results-pagination", "style"),
    Output("search-page-info", "children"),
    Output("search-page-prev", "disabled"),
    Output("search-page-next", "disabled"),
    Output("search-results-page", "data"),
    Output("homepage-summary", "children", allow_duplicate=True),
    Output("facet-filters", "children"),
    Output("download-metadata-container", "style"),
    Input("search-results-store", "data"),
    Input({"type": "facet-filter", "field": ALL}, "value"),
    Input("search-page-prev", "n_clicks"),
    Input("search-page-next", "n_clicks"),
    State({"type": "facet-filter", "field": ALL}, "id"),
    State("search-results-page", "data"),
    prevent_initial_call="initial_duplicate",
)
def render_filtered_results(
    store_data, selected_values, prev_clicks, next_clicks, filter_ids, current_page
):
    """Render the results table and pagination based on search data and selected filters."""
    ctx = callback_context
    trigger = ctx.triggered[0]["prop_id"] if ctx.triggered else None

    if not store_data or not store_data.get("query"):
        return (
            "",
            {"display": "none"},
            "Page 1 of 1",
            True,
            True,
            1,
            dash.no_update,
            html.Div(),
            {"display": "none"},
        )

    selected_filters = {}
    if selected_values and filter_ids:
        for values, filter_id in zip(selected_values, filter_ids):
            if values:
                cleaned = [str(v).strip() for v in values if v not in (None, "")]
                if cleaned:
                    selected_filters[filter_id["field"]] = cleaned

    query = store_data.get("query")

    if not query:
        return (
            "",
            {"display": "none"},
            "Page 1 of 1",
            True,
            True,
            1,
            dash.no_update,
            html.Div(),
            {"display": "none"},
        )

    # Determine search mode
    search_mode = store_data.get("search_mode", "targets")
    is_dataset_mode = search_mode == "datasets"

    # Always use stored results for search queries - filter client-side
    all_results = store_data.get("results", [])
    original_facets = store_data.get("facets", {})

    # Apply client-side filtering
    def matches_filters(result, filters):
        for field, selected_values in filters.items():
            if not selected_values:
                continue
            # For dataset mode, map canonical target field names to dataset field names
            actual_field = TARGET_TO_DATASET_FIELD.get(field, field) if is_dataset_mode else field
            result_value = result.get(actual_field)
            if result_value is None:
                return False
            # Handle both list and scalar values
            if isinstance(result_value, list):
                if not any(v in selected_values for v in result_value):
                    return False
            else:
                if result_value not in selected_values:
                    return False
        return True

    if selected_filters:
        filtered_results = [r for r in all_results if matches_filters(r, selected_filters)]
    else:
        filtered_results = all_results

    # Recalculate facet counts based on filtered results
    if selected_filters:
        facets = _recalculate_facets(filtered_results, original_facets, search_mode)
    else:
        facets = original_facets

    # Pagination: 10 records per page
    page_size = 10
    total_results = len(filtered_results)
    total_pages = max(1, (total_results + page_size - 1) // page_size)

    # Handle page navigation
    triggered_prop = trigger or ""
    current_page_number = current_page or 1

    if "search-results-store" in triggered_prop or "facet-filter" in triggered_prop:
        # Reset to page 1 on new search or filter change
        current_page_number = 1
    elif "search-page-prev" in triggered_prop:
        current_page_number = max(1, current_page_number - 1)
    elif "search-page-next" in triggered_prop:
        current_page_number = min(total_pages, current_page_number + 1)

    # Ensure page is within bounds
    if current_page_number > total_pages:
        current_page_number = total_pages
    if current_page_number < 1:
        current_page_number = 1

    # Slice results for current page
    start_idx = (current_page_number - 1) * page_size
    end_idx = start_idx + page_size
    results = filtered_results[start_idx:end_idx]

    if not results and total_results == 0:
        filters_children = _build_filter_controls(facets, selected_filters, FACET_FIELDS, search_mode)
        no_results_msg = "No datasets found. Try another search term." if is_dataset_mode else "No results found. Try another target name."
        return (
            dbc.Alert(
                [
                    html.I(className="bi bi-info-circle me-2"),
                    no_results_msg,
                ],
                color="info",
                className="shadow-sm",
                style={"borderRadius": "10px"},
            ),
            {"display": "none"},
            f"Page 1 of {total_pages} ({total_results} results)",
            True,
            True,
            current_page_number,
            "",
            filters_children,
            {"display": "none"},
        )

    if is_dataset_mode:
        table = _render_dataset_search_results(results)
    else:
        table = _render_search_results(results)
    filters_children = _build_filter_controls(facets, selected_filters, FACET_FIELDS, search_mode)

    pagination_style = {"display": "flex"} if total_pages > 1 else {"display": "none"}
    page_info = f"Page {current_page_number} of {total_pages} ({total_results} results)"
    prev_disabled = current_page_number <= 1
    next_disabled = current_page_number >= total_pages

    return (
        table,
        pagination_style,
        page_info,
        prev_disabled,
        next_disabled,
        current_page_number,
        "",
        filters_children,
        {"display": "flex", "justifyContent": "flex-end"},
    )


@callback(
    Output("download-metadata", "data"),
    Input("download-metadata-btn", "n_clicks"),
    State("search-results-store", "data"),
    State({"type": "facet-filter", "field": ALL}, "value"),
    State({"type": "facet-filter", "field": ALL}, "id"),
    prevent_initial_call=True,
)
def download_metadata(n_clicks, store_data, selected_values, filter_ids):
    """Download all search results as CSV when the button is clicked."""
    if not n_clicks or not store_data or not store_data.get("query"):
        return dash.no_update

    query = store_data.get("query")
    search_mode = store_data.get("search_mode", "targets")

    # Build filters from current selection
    selected_filters = {}
    if selected_values and filter_ids:
        for values, filter_id in zip(selected_values, filter_ids):
            if values:
                cleaned = [str(v).strip() for v in values if v not in (None, "")]
                if cleaned:
                    selected_filters[filter_id["field"]] = cleaned

    # Fetch all results by paginating through the API
    all_results = fetch_all_search_results(
        query=query,
        filters=selected_filters or None,
        search_mode=search_mode,
    )

    if not all_results:
        return dash.no_update

    # Define columns to include in the CSV based on search mode
    if search_mode == "datasets":
        columns = [
            "dataset_id",
            "study_title",
            "study_year",
            "data_modalities",
        ]
    else:
        columns = [
            "perturbed_target_symbol",
            "n_experiments",
            "n_sig_perturb_pairs_up",
            "n_sig_perturb_pairs_down",
            "n_sig_crispr",
            "n_mave",
            "top_gsea_terms",
            "data_modalities",
            "tissues_tested",
            "cell_types_tested",
            "cell_lines_tested",
            "sex_tested",
            "developmental_stages_tested",
            "diseases_tested",
            "license",
        ]

    # Build CSV content
    csv_lines = [",".join(columns)]
    for record in all_results:
        row_values = []
        for col in columns:
            value = record.get(col, "")
            if isinstance(value, list):
                value = "; ".join(str(v) for v in value)
            elif value is None:
                value = ""
            else:
                value = str(value)
            # Escape quotes and wrap in quotes if contains comma or quote
            if "," in value or '"' in value or "\n" in value:
                value = '"' + value.replace('"', '""') + '"'
            row_values.append(value)
        csv_lines.append(",".join(row_values))

    csv_content = "\n".join(csv_lines)

    # Generate filename with query and search mode
    safe_query = "".join(c if c.isalnum() or c in "-_" else "_" for c in query[:30])
    filename = f"perturbation_catalogue_{search_mode}_{safe_query}.csv"

    return dict(content=csv_content, filename=filename, type="text/csv")
