"""Home page for Perturbation Search"""

import math
import dash
from dash import dcc, html, Input, Output, State, callback, ALL, callback_context
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from utils import (
    COLORS,
    DATA_MODALITIES_COLOURS,
    FACET_FIELDS,
    fetch_search_results,
    get_landing_page_summary,
    results_store,
    format_value,
)

SEARCH_RESULTS_PAGE_SIZE = 15


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
                html.Th("Target name", className="fw-semibold"),
                html.Th(
                    "Total number of experiments", className="fw-semibold text-center"
                ),
                html.Th(
                    "Significant Perturb-Seq gene pairs",
                    className="fw-semibold text-center",
                ),
                html.Th(
                    "Number of Significant CRISPR experiments",
                    className="fw-semibold text-center",
                ),
                html.Th(
                    "Number of MAVE experiments", className="fw-semibold text-center"
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

        n_experiments = _format_count(record.get("n_experiments"))
        n_sig_up = _format_count(record.get("n_sig_perturb_pairs_up"))
        n_sig_down = _format_count(record.get("n_sig_perturb_pairs_down"))
        n_sig_crispr = _format_count(record.get("n_sig_crispr"))
        n_mave = _format_count(record.get("n_mave"))
        data_modalities = record.get("data_modalities") or []

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
                            href=f"/perturbation-catalogue/target/{symbol}",
                            className="text-decoration-none fw-semibold",
                            style={"color": COLORS["primary"]},
                        )
                    ),
                    html.Td(n_experiments, className="text-center"),
                    html.Td(
                        [
                            html.Div(
                                [
                                    html.Span(
                                        [
                                            html.Span(
                                                "↑",
                                                style={
                                                    "color": "#1f9d55",
                                                    "marginRight": "0.25rem",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(
                                                n_sig_up,
                                                style={
                                                    "color": "#1f9d55",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                        ],
                                        className="d-inline-flex align-items-center me-3",
                                    ),
                                    html.Span(
                                        [
                                            html.Span(
                                                "↓",
                                                style={
                                                    "color": "#c53030",
                                                    "marginRight": "0.25rem",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(
                                                n_sig_down,
                                                style={
                                                    "color": "#c53030",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                        ],
                                        className="d-inline-flex align-items-center",
                                    ),
                                ],
                                className="d-flex justify-content-center flex-wrap",
                            )
                        ],
                        className="text-center",
                    ),
                    html.Td(n_sig_crispr, className="text-center"),
                    html.Td(n_mave, className="text-center"),
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


def _build_filter_controls(facets, selected_filters=None):
    """Build filter controls for facet fields."""
    if not facets:
        return _filter_placeholder()

    normalized_selected_filters = {}
    if selected_filters:
        for field, values in selected_filters.items():
            normalized_selected_filters[field] = [
                str(v).strip() for v in values if v not in (None, "")
            ]
    else:
        normalized_selected_filters = {}

    # Icon mapping for facet fields
    field_icons = {
        "license": "bi-award-fill",
        "data_modalities": "bi-database",
        "tissues_tested": "bi-universal-access-circle",
        "cell_types_tested": "bi-puzzle",
        "cell_lines_tested": "bi-puzzle-fill",
        "diseases_tested": "bi-virus2",
        "sex_tested": "bi-gender-ambiguous",
        "developmental_stages_tested": "bi-graph-up-arrow",
    }

    controls = []
    for field in FACET_FIELDS:
        values = facets.get(field, [])
        if not values:
            continue

        # Custom display names for specific fields
        if field == "license":
            display_name = "License"
        else:
            display_name = field.replace("_", " ").title()
        options = []
        option_value_map = {}
        for item in values:
            raw_value = item.get("value")
            count = item.get("count", 0)
            if raw_value is None or count <= 0:
                continue
            value = str(raw_value).strip()
            if not value:
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

        # Build header with icon
        icon_class = field_icons.get(field)
        if icon_class:
            header_content = html.Div(
                [
                    html.I(
                        className=f"bi {icon_class} me-2",
                    ),
                    html.Span(display_name),
                ],
                className="d-flex align-items-center",
            )
        else:
            header_content = display_name

        controls.append(
            dbc.Card(
                [
                    dbc.CardHeader(
                        header_content,
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
                                                            dbc.Input(
                                                                id="search-input",
                                                                placeholder="Search by target name...",
                                                                type="text",
                                                                value="",
                                                                debounce=True,
                                                                className="form-control-lg",
                                                                style={
                                                                    "borderRadius": "8px 0 0 8px"
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
                                                    )
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
                                type="default",
                                children=html.Div(
                                    id="search-results",
                                    className="mt-4",
                                    children=[
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
    Output("search-results-store", "data"),
    Output("homepage-summary", "children", allow_duplicate=True),
    Output("facet-filters", "style"),
    Input("search-input", "value"),
    Input("search-button", "n_clicks"),
    prevent_initial_call="initial_duplicate",
)
def update_search_results(query, _):
    """Fetch fuzzy search results and store them."""
    search_term = (query or "").strip() if query else ""

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
            },
            summary_content,
            filters_style,
        )

    data = fetch_search_results(
        query=search_term if search_term else None,
        page=1,
        size=SEARCH_RESULTS_PAGE_SIZE,
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
    }

    summary_content = ""
    filters_style = {
        "position": "relative",
        "zIndex": 50,
        "overflow": "visible",
        "display": "block",
    }

    return store_payload, summary_content, filters_style


@callback(
    Output("search-results-table", "children"),
    Output("search-results-pagination", "style"),
    Output("search-page-info", "children"),
    Output("search-page-prev", "disabled"),
    Output("search-page-next", "disabled"),
    Output("search-results-page", "data"),
    Output("homepage-summary", "children", allow_duplicate=True),
    Output("facet-filters", "children"),
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
        )

    triggered_prop = trigger or ""
    requested_page = current_page or 1

    if (
        "search-results-store" in triggered_prop
        or not triggered_prop
        or "facet-filter" in triggered_prop
    ):
        requested_page = 1
    elif "search-page-prev" in triggered_prop:
        requested_page = max(1, requested_page - 1)
    elif "search-page-next" in triggered_prop:
        requested_page = requested_page + 1

    api_response = None

    if "search-results-store" in triggered_prop or not triggered_prop:
        api_response = {
            "results": store_data.get("results", []),
            "facets": store_data.get("facets", {}),
            "total": store_data.get("total", 0),
            "page": store_data.get("page", 1),
            "size": store_data.get("size", SEARCH_RESULTS_PAGE_SIZE),
            "total_pages": store_data.get("total_pages"),
        }
    else:
        api_response = fetch_search_results(
            query=query,
            filters=selected_filters or None,
            page=requested_page,
            size=SEARCH_RESULTS_PAGE_SIZE,
        )

    total_results = api_response.get("total", 0)
    page_size = (
        api_response.get("size", SEARCH_RESULTS_PAGE_SIZE) or SEARCH_RESULTS_PAGE_SIZE
    )
    total_pages = api_response.get("total_pages")
    if not total_pages:
        total_pages = max(1, math.ceil(total_results / page_size)) if page_size else 1

    current_page_number = api_response.get("page", requested_page)
    if current_page_number < 1:
        current_page_number = 1

    if current_page_number > total_pages and total_pages > 0:
        current_page_number = total_pages
        api_response = fetch_search_results(
            query=query,
            filters=selected_filters or None,
            page=current_page_number,
            size=SEARCH_RESULTS_PAGE_SIZE,
        )
        total_results = api_response.get("total", 0)
        page_size = (
            api_response.get("size", SEARCH_RESULTS_PAGE_SIZE)
            or SEARCH_RESULTS_PAGE_SIZE
        )
        total_pages = api_response.get("total_pages") or (
            max(1, math.ceil(total_results / page_size)) if page_size else 1
        )

    results = api_response.get("results", [])
    facets = api_response.get("facets", {})

    if not results:
        filters_children = _build_filter_controls(facets, selected_filters)
        return (
            dbc.Alert(
                [
                    html.I(className="bi bi-info-circle me-2"),
                    "No results found. Try another target name.",
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
        )

    table = _render_search_results(results)
    filters_children = _build_filter_controls(facets, selected_filters)

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
    )
