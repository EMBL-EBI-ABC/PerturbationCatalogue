"""Home page for Perturbation Search"""

import dash
from dash import dcc, html, Input, Output, State, callback, callback_context
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from urllib.parse import quote
from utils import (
    COLORS,
    get_landing_page_summary,
)

# Register this page with Dash Pages
dash.register_page(__name__, path="/")


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


# Layout for home page
layout = html.Div(
    [
        dcc.Location(id="home-redirect", refresh=True),
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
                                                            dbc.InputGroupText(
                                                                [
                                                                    html.Span(
                                                                        html.I(
                                                                            className="bi bi-question-circle"
                                                                        ),
                                                                        id="search-mode-info-icon",
                                                                        style={
                                                                            "cursor": "pointer",
                                                                            "color": "#6c757d",
                                                                        },
                                                                    ),
                                                                    dbc.Popover(
                                                                        [
                                                                            dbc.PopoverHeader(
                                                                                "Search Mode"
                                                                            ),
                                                                            dbc.PopoverBody(
                                                                                [
                                                                                    html.Strong(
                                                                                        "Targets:"
                                                                                    ),
                                                                                    " Search for gene targets across all experiments. Results show aggregated data for each target gene, including the number of datasets, significant perturbation effects, and pathway enrichment across Perturb-seq, CRISPR, and MAVE experiments.",
                                                                                    html.Br(),
                                                                                    html.Br(),
                                                                                    html.Strong(
                                                                                        "Datasets:"
                                                                                    ),
                                                                                    " Search for datasets by metadata fields such as dataset ID, study title, tissue, cell type, disease, and more. Results show individual experiments with their study details, publication year, and data modality.",
                                                                                ]
                                                                            ),
                                                                        ],
                                                                        target="search-mode-info-icon",
                                                                        trigger="click",
                                                                        placement="bottom",
                                                                    ),
                                                                ],
                                                                style={
                                                                    "backgroundColor": "#f8f9fa",
                                                                    "borderRadius": "8px 0 0 8px",
                                                                    "borderRight": "none",
                                                                },
                                                            ),
                                                            dbc.Select(
                                                                id="search-mode-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": "Targets",
                                                                        "value": "targets",
                                                                    },
                                                                    {
                                                                        "label": "Datasets",
                                                                        "value": "datasets",
                                                                    },
                                                                ],
                                                                value="targets",
                                                                className="form-select-lg",
                                                                style={
                                                                    "borderRadius": "0",
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
            ],
            className="content-container",
        ),
    ]
)


@callback(
    Output("home-redirect", "href"),
    Input("search-input", "value"),
    Input("search-button", "n_clicks"),
    Input("search-example-1", "n_clicks"),
    Input("search-example-2", "n_clicks"),
    Input("search-example-3", "n_clicks"),
    State("search-mode-dropdown", "value"),
    prevent_initial_call=True,
)
def redirect_search(query, search_clicks, ex1, ex2, ex3, mode):
    """Redirect search to the targets or datasets browse page."""
    ctx = callback_context
    if not ctx.triggered:
        return dash.no_update

    trigger = ctx.triggered[0]["prop_id"].split(".")[0]

    if trigger == "search-example-1":
        query = "SUMO1"
    elif trigger == "search-example-2":
        query = "retina"
    elif trigger == "search-example-3":
        query = "acute myeloid leukemia"

    if not query or not query.strip():
        return dash.no_update

    mode = mode or "targets"
    page = "targets" if mode == "targets" else "datasets"
    return f"/perturbation-catalogue/{page}?q={quote(query.strip())}"


@callback(
    Output("homepage-summary", "children"),
    Input("home-redirect", "href"),
    prevent_initial_call=False,
)
def load_summary(_):
    """Load the summary dashboard on page mount."""
    summary_data, error = get_landing_page_summary()
    return _build_summary_component(summary_data, error)
