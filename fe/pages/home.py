"""Home page for Perturbation Search"""
import dash
from dash import dcc, html, Input, Output, State, callback_context, ALL, callback
import dash_bootstrap_components as dbc
import json
from utils import (
    COLORS, FACET_FIELDS, fetch_search_results, results_store,
    format_value, format_number, create_chart_component,
    format_rnaseq_stats, format_crispr_stats, format_mave_stats
)

# Register this page with Dash Pages
dash.register_page(__name__, path="/")

# Layout for home page
layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H1("Perturbation Search", className="text-center mb-4",
                   style={"color": COLORS["primary"], "fontWeight": "700", "letterSpacing": "-0.5px"}),
            html.P("Search and explore perturbation data across multiple datasets",
                  className="text-center text-muted mb-5", style={"fontSize": "1.1rem"})
        ], width=12)
    ]),
    dbc.Row([
        dbc.Col([
            html.Div([
                html.Div([
                    dbc.InputGroup([
                        dbc.Input(
                            id="search-input",
                            placeholder="Search by target symbol (e.g., BRCA1, TP53)...",
                            type="text",
                            value="",
                            debounce=True,
                            className="form-control-lg",
                            style={"borderRadius": "8px 0 0 8px"}
                        ),
                        dbc.Button(
                            [html.I(className="bi bi-search me-2"), "Search"],
                            id="search-button",
                            color="primary",
                            n_clicks=0,
                            className="btn-lg",
                            style={
                                "backgroundColor": COLORS["primary"],
                                "borderColor": COLORS["primary"],
                                "borderRadius": "0 8px 8px 0",
                                "paddingLeft": "2rem",
                                "paddingRight": "2rem"
                            }
                        )
                    ], className="search-input-group", style={"maxWidth": "600px", "margin": "0 auto"})
                ], className="mb-4"),
                html.Div(id="pagination-container", children=[
                    dbc.Button([
                        html.I(className="bi bi-arrow-left me-1"),
                        "Previous"
                    ], id="prev-page", color="secondary",
                              disabled=True, className="me-2", size="sm"),
                    html.Span("", id="pagination-info", className="mx-3"),
                    dbc.Button([
                        "Next",
                        html.I(className="bi bi-arrow-right ms-1")
                    ], id="next-page", color="secondary",
                              disabled=True, className="ms-2", size="sm")
                ], className="pagination-container", style={"display": "none", "margin": "0 auto 2rem"})
            ], style={"textAlign": "center"})
        ], width=12)
    ], className="mb-5"),
    # Charts section - filters displayed as charts
    dbc.Row([
        dbc.Col([
            html.H4([
                html.I(className="bi bi-bar-chart me-2"),
                "Filters"
            ], className="mb-2", style={"color": COLORS["primary"], "fontWeight": "600", "fontSize": "1.25rem"}),
            html.Div(id="chart-container"),
            dbc.Button(
                "Clear All Filters",
                id="clear-filters",
                color="outline-secondary",
                className="mt-2 w-100",
                size="sm",
                style={"borderRadius": "6px"}
            )
        ], width=12, className="mb-3")
    ]),
    # Results table section - only shown when filters/search are active
    dbc.Row([
        dbc.Col([
            dcc.Loading(
                id="results-loading",
                type="default",
                children=[
                    html.Div(id="results-container")
                ]
            )
        ], width=12)
    ]),

    # Store for current filters and search state
    dcc.Store(id="current-filters", data={}),
    dcc.Store(id="current-page", data=1),
    dcc.Store(id="current-query", data=""),
], fluid=True)

# Callbacks for home page
@callback(
    Output("chart-container", "children"),
    Output("current-filters", "data"),
    Input("search-input", "value"),
    Input({"type": "filter-dropdown", "field": ALL}, "value"),
    Input("clear-filters", "n_clicks"),
    State("current-filters", "data"),
    prevent_initial_call=False
)
def update_charts(search_query, filter_values, clear_clicks, current_filters):
    """Update chart components based on current search results"""
    ctx = callback_context

    # Initialize filters if not exists
    if not current_filters:
        current_filters = {field: [] for field in FACET_FIELDS}

    # Handle clear filters button
    if ctx.triggered and "clear-filters" in ctx.triggered[0]["prop_id"]:
        current_filters = {field: [] for field in FACET_FIELDS}
        # Fetch fresh facets (even with empty query to get all facets)
        # Use size=1 to minimize data transfer since we only need facets
        data = fetch_search_results(query=search_query or None, filters=current_filters, page=1, size=1)
    else:
        # Update filters from dropdown values
        if ctx.triggered:
            for prop in ctx.triggered:
                if "filter-dropdown" in prop["prop_id"]:
                    try:
                        prop_id = json.loads(prop["prop_id"].split(".")[0])
                        field = prop_id["field"]
                        value = prop["value"] or []
                        current_filters[field] = value
                    except:
                        pass

        # Fetch data with current filters (even with empty query to get all facets)
        # Use size=1 to minimize data transfer since we only need facets
        data = fetch_search_results(query=search_query or None, filters=current_filters, page=1, size=1)

    facets = data.get("facets", {})

    # Create chart components - arrange in a grid
    chart_components = []
    for field in FACET_FIELDS:
        facet_values = facets.get(field, [])
        selected = current_filters.get(field, [])
        chart_components.append(
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        create_chart_component(field, facet_values, selected)
                    ], style={"padding": "0.75rem"})
                ], className="h-100", style={"borderRadius": "8px", "boxShadow": "0 1px 4px rgba(0,0,0,0.08)"})
            ], width={"size": 4, "lg": 4, "md": 6, "sm": 12, "xs": 12}, className="mb-2")
        )

    # Wrap charts in a row container for responsive grid layout
    # Bootstrap will automatically wrap columns when they exceed 12 units
    chart_container = dbc.Row(chart_components, className="g-2")

    return chart_container, current_filters

@callback(
    Output("results-container", "children"),
    Output("pagination-container", "children"),
    Output("current-page", "data"),
    Input("search-input", "value"),
    Input("search-button", "n_clicks"),
    Input("current-filters", "data"),
    Input("current-page", "data"),
    State("current-query", "data"),
    prevent_initial_call=False
)
def update_results(search_query, search_clicks, filters, page, stored_query):
    """Update search results based on query and filters - show as table"""
    ctx = callback_context

    # Check if filters or search are active
    has_search = search_query and search_query.strip()
    has_filters = filters and any(filters.get(field, []) for field in FACET_FIELDS)
    should_show_table = has_search or has_filters

    # Determine if we should reset to page 1
    if ctx.triggered:
        trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if trigger_id in ["search-input", "search-button"] or "current-filters" in trigger_id:
            if search_query != stored_query or filters != {}:
                page = 1

    # Only fetch data if we should show the table
    if not should_show_table:
        return html.Div([
            dbc.Alert(
                [
                    html.I(className="bi bi-info-circle me-2"),
                    "Use the search box or select filters above to see results."
                ],
                color="info",
                className="text-center",
                style={"borderRadius": "8px", "marginTop": "2rem"}
            )
        ]), html.Div(), page

    # Fetch data - increase size for table view
    data = fetch_search_results(query=search_query, filters=filters, page=page, size=20)
    results = data.get("results", [])
    total = data.get("total", 0)
    total_pages = data.get("total_pages", 0)

    # Store results by target symbol for navigation to details page
    results_store.clear()
    results_store.update({r.get("perturbed_target_symbol"): r for r in results})

    # Create results table
    if results:
        table_rows = [
            html.Tr([
                html.Th("Perturbed Target Symbol", style={"padding": "12px", "border": "1px solid #ddd",
                       "backgroundColor": COLORS["primary"], "color": "white", "fontWeight": "600"}),
                html.Th("Total Experiments", style={"padding": "12px", "border": "1px solid #ddd",
                       "backgroundColor": COLORS["primary"], "color": "white", "fontWeight": "600"}),
                html.Th("RNA-seq Statistics", style={"padding": "12px", "border": "1px solid #ddd",
                       "backgroundColor": COLORS["primary"], "color": "white", "fontWeight": "600"}),
                html.Th("CRISPR Screen Statistics", style={"padding": "12px", "border": "1px solid #ddd",
                       "backgroundColor": COLORS["primary"], "color": "white", "fontWeight": "600"}),
                html.Th("MAVE Statistics", style={"padding": "12px", "border": "1px solid #ddd",
                       "backgroundColor": COLORS["primary"], "color": "white", "fontWeight": "600"})
            ])
        ]

        for result in results:
            perturbed_target = result.get("perturbed_target_symbol", "N/A")
            n_experiments = format_value(result.get("n_experiments", "N/A"))
            rnaseq_stats = format_rnaseq_stats(result)
            crispr_stats = format_crispr_stats(result)
            mave_stats = format_mave_stats(result)

            # Create link for perturbed_target_symbol
            target_link = dcc.Link(
                perturbed_target,
                href=f"/target/{perturbed_target}",
                style={
                    "color": COLORS["primary"],
                    "textDecoration": "none",
                    "fontWeight": "500"
                },
                className="target-link"
            )

            table_rows.append(
                html.Tr([
                    html.Td(target_link, style={"padding": "12px", "border": "1px solid #ddd"}),
                    html.Td(n_experiments, style={"padding": "12px", "border": "1px solid #ddd"}),
                    html.Td(rnaseq_stats, style={"padding": "12px", "border": "1px solid #ddd"}),
                    html.Td(crispr_stats, style={"padding": "12px", "border": "1px solid #ddd"}),
                    html.Td(mave_stats, style={"padding": "12px", "border": "1px solid #ddd"})
                ], style={"cursor": "pointer"})
            )

        table = dbc.Table([
            html.Thead(table_rows[0]),
            html.Tbody(table_rows[1:])
        ], bordered=True, hover=True, responsive=True, striped=True,
           className="table align-middle", style={"fontSize": "0.9rem", "marginTop": "1rem"})

        results_display = html.Div([
            html.H5(f"Results ({total} total)", className="mb-3", style={"color": COLORS["primary"], "fontWeight": "600"}),
            table
        ])
    else:
        results_display = html.Div([
            dbc.Alert(
                [
                    html.I(className="bi bi-info-circle me-2"),
                    "No results found. Try adjusting your search terms or filters."
                ],
                color="info",
                className="text-center",
                style={"borderRadius": "8px", "marginTop": "2rem"}
            )
        ])

    # Update pagination
    prev_disabled = page == 1 or total_pages == 0
    next_disabled = page >= total_pages or total_pages == 0

    # Page info
    if total > 0:
        pagination_info = f"Page {page} of {total_pages} ({total} total results)"
    else:
        pagination_info = "No results"

    # Create updated pagination container with buttons
    if total > 0:
        pagination = html.Div([
            dbc.Button([
                html.I(className="bi bi-arrow-left me-1"),
                "Previous"
            ], id="prev-page", color="secondary",
                      disabled=prev_disabled, className="me-2", size="sm",
                      style={"borderRadius": "6px"}),
            html.Span(pagination_info, id="pagination-info", className="mx-3",
                     style={"fontWeight": "500", "color": "#495057"}),
            dbc.Button([
                "Next",
                html.I(className="bi bi-arrow-right ms-1")
            ], id="next-page", color="secondary",
                      disabled=next_disabled, className="ms-2", size="sm",
                      style={"borderRadius": "6px"})
        ], className="pagination-container", style={"display": "inline-flex", "margin": "2rem auto"})
    else:
        pagination = html.Div()

    return results_display, pagination, page

@callback(
    Output("current-page", "data", allow_duplicate=True),
    Input("prev-page", "n_clicks"),
    Input("next-page", "n_clicks"),
    State("current-page", "data"),
    prevent_initial_call=True
)
def update_page(prev_clicks, next_clicks, current_page):
    """Handle pagination navigation"""
    ctx = callback_context
    if not ctx.triggered:
        return current_page

    trigger_id = ctx.triggered[0]["prop_id"]
    if "prev-page" in trigger_id and current_page > 1:
        return current_page - 1
    elif "next-page" in trigger_id:
        return current_page + 1

    return current_page

@callback(
    Output("current-query", "data"),
    Input("search-input", "value"),
    prevent_initial_call=False
)
def update_stored_query(query):
    """Update stored query value"""
    return query or ""