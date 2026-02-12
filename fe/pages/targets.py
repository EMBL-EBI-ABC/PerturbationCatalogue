"""Browse targets page for Perturbation Catalogue."""

import dash
from dash import dcc, html, Input, Output, State, callback, ALL, callback_context
import dash_bootstrap_components as dbc
from urllib.parse import parse_qs
from utils import (
    COLORS,
    FACET_FIELDS,
    fetch_search_results,
    fetch_all_search_results,
)
from components.search_table import (
    render_targets_table,
    filter_placeholder,
    build_filter_controls,
)

dash.register_page(__name__, path="/targets", name="Targets", title="Browse Targets")

PAGE_SIZE = 10


layout = html.Div(
    [
        dcc.Location(id="targets-url", refresh=False),
        dbc.Container(
            [
                html.Div(
                    [
                        html.H2("Browse Targets", className="mb-3 mt-4"),
                        html.P(
                            "Explore all gene targets across perturbation experiments. Use the search bar to filter by gene name or metadata, and the facets to narrow down results.",
                            className="text-muted mb-3",
                        ),
                        dbc.InputGroup(
                            [
                                dbc.Input(
                                    id="targets-search-input",
                                    placeholder="Search targets by gene name, tissue, disease...",
                                    type="text",
                                    value="",
                                    debounce=True,
                                    className="form-control-lg",
                                    style={"borderRadius": "8px 0 0 8px"},
                                ),
                                dbc.Button(
                                    [
                                        html.I(className="bi bi-search me-2"),
                                        "Search",
                                    ],
                                    id="targets-search-button",
                                    color="primary",
                                    n_clicks=0,
                                    className="btn-lg",
                                    style={
                                        "backgroundColor": COLORS["primary"],
                                        "borderColor": COLORS["primary"],
                                        "borderRadius": "0 8px 8px 0",
                                        "paddingLeft": "2rem",
                                        "paddingRight": "2rem",
                                    },
                                ),
                            ],
                            className="mb-4",
                            style={"maxWidth": "700px"},
                        ),
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(
                                id="targets-facets",
                                className="mt-2",
                                style={
                                    "position": "relative",
                                    "zIndex": 50,
                                    "overflow": "visible",
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
                                id="targets-results-loading",
                                type="circle",
                                color=COLORS["primary"],
                                children=html.Div(
                                    id="targets-results",
                                    className="mt-2",
                                    children=[
                                        html.Div(
                                            id="targets-download-container",
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
                                                    id="targets-download-btn",
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
                                                dcc.Download(id="targets-download"),
                                            ],
                                        ),
                                        html.Div(id="targets-results-table"),
                                        html.Div(
                                            id="targets-pagination",
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
                                                    id="targets-page-prev",
                                                    color="secondary",
                                                    size="sm",
                                                    disabled=True,
                                                    style={"borderRadius": "6px"},
                                                    className="me-3",
                                                ),
                                                html.Span(
                                                    "Page 1 of 1",
                                                    id="targets-page-info",
                                                    className="fw-semibold me-3",
                                                ),
                                                dbc.Button(
                                                    [
                                                        "Next",
                                                        html.I(
                                                            className="bi bi-arrow-right ms-1"
                                                        ),
                                                    ],
                                                    id="targets-page-next",
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
                                    "targets-results-table": "children",
                                    "targets-pagination": "style",
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
                    className="py-2 g-4",
                ),
                dcc.Store(
                    id="targets-results-store",
                    data={
                        "results": [],
                        "facets": {},
                        "query": "",
                        "total": 0,
                        "total_pages": 0,
                        "current_page": 1,
                    },
                ),
                dcc.Store(id="targets-page-store", data=1),
            ],
            className="content-container",
        ),
    ]
)


@callback(
    Output("targets-search-input", "value"),
    Input("targets-url", "search"),
)
def populate_search_from_url(search):
    """Populate search input from URL query parameter on initial load."""
    if search:
        params = parse_qs(search.lstrip("?"))
        query = params.get("q", [None])[0]
        if query:
            return query[:500]
    return ""


@callback(
    Output("targets-results-store", "data"),
    Input("targets-search-input", "value"),
    Input("targets-search-button", "n_clicks"),
)
def load_targets_data(query, _n_clicks):
    """Fetch target search results from backend."""
    search_term = (query or "").strip() if query else ""

    data = fetch_search_results(
        query=search_term or None,
        page=1,
        size=PAGE_SIZE,
        search_mode="targets",
    )
    return {
        "results": data.get("results", []),
        "facets": data.get("facets", {}),
        "query": search_term,
        "total": data.get("total", 0),
        "total_pages": data.get("total_pages", 0),
        "current_page": 1,
    }


@callback(
    Output("targets-results-table", "children"),
    Output("targets-pagination", "style"),
    Output("targets-page-info", "children"),
    Output("targets-page-prev", "disabled"),
    Output("targets-page-next", "disabled"),
    Output("targets-page-store", "data"),
    Output("targets-facets", "children"),
    Output("targets-download-container", "style"),
    Input("targets-results-store", "data"),
    Input({"type": "targets-facet-filter", "field": ALL}, "value"),
    Input("targets-page-prev", "n_clicks"),
    Input("targets-page-next", "n_clicks"),
    State({"type": "targets-facet-filter", "field": ALL}, "id"),
    State("targets-page-store", "data"),
)
def render_targets_results(
    store_data, selected_values, prev_clicks, next_clicks, filter_ids, current_page
):
    """Render target results with server-side or client-side filtering and pagination."""
    ctx = callback_context
    trigger = ctx.triggered[0]["prop_id"] if ctx.triggered else None

    if not store_data or not store_data.get("results"):
        return (
            filter_placeholder("No targets found. Try a search query."),
            {"display": "none"},
            "Page 1 of 1",
            True,
            True,
            1,
            filter_placeholder("Search to see available filters."),
            {"display": "none"},
        )

    # Build selected filters dict
    selected_filters = {}
    if selected_values and filter_ids:
        for values, filter_id in zip(selected_values, filter_ids):
            if values:
                cleaned = [str(v).strip() for v in values if v not in (None, "")]
                if cleaned:
                    selected_filters[filter_id["field"]] = cleaned

    triggered_prop = trigger or ""

    return _render_server_side(
        store_data, selected_filters, triggered_prop, current_page
    )


def _render_server_side(store_data, selected_filters, triggered_prop, current_page):
    """Handle rendering with server-side pagination and filtering."""
    current_page_number = store_data.get("current_page", 1) or 1
    query = store_data.get("query") or None

    if "targets-results-store" in triggered_prop:
        # New search or initial load
        if selected_filters:
            # Active filters: re-fetch with both query and filters
            current_page_number = 1
            data = fetch_search_results(
                query=query,
                filters=selected_filters,
                page=1,
                size=PAGE_SIZE,
                search_mode="targets",
            )
            results = data.get("results", [])
            facets = data.get("facets", {})
            total = data.get("total", 0)
            total_pages = data.get("total_pages", 0)
        else:
            # No filters: use store data directly
            results = store_data.get("results", [])
            facets = store_data.get("facets", {})
            total = store_data.get("total", 0)
            total_pages = store_data.get("total_pages", 0)
            current_page_number = store_data.get("current_page", 1)
    elif "targets-page-prev" in triggered_prop:
        current_page_number = max(1, (current_page or 1) - 1)
        data = fetch_search_results(
            query=query,
            filters=selected_filters or None,
            page=current_page_number,
            size=PAGE_SIZE,
            search_mode="targets",
        )
        results = data.get("results", [])
        facets = data.get("facets", {})
        total = data.get("total", 0)
        total_pages = data.get("total_pages", 0)
    elif "targets-page-next" in triggered_prop:
        old_total_pages = store_data.get("total_pages", 1)
        current_page_number = min(old_total_pages, (current_page or 1) + 1)
        data = fetch_search_results(
            query=query,
            filters=selected_filters or None,
            page=current_page_number,
            size=PAGE_SIZE,
            search_mode="targets",
        )
        results = data.get("results", [])
        facets = data.get("facets", {})
        total = data.get("total", 0)
        total_pages = data.get("total_pages", 0)
    elif "facet-filter" in triggered_prop:
        # Facet changed: reset to page 1 with new filters
        current_page_number = 1
        data = fetch_search_results(
            query=query,
            filters=selected_filters or None,
            page=1,
            size=PAGE_SIZE,
            search_mode="targets",
        )
        results = data.get("results", [])
        facets = data.get("facets", {})
        total = data.get("total", 0)
        total_pages = data.get("total_pages", 0)
    else:
        # Fallback: use store data
        results = store_data.get("results", [])
        facets = store_data.get("facets", {})
        total = store_data.get("total", 0)
        total_pages = store_data.get("total_pages", 0)
        current_page_number = store_data.get("current_page", 1)

    total_pages = max(1, total_pages)
    if current_page_number > total_pages:
        current_page_number = total_pages
    if current_page_number < 1:
        current_page_number = 1

    if not results and total == 0:
        filters_children = build_filter_controls(
            facets, selected_filters, FACET_FIELDS, "targets", id_prefix="targets"
        )
        return (
            dbc.Alert(
                [
                    html.I(className="bi bi-info-circle me-2"),
                    "No targets found matching your filters.",
                ],
                color="info",
                className="shadow-sm",
                style={"borderRadius": "10px"},
            ),
            {"display": "none"},
            f"Page 1 of 1 ({total} results)",
            True,
            True,
            current_page_number,
            filters_children,
            {"display": "none"},
        )

    table = render_targets_table(results)
    filters_children = build_filter_controls(
        facets, selected_filters, FACET_FIELDS, "targets", id_prefix="targets"
    )

    pagination_style = {"display": "flex"} if total_pages > 1 else {"display": "none"}
    page_info = f"Page {current_page_number} of {total_pages} ({total} results)"
    prev_disabled = current_page_number <= 1
    next_disabled = current_page_number >= total_pages

    return (
        table,
        pagination_style,
        page_info,
        prev_disabled,
        next_disabled,
        current_page_number,
        filters_children,
        {"display": "flex", "justifyContent": "flex-end"},
    )


@callback(
    Output("targets-download", "data"),
    Input("targets-download-btn", "n_clicks"),
    State("targets-results-store", "data"),
    State({"type": "targets-facet-filter", "field": ALL}, "value"),
    State({"type": "targets-facet-filter", "field": ALL}, "id"),
    prevent_initial_call=True,
)
def download_targets_metadata(n_clicks, store_data, selected_values, filter_ids):
    """Download all target results as CSV."""
    if not n_clicks or not store_data:
        return dash.no_update

    query = store_data.get("query")

    selected_filters = {}
    if selected_values and filter_ids:
        for values, filter_id in zip(selected_values, filter_ids):
            if values:
                cleaned = [str(v).strip() for v in values if v not in (None, "")]
                if cleaned:
                    selected_filters[filter_id["field"]] = cleaned

    all_results = fetch_all_search_results(
        query=query if query else None,
        filters=selected_filters or None,
        search_mode="targets",
    )

    if not all_results:
        return dash.no_update

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
            # Sanitize formula injection characters
            if value and value[0:1] in ("=", "+", "-", "@", "\t", "\r"):
                value = "'" + value
            if "," in value or '"' in value or "\n" in value:
                value = '"' + value.replace('"', '""') + '"'
            row_values.append(value)
        csv_lines.append(",".join(row_values))

    csv_content = "\n".join(csv_lines)
    safe_query = "".join(
        c if c.isalnum() or c in "-_" else "_" for c in (query or "all")[:30]
    )
    filename = f"perturbation_catalogue_targets_{safe_query}.csv"

    return dict(content=csv_content, filename=filename, type="text/csv")
