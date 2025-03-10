"""A reusable and configurable class for table and details views populated by the Elastic API."""

from collections import namedtuple

import dash
from dash import html, dcc, Output, Input, State
import dash_bootstrap_components as dbc
import requests

Column = namedtuple(
    "ColumnDefinition",
    [
        "field_name",
        "display_name",
        "display_table",
        "display_details",
        "filterable",
        "sortable",
        "default_sort",
    ],
    defaults=[None, None, lambda x: x, None, False, False, False],
)


class ElasticTable:
    def __init__(
        self,
        api_endpoint,
        columns,
        details_button_name,
        details_button_link,
        title=None,
    ):
        self.api_endpoint = api_endpoint
        self.columns = columns
        self.details_button_name = details_button_name
        self.details_button_link = details_button_link
        self.title = title

    # Table view.

    def _get_table_columns(self):
        """Returns columns that should be displayed in the table view."""
        return [col for col in self.columns if col.display_table]

    def _fetch_data(self, q, size, page, sort_data, filter_values):
        """Fetches data from the API endpoint with the given parameters."""
        params = {
            "q": q,
            "size": size,
            "start": (page - 1) * size if page else 0,
            **{
                field_name: value
                for field_name, value in zip(
                    [col.field_name for col in self.columns if col.filterable],
                    filter_values,
                )
                if value
            },
            "sort_field": sort_data.get("field"),
            "sort_order": sort_data.get("order"),
        }

        response = requests.get(
            self.api_endpoint,
            params={k: v for k, v in params.items() if v is not None},
        )

        return response, params

    def _create_table_header(self, col, current_sort):
        """Creates a table header cell with optional sort indicator."""
        if not col.sortable:
            return html.Th(col.display_name)

        is_sorted = current_sort.get("field") == col.field_name
        order = current_sort.get("order") if is_sorted else None
        icon = {"asc": "sort-up", "desc": "sort-down"}.get(order, "")

        content = [col.display_name]
        if icon:
            content.append(html.I(className=f"bi bi-{icon}"))

        return html.Th(
            html.Div(
                content,
                className="d-flex align-items-center gap-1",
                style={
                    "cursor": "pointer",
                    "color": "var(--custom-color)" if is_sorted else "inherit",
                },
            ),
            id={"type": "sort-header-container", "field": col.field_name},
        )

    def _create_table(self, data, sort_data):
        """Creates a Dash Bootstrap table from the provided data."""
        if not data:
            return html.Div("No data found.")

        return dbc.Table(
            [
                html.Thead(
                    html.Tr(
                        [
                            self._create_table_header(col, sort_data)
                            for col in self._get_table_columns()
                        ]
                    )
                ),
                html.Tbody(
                    [
                        html.Tr(
                            [
                                html.Td(
                                    col.display_table(row.get(col.field_name, "N/A"))
                                )
                                for col in self._get_table_columns()
                            ]
                        )
                        for row in data
                    ]
                ),
            ],
            bordered=True,
            hover=True,
            responsive=True,
        )

    def table_layout(self):
        """Returns the complete layout for the table view."""
        default_sort_column = next(
            (col for col in self.columns if col.default_sort), None
        )

        # Prepare filter column components
        filter_columns = [
            dbc.Card(
                dbc.CardBody(
                    [
                        html.H5(col.display_name),
                        dbc.Checklist(
                            id=col.field_name,
                            options=[],
                            className="w-100 elastic-table-filter-checklist",
                        ),
                        dbc.Button(
                            "Clear",
                            id=f"clear-{col.field_name}",
                            color="success",
                            size="sm",
                            className="mt-2 elastic-table-filter-clear",
                        ),
                    ]
                )
            )
            for col in self._get_table_columns()
            if col.filterable
        ]

        # Table title block
        title_block = []
        if self.title:
            title_block = [html.H2(self.title, className="ps-4 pt-4 mb-0")]

        return html.Div(
            [
                *title_block,
                dbc.Row(
                    [
                        # Filters sidebar
                        dbc.Col(
                            filter_columns,
                            xs=12,
                            md=3,
                            className="pt-4 ps-4 d-flex flex-column gap-3",
                        ),
                        # Main content area
                        dbc.Col(
                            [
                                # Search field
                                dcc.Input(
                                    id="search",
                                    type="text",
                                    placeholder="Search...",
                                    className="mb-3 w-100",
                                ),
                                # Store to track debounce state
                                dcc.Store(id="search-input-value", data=""),
                                # Current sort direction store
                                dcc.Store(
                                    id="sort-store",
                                    data={
                                        "field": default_sort_column.field_name,
                                        "order": default_sort_column.default_sort,
                                    },
                                ),
                                # Main table with spinner
                                dbc.Spinner(
                                    html.Div(
                                        id="data-table",
                                        className="w-100 overflow-auto",
                                        style={"minHeight": "100px"},
                                    ),
                                    color="primary",
                                    spinner_style={"width": "48px", "height": "48px"},
                                ),
                                # Pagination controls
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            [
                                                html.Label("Items per page:"),
                                                dcc.Dropdown(
                                                    id="size",
                                                    options=[
                                                        {"label": str(i), "value": i}
                                                        for i in [20, 50, 100, 200]
                                                    ],
                                                    value=20,
                                                    clearable=False,
                                                    style={"width": "70px"},
                                                ),
                                                html.Span(id="pagination-info"),
                                            ],
                                            width="auto",
                                            className="d-flex align-items-center gap-2 me-3",
                                        ),
                                        dbc.Col(
                                            dbc.Pagination(
                                                id="pagination",
                                                max_value=1,
                                                active_page=1,
                                                first_last=True,
                                                previous_next=True,
                                                fully_expanded=False,
                                                className="mt-3",
                                            ),
                                            width="auto",
                                        ),
                                    ],
                                    className="justify-content-end g-0 mt-0",
                                ),
                            ],
                            className="pt-4 pe-4 ps-3",
                        ),
                    ],
                    className="g-0",
                ),
            ],
            className="w-100",
        )

    # Details view.

    def _get_detail(self, record_id):
        """Fetches details for a specific item by its URN."""
        response = requests.get(f"{self.api_endpoint}/{record_id}")
        return response.json()

    def _format_title(self, data):
        """Formats the title section of the details view."""
        title_field = next(
            col.field_name for col in self.columns if col.display_details == "title"
        )
        return html.H4(
            data.get(title_field, "N/A"),
            className="card-title",
        )

    def _format_subtitle(self, data):
        """Formats the subtitle section of the details view."""
        subtitle_field = next(
            col.field_name for col in self.columns if col.display_details == "subtitle"
        )
        return html.P(
            data.get(subtitle_field, "N/A"),
            className="card-text",
        )

    def _format_button(self, urn):
        """Creates a button for the details view."""
        return html.A(
            self.details_button_name,
            href=self.details_button_link(urn),
            className="btn btn-primary mb-3",
        )

    def _format_item(self, data, col):
        """Formats a single item in the details view."""
        if col.display_details == "text":
            return html.P(
                [
                    html.Strong(f"{col.display_name}: "),
                    str(data.get(col.field_name, "N/A")),
                ],
                className="card-text",
            )
        else:
            return html.P(
                [
                    html.Strong(f"{col.display_name}: "),
                    html.A(
                        data.get(col.field_name, "N/A"),
                        href=data.get(col.field_name, "#"),
                        target="_blank",
                    ),
                ],
                className="card-text",
            )

    def details_layout(self, record_id):
        """Returns the complete layout for the details view."""
        data = self._get_detail(record_id).get("results", [{}])

        if not data:
            return dbc.Alert(
                "Error: Unable to fetch data from the API.",
                color="danger",
            )

        data = data[0]

        return html.Div(
            dbc.Row(
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                self._format_title(data),
                                self._format_subtitle(data),
                                self._format_button(record_id),
                                *[
                                    self._format_item(data, col)
                                    for col in self.columns
                                    if col.display_details in ("text", "link")
                                ],
                            ]
                        )
                    ),
                    width={"size": 8},
                    className="mx-auto",
                ),
                className="mt-4 g-0",
            ),
            className="w-100",
        )

    # Callbacks.

    def register_callbacks(self, app):
        """Registers all necessary Dash callbacks for the table functionality."""

        # Clientside callback for text input debounce.
        app.clientside_callback(
            """
        function(value, oldValue) {
            if (value === oldValue) {
                return window.dash_clientside.no_update;
            }

            // If search is empty or being cleared, update immediately without debounce
            if (!value || value === '') {
                return value;
            }

            // Clear any existing timeouts
            if (window.searchDebounceTimeout) {
                clearTimeout(window.searchDebounceTimeout);
            }

            return new Promise((resolve) => {
                window.searchDebounceTimeout = setTimeout(() => {
                    resolve(value);
                }, 750); // Debounce in ms
            });
        }
        """,
            Output("search-input-value", "data"),
            Input("search", "value"),
            State("search-input-value", "data"),
        )

        @app.callback(
            Output("sort-store", "data"),
            Input({"type": "sort-header-container", "field": dash.ALL}, "n_clicks"),
            State("sort-store", "data"),
        )
        def update_sort(_, current_sort):
            """Updates the sort direction when a column header is clicked."""
            if not dash.callback_context.triggered:
                return dash.no_update

            field = eval(dash.callback_context.triggered[0]["prop_id"].split(".")[0])[
                "field"
            ]
            if current_sort["field"] != field:
                return {"field": field, "order": "asc"}
            return {
                "field": field,
                "order": "desc" if current_sort["order"] == "asc" else "asc",
            }

        @app.callback(
            [
                Output("data-table", "children"),
                *[
                    Output(col.field_name, "options")
                    for col in self._get_table_columns()
                    if col.filterable
                ],
                Output("pagination", "max_value"),
                Output("pagination-info", "children"),
            ],
            [
                Input("search-input-value", "data"),
                Input("size", "value"),
                Input("pagination", "active_page"),
                Input("sort-store", "data"),
                *[
                    Input(col.field_name, "value")
                    for col in self._get_table_columns()
                    if col.filterable
                ],
            ],
        )
        def fetch_data(q, size, page, sort_data, *filter_values):
            """Fetches and refreshes the table data based on current filters and parameters."""
            response, params = self._fetch_data(q, size, page, sort_data, filter_values)

            if response.status_code != 200:
                return (
                    html.Div("Error fetching data."),
                    *[[]]
                    * len([col for col in self._get_table_columns() if col.filterable]),
                    1,
                    "",
                )

            data = response.json()
            total = data.get("total", 0)
            start = params["start"]

            filter_options = [
                [
                    {"label": f"{b['key']} ({b['doc_count']})", "value": b["key"]}
                    for b in data.get("aggregations", {})
                    .get(col.field_name, {})
                    .get("buckets", [])
                ]
                for col in self._get_table_columns()
                if col.filterable
            ]

            return (
                self._create_table(data.get("results", []), sort_data),
                *filter_options,
                max(1, (total + size - 1) // size) if total > 0 else 1,
                (
                    f"{start + 1} – {min(start + size, total)} of {total}"
                    if total > 0
                    else ""
                ),
            )

        @app.callback(
            [
                Output(col.field_name, "value")
                for col in self._get_table_columns()
                if col.filterable
            ],
            [
                Input(f"clear-{col.field_name}", "n_clicks")
                for col in self._get_table_columns()
                if col.filterable
            ],
        )
        def clear_filters(*_):
            """Clears the selected filter values when a clear button is clicked."""
            triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
            return [
                [] if f"clear-{col.field_name}" == triggered else dash.no_update
                for col in self._get_table_columns()
                if col.filterable
            ]
