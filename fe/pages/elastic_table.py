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
        id,
        api_endpoint,
        columns,
        details_button_name,
        details_button_link,
        title=None,
        description=None,
        default_page_size=20,
    ):
        # A globally unique DOM prefix, based on the ID, to distinguish this table from all other ElasticTable instances.
        self.dom_prefix = f"elastic-table-{id}"
        self.api_endpoint = api_endpoint
        self.columns = columns
        self.details_button_name = details_button_name
        self.details_button_link = details_button_link
        self.title = title
        self.description = description
        self.default_page_size = default_page_size

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
            id={
                "type": f"{self.dom_prefix}-sort-header-container",
                "field": col.field_name,
            },
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
                            id=f"{self.dom_prefix}-filter-{col.field_name}",
                            options=[],
                            className="w-100 elastic-table-filter-checklist",
                        ),
                        dbc.Button(
                            "Clear",
                            id=f"{self.dom_prefix}-clear-{col.field_name}",
                            color="success",
                            size="sm",
                            className="mt-2 elastic-table-filter-clear",
                        ),
                    ]
                )
            )
            for col in self.columns
            if col.filterable
        ]

        # Table title block
        title_block = []
        if self.title:
            title_block = [html.H2(self.title, className="ps-4 pt-4 mb-0")]

        # Table description block
        description_block = []
        if self.description:
            description_block = [html.H5(self.description, className="ps-4 pt-3 mb-0")]

        return html.Div(
            [
                *title_block,
                *description_block,
                dbc.Row(
                    [
                        # Filters sidebar
                        dbc.Col(
                            filter_columns,
                            xs=12,
                            md=2,
                            className="pt-4 ps-4 d-flex flex-column gap-3",
                        ),
                        # Main content area
                        dbc.Col(
                            [
                                # Search field
                                dcc.Input(
                                    id=f"{self.dom_prefix}-search",
                                    type="text",
                                    placeholder="Search...",
                                    className="mb-3 w-100",
                                ),
                                # Store to track debounce state
                                dcc.Store(
                                    id=f"{self.dom_prefix}-search-input-value", data=""
                                ),
                                # Current sort direction store
                                dcc.Store(
                                    id=f"{self.dom_prefix}-sort-store",
                                    data={
                                        "field": default_sort_column.field_name,
                                        "order": default_sort_column.default_sort,
                                    },
                                ),
                                # Main table with spinner
                                dbc.Spinner(
                                    html.Div(
                                        id=f"{self.dom_prefix}-data-table",
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
                                                    id=f"{self.dom_prefix}-size",
                                                    options=[
                                                        {"label": str(i), "value": i}
                                                        for i in [10, 20, 50, 100, 200]
                                                    ],
                                                    value=self.default_page_size,
                                                    clearable=False,
                                                    style={"width": "70px"},
                                                ),
                                                html.Span(
                                                    id=f"{self.dom_prefix}-pagination-info"
                                                ),
                                            ],
                                            width="auto",
                                            className="d-flex align-items-center gap-2 me-3",
                                        ),
                                        dbc.Col(
                                            dbc.Pagination(
                                                id=f"{self.dom_prefix}-pagination",
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
            className="w-100 mb-4",
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
        elif col.display_details == "link":
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
        else:
            return col.display_details(data.get(col.field_name, "N/A"))

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
                                    if col.display_details
                                    and col.display_details not in ("title", "subtitle")
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
            function(value, state) {
                // Special case for initial load - check if this is the first callback execution
                if (!window.initialLoadComplete<id>) {
                    window.initialLoadComplete<id> = true;
                    return value;  // Allow the value through on first load
                }

                if (value === state.search) {
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
            """.replace(
                "<id>", self.id
            ),
            Output(f"{self.dom_prefix}-search", "value"),
            Input(f"{self.dom_prefix}-search", "value"),
            State(f"{self.dom_prefix}-state", "data"),
        )

        # Callback 1: Update state when any input changes
        @app.callback(
            [
                Output(f"{self.dom_prefix}-state", "data"),
                *[
                    Output(f"{self.dom_prefix}-filter-{col.field_name}", "value")
                    for col in self.columns
                    if col.filterable
                ],
                Output(f"{self.dom_prefix}-pagination", "active_page"),
            ],
            [
                Input(f"{self.dom_prefix}-search", "value"),
                Input(f"{self.dom_prefix}-size", "value"),
                Input(f"{self.dom_prefix}-pagination", "active_page"),
                Input(
                    {
                        "type": f"{self.dom_prefix}-sort-header-container",
                        "field": dash.ALL,
                    },
                    "n_clicks",
                ),
                *[
                    Input(f"{self.dom_prefix}-filter-{col.field_name}", "value")
                    for col in self.columns
                    if col.filterable
                ],
                *[
                    Input(f"{self.dom_prefix}-clear-{col.field_name}", "n_clicks")
                    for col in self.columns
                    if col.filterable
                ],
            ],
            [State(f"{self.dom_prefix}-state", "data")],
        )
        def update_state(search, size, page, sort_clicks, *args):
            """Updates the state store based on user interactions."""
            ctx = dash.callback_context
            if not ctx.triggered:
                return dash.no_update

            # Split args into filter values and clear button clicks
            filterable_cols = [col for col in self.columns if col.filterable]
            num_filterable = len(filterable_cols)
            filter_values = args[:num_filterable]
            clear_clicks = args[num_filterable : num_filterable * 2]
            current_state = args[-1]  # Last argument is the current state

            triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

            # By default, no pagination reset is required.
            pagination_page_update = dash.no_update

            # Initialize updated filter values to current values
            updated_filter_values = list(filter_values)

            # Handle sort header clicks
            if "sort-header-container" in triggered_id:
                field = eval(triggered_id)["field"]
                current_sort = current_state["sort"]

                if current_sort["field"] != field:
                    current_state["sort"] = {"field": field, "order": "asc"}
                else:
                    current_state["sort"] = {
                        "field": field,
                        "order": "desc" if current_sort["order"] == "asc" else "asc",
                    }
                pagination_page_update = 1

            # Handle clear filter buttons
            elif "clear" in triggered_id:
                for i, col in enumerate(filterable_cols):
                    if f"{self.dom_prefix}-clear-{col.field_name}" == triggered_id:
                        current_state["filters"][i] = []
                        updated_filter_values[i] = []  # Update the filter value
                        pagination_page_update = 1
                        break

            # Handle normal inputs
            else:
                # Update search
                if f"{self.dom_prefix}-search" == triggered_id:
                    current_state["search"] = search
                    # Reset to page 1 when search changes, except when this is initial load.
                    if not current_state["initial_load"]:
                        pagination_page_update = 1
                    else:
                        current_state["initial_load"] = False

                # Update page size
                elif f"{self.dom_prefix}-size" == triggered_id:
                    current_state["size"] = size
                    pagination_page_update = 1

                # Update page number
                elif f"{self.dom_prefix}-pagination" == triggered_id:
                    current_state["page"] = page

                # Update filters
                else:
                    for i, col in enumerate(filterable_cols):
                        if f"{self.dom_prefix}-filter-{col.field_name}" == triggered_id:
                            current_state["filters"][i] = filter_values[i]
                            pagination_page_update = 1
                            break

            # If we need to update the page, update it also in the state.
            if pagination_page_update == 1:
                current_state["page"] = pagination_page_update

            return [current_state, *updated_filter_values, pagination_page_update]

        # Callback 2: Fetch data when state changes
        @app.callback(
            [
                Output(f"{self.dom_prefix}-data-table", "children"),
                *[
                    Output(f"{self.dom_prefix}-filter-{col.field_name}", "options")
                    for col in self.columns
                    if col.filterable
                ],
                Output(f"{self.dom_prefix}-pagination", "max_value"),
                Output(f"{self.dom_prefix}-pagination-info", "children"),
            ],
            [Input(f"{self.dom_prefix}-state", "data")],
        )
        def fetch_data_from_state(state):
            """Fetches and refreshes the table data based on current state."""
            response, params = self._fetch_data(state)

            if response.status_code != 200:
                return (
                    html.Div("Error fetching data."),
                    *[[]] * len([col for col in self.columns if col.filterable]),
                    1,
                    "",
                )

            data = response.json()
            total = data.get("total", 0)
            size = state.get("size", self.default_page_size)
            start = params["start"]

            filter_options = [
                [
                    {"label": f"{b['key']} ({b['doc_count']})", "value": b["key"]}
                    for b in data.get("aggregations", {})
                    .get(col.field_name, {})
                    .get("buckets", [])
                ]
                for col in self.columns
                if col.filterable
            ]

            max_pages = math.ceil(total / size)

            return (
                self._create_table(data.get("results", []), state.get("sort", {})),
                *filter_options,
                max_pages,
                (
                    f"{start + 1} – {min(start + size, total)} of {total}"
                    if total > 0
                    else ""
                ),
            )
