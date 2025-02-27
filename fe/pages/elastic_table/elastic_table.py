import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
import requests


class ElasticTable:
    def __init__(
        self,
        api_endpoint,
        columns,
        filter_fields,
        sortable_columns,
        default_sort_field,
        default_sort_order="desc",
    ):
        self.api_endpoint = api_endpoint
        self.columns = columns
        self.filter_fields = filter_fields
        self.sortable_columns = sortable_columns
        self.default_sort_field = default_sort_field
        self.default_sort_order = default_sort_order

    def table_layout(self):
        return html.Div(
            [
                dcc.Interval(
                    id="elastic-table-timer",
                    interval=1000,
                    max_intervals=1,
                    disabled=True,
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dbc.Card(
                                    dbc.CardBody(
                                        [
                                            html.H5(field.title),
                                            dbc.Checklist(
                                                id=field.id,
                                                options=[],
                                                className="w-100 elastic-table-filter-checklist",
                                            ),
                                            dbc.Button(
                                                "Clear",
                                                id=f"clear-{field.id}",
                                                color="success",
                                                size="sm",
                                                className="mt-2 elastic-table-filter-clear",
                                            ),
                                        ],
                                    )
                                )
                                for field in self.filter_fields
                            ],
                            width=2,
                            lg=2,
                            md=4,
                            sm=12,
                            xs=12,
                            className="pt-4 ps-4 d-flex flex-column gap-3",
                        ),
                        dbc.Col(
                            [
                                dcc.Input(
                                    id="search",
                                    type="text",
                                    placeholder="Search...",
                                    className="mb-3 w-100",
                                ),
                                dcc.Store(
                                    id="sort-store",
                                    data={
                                        "field": self.default_sort_field,
                                        "order": self.default_sort_order,
                                    },
                                ),
                                dbc.Spinner(
                                    html.Div(
                                        id="data-table",
                                        className="w-100 overflow-auto",
                                        style={"minHeight": "100px"},
                                    ),
                                    color="primary",
                                    spinner_style={"width": "48px", "height": "48px"},
                                ),
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
                                    className="g-0 mt-0 justify-content-end",
                                ),
                            ],
                            width=10,
                            lg=10,
                            md=8,
                            sm=12,
                            xs=12,
                            className="pt-4 pe-4 ps-3",
                        ),
                    ],
                    className="g-0",
                ),
            ],
            style={"width": "100%"},
        )

    def create_table_header(self, column_name, field_name, current_sort):
        if field_name not in self.sortable_columns and field_name is not None:
            return html.Th(column_name)

        is_sorted = current_sort.get("field") == field_name
        order = current_sort.get("order") if is_sorted else None
        icon = "sort-up" if order == "asc" else "sort-down" if order == "desc" else ""

        return (
            html.Th(
                html.Div(
                    (
                        [column_name, html.I(className=f"bi bi-{icon}")]
                        if icon
                        else column_name
                    ),
                    className="d-flex align-items-center gap-1",
                    style={
                        "cursor": "pointer",
                        "color": "var(--custom-color)" if is_sorted else "inherit",
                    },
                ),
                id={"type": "sort-header-container", "field": field_name},
            )
            if field_name is not None
            else html.Th(column_name)
        )

    def create_table(self, data, sort_data):
        if not data:
            return html.Div("No data found.")

        rows = [
            html.Tr(
                [
                    (
                        html.Td(
                            html.A(
                                row["urn"],
                                href=f"{self.get_detail_url(row['urn'])}",
                                className="text-decoration-none text-nowrap",
                            )
                        )
                        if column[0] == "URN"
                        else html.Td(row.get(column[1], "N/A") if column[1] else "N/A")
                    )
                    for column in self.columns
                ]
            )
            for row in data
        ]

        return dbc.Table(
            [
                html.Thead(
                    html.Tr(
                        [
                            self.create_table_header(col[0], col[1], sort_data)
                            for col in self.columns
                        ]
                    )
                ),
                html.Tbody(rows),
            ],
            bordered=True,
            hover=True,
            responsive=True,
        )

    def get_detail_url(self, urn):
        return f"/data-portal/{urn}"

    def fetch_data(self, q, size, page, sort_data, filter_values):
        params = {
            "q": q,
            "size": size,
            "start": (page - 1) * size if page else 0,
            **{f.id: v for f, v in zip(self.filter_fields, filter_values) if v},
            "sort_field": sort_data.get("field"),
            "sort_order": sort_data.get("order"),
        }

        response = requests.get(
            self.api_endpoint,
            params={k: v for k, v in params.items() if v is not None},
        )

        return response, params

    def register_callbacks(self, app):
        @app.callback(
            dash.Output("sort-store", "data"),
            dash.Input(
                {"type": "sort-header-container", "field": dash.ALL}, "n_clicks"
            ),
            dash.State("sort-store", "data"),
        )
        def update_sort(_, current_sort):
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
            dash.Output("elastic-table-timer", "n_intervals"),
            dash.Output("elastic-table-timer", "disabled"),
            dash.Input("search", "value"),
            prevent_initial_call=True,
        )
        def start_timer(value):
            return 0, False

        @app.callback(
            [
                dash.Output("data-table", "children"),
                *[dash.Output(f.id, "options") for f in self.filter_fields],
                dash.Output("pagination", "max_value"),
                dash.Output("pagination-info", "children"),
            ],
            [
                dash.State("search", "value"),
                dash.Input("elastic-table-timer", "n_intervals"),
                dash.Input("size", "value"),
                dash.Input("pagination", "active_page"),
                dash.Input("sort-store", "data"),
                *[dash.Input(f.id, "value") for f in self.filter_fields],
            ],
            dash.State("search", "value"),
        )
        def fetch_data(q, n_intervals, size, page, sort_data, *filter_values):
            if n_intervals is not None and n_intervals != 1:
                return dash.no_update

            response, params = self.fetch_data(q, size, page, sort_data, filter_values)

            if response.status_code != 200:
                return (
                    html.Div("Error fetching data."),
                    *[[]] * len(self.filter_fields),
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
                    .get(f.id, {})
                    .get("buckets", [])
                ]
                for f in self.filter_fields
            ]

            return (
                self.create_table(data.get("results", []), sort_data),
                *filter_options,
                max(1, (total + size - 1) // size) if total > 0 else 1,
                (
                    f"{start + 1} – {min(start + size, total)} of {total}"
                    if total > 0
                    else ""
                ),
            )

        @app.callback(
            [dash.Output(f.id, "value") for f in self.filter_fields],
            [dash.Input(f"clear-{f.id}", "n_clicks") for f in self.filter_fields],
        )
        def clear_filters(*_):
            triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
            return [
                [] if f"clear-{f.id}" == triggered else dash.no_update
                for f in self.filter_fields
            ]

    def get_detail(self, urn):
        response = requests.get(f"{self.api_endpoint}/{urn}")

        if response.status_code != 200:
            return html.Div(
                "Error: Unable to fetch data from the API.",
                className="alert alert-danger",
            )

        return response.json()
