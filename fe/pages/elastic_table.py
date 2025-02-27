from collections import namedtuple

import dash
from dash import html, dcc
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
    ):
        self.api_endpoint = api_endpoint
        self.columns = columns
        self.details_button_name = details_button_name
        self.details_button_link = details_button_link

    # Table view.

    def _get_table_columns(self):
        return [col for col in self.columns if col.display_table]

    def _fetch_data(self, q, size, page, sort_data, filter_values):
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
        if not col.sortable:
            return html.Th(col.display_name)

        is_sorted = current_sort.get("field") == col.field_name
        order = current_sort.get("order") if is_sorted else None
        icon = {"asc": "sort-up", "desc": "sort-down"}.get(order, "")
        content = (
            [col.display_name, html.I(className=f"bi bi-{icon}")]
            if icon
            else col.display_name
        )

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
        default_sort_column = next(
            (col for col in self.columns if col.default_sort), None
        )
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
                                        ],
                                    )
                                )
                                for col in self._get_table_columns()
                                if col.filterable
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
                                        "field": default_sort_column.field_name,
                                        "order": default_sort_column.default_sort,
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

    # Details view.

    def _get_detail(self, urn):
        response = requests.get(f"{self.api_endpoint}/{urn}")
        return response.json()

    def _format_title(self, data):
        title_field = [
            col.field_name for col in self.columns if col.display_details == "title"
        ][0]
        return html.H4(
            data.get(title_field, "N/A"),
            className="card-title",
        )

    def _format_subtitle(self, data):
        subtitle_field = [
            col.field_name for col in self.columns if col.display_details == "subtitle"
        ][0]
        return html.P(
            data.get(subtitle_field, "N/A"),
            className="card-text",
        )

    def _format_button(self, urn):
        return html.Div(
            [
                html.A(
                    self.details_button_name,
                    href=self.details_button_link(urn),
                    className="btn btn-primary mb-3",
                )
            ]
        )

    def _format_item(self, data, col):
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

    def details_layout(self, urn):
        data = self._get_detail(urn).get("results", [{}])

        if not data:
            return html.Div(
                "Error: Unable to fetch data from the API.",
                className="alert alert-danger",
            )

        data = data[0]

        return dbc.Container(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Card(
                                [
                                    dbc.CardBody(
                                        [
                                            self._format_title(data),
                                            self._format_subtitle(data),
                                            self._format_button(urn),
                                            *[
                                                self._format_item(data, col)
                                                for col in self.columns
                                                if col.display_details
                                                in ("text", "link")
                                            ],
                                        ]
                                    )
                                ]
                            ),
                            width={"size": 8, "offset": 2},
                        )
                    ]
                )
            ],
            className="mt-4",
        )

    # Callbacks.

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
                *[
                    dash.Output(col.field_name, "options")
                    for col in self._get_table_columns()
                    if col.filterable
                ],
                dash.Output("pagination", "max_value"),
                dash.Output("pagination-info", "children"),
            ],
            [
                dash.State("search", "value"),
                dash.Input("elastic-table-timer", "n_intervals"),
                dash.Input("size", "value"),
                dash.Input("pagination", "active_page"),
                dash.Input("sort-store", "data"),
                *[
                    dash.Input(col.field_name, "value")
                    for col in self._get_table_columns()
                    if col.filterable
                ],
            ],
            dash.State("search", "value"),
        )
        def fetch_data(q, n_intervals, size, page, sort_data, *filter_values):
            if n_intervals is not None and n_intervals != 1:
                return dash.no_update

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
                dash.Output(col.field_name, "value")
                for col in self._get_table_columns()
                if col.filterable
            ],
            [
                dash.Input(f"clear-{col.field_name}", "n_clicks")
                for col in self._get_table_columns()
                if col.filterable
            ],
        )
        def clear_filters(*_):
            triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
            return [
                [] if f"clear-{col.field_name}" == triggered else dash.no_update
                for col in self._get_table_columns()
                if col.filterable
            ]
