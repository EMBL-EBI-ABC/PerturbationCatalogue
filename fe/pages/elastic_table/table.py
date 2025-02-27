from dash import dcc, html
import dash_bootstrap_components as dbc


def layout(filter_fields):
    return html.Div(
        [
            dcc.Interval(
                id="elastic-table-timer", interval=1000, max_intervals=1, disabled=True
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
                            for field in filter_fields
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
                                data={"field": "publicationYear", "order": "desc"},
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
