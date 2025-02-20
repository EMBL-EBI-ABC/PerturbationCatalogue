from dash import dcc, html
import dash_bootstrap_components as dbc


def layout(filter_fields):
    return html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Card(
                                dbc.CardBody(
                                    [
                                        html.H5(field.title, className="card-title"),
                                        dbc.RadioItems(
                                            id=field.id,
                                            options=[],
                                            inline=True,
                                            style={"width": "100%"},
                                        ),
                                        dbc.Button(
                                            "Clear",
                                            id=f"clear-{field.id}",
                                            color="success",
                                            size="sm",
                                            style={"margin-top": "10px"},
                                        ),
                                    ]
                                )
                            )
                            for field in filter_fields
                        ],
                        width=2,
                        lg=2,
                        md=4,
                        sm=12,
                        xs=12,
                        style={
                            "padding": "30px 0px 0px 20px",
                            "display": "flex",
                            "flexDirection": "column",
                            "gap": "12px",
                        },
                    ),
                    dbc.Col(
                        [
                            dcc.Input(
                                id="search",
                                type="text",
                                placeholder="Search...",
                                debounce=True,
                                style={
                                    "width": "100%",
                                    "margin-bottom": "20px",
                                },
                            ),
                            dcc.Store(
                                id="sort-store",
                                data={"field": "publicationYear", "order": "desc"},
                            ),
                            dbc.Spinner(
                                html.Div(
                                    id="data-table",
                                    style={
                                        "overflowY": "auto",
                                        "width": "100%",
                                        "minHeight": "100px",
                                    },
                                ),
                                color="primary",
                                spinner_style={"width": "48px", "height": "48px"},
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        dbc.Pagination(
                                            id="pagination",
                                            max_value=1,
                                            active_page=1,
                                            first_last=True,
                                            fully_expanded=False,
                                            previous_next=True,
                                            style={
                                                "margin-top": "15px",
                                                "color": "white",
                                            },
                                        ),
                                        width="auto",
                                        className="me-auto",
                                    ),
                                    dbc.Col(
                                        [
                                            html.Label("Items per page:"),
                                            dcc.Dropdown(
                                                id="size",
                                                options=[
                                                    {"label": str(i), "value": i}
                                                    for i in [10, 50, 100, 200]
                                                ],
                                                value=10,
                                                clearable=False,
                                                style={"width": "70px"},
                                            ),
                                            html.Span(id="pagination-info"),
                                        ],
                                        width="auto",
                                        className="d-flex align-items-center gap-2",
                                    ),
                                ],
                                className="g-0 mt-0",
                            ),
                        ],
                        width=10,
                        lg=10,
                        md=8,
                        sm=12,
                        xs=12,
                        style={"padding": "32px 25px 25px 25px"},
                    ),
                ],
                className="g-0",
            ),
        ]
    )
