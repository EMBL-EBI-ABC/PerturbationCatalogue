import dash
from dash import html, callback, Output, Input, dcc
from time import sleep
import dash_bootstrap_components as dbc

dash.register_page(
    __name__,
    path="/dwh_interface",
    name="DWH Interface",
    button="DWH Interface",
    description="Data warehousing advance query interfaces and presentations",
    icon="bi-database-fill-gear",
)


def build_table() -> dbc.Table:
    table_header = [html.Thead(html.Tr([html.Th("GeName"), html.Th("GeneProb")]))]
    table_body = [
        html.Tbody([html.Tr([html.Td("YRDC"), html.Td("0.99930764779822245")])])
    ]
    return dbc.Table(
        table_header + table_body,
        striped=True,
        bordered=True,
        hover=True,
        responsive=True,
    )


tab1_content = dbc.Card(
    dbc.CardBody(
        [
            html.P("Write SQL query in a box below", className="card-text"),
            dbc.Textarea(
                id="sql_query",
                size="lg",
                className="mb-3",
                rows=15,
                debounce=True,
                style={
                    "backgroundColor": "#000",
                    "border": "1px solid #000",
                    "color": "#00ff00",
                    "padding": "8px",
                    "fontFamily": "courier new",
                },
            ),
            dbc.Spinner(html.Span(id="validation_text")),
            dbc.Button(
                "Run Query",
                color="success",
                style={"marginTop": "15px", "marginBottom": "15px"},
                id="run-sql-query-button",
            ),
            dbc.Spinner(html.Div(id="sql-query-table")),
        ]
    ),
    className="mt-3",
)

tab2_content = dbc.Card(
    dbc.CardBody(
        [
            html.H3("Choose fields", className="card-text"),
            dbc.Container(
                dbc.Row(
                    [
                        dbc.Col(
                            dcc.Dropdown(
                                ["GeneName", "GeneProb"],
                                "GeneName",
                            ),
                            md=4,
                            id="field-1",
                        ),
                        dbc.Col(
                            dcc.Dropdown(
                                ["GeneName", "GeneProb"],
                                "GeneName",
                            ),
                            md=4,
                            id="field-2",
                        ),
                    ]
                )
            ),
            html.P(
                html.Span(
                    [
                        "Accept fields? ",
                        dbc.Button(
                            "Yes",
                            color="success",
                            style={"marginRight": "5px"},
                            id="fields-yes-button",
                        ),
                        dbc.Button("No", color="danger"),
                    ]
                ),
                style={"marginTop": "15px", "marginBottom": "15px"},
            ),
            html.H3("Choose filters"),
            dbc.Container(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(
                                    [
                                        "Primary Disease",
                                        "Age Category",
                                        "Sex",
                                        "GeneProb",
                                        "Essential Gene",
                                    ],
                                    "Primary Disease",
                                ),
                                md=4,
                                id="filter-1-a",
                            ),
                            dbc.Col(
                                dcc.Dropdown(["=", "!=", ">=", "<="], "="),
                                md=4,
                                id="filter-1-b",
                            ),
                            dbc.Col(
                                dbc.Input(placeholder="Type something..."),
                                md=4,
                                id="filter-1-c",
                            ),
                        ],
                        style={"marginBottom": "5px"},
                    ),
                    dbc.Row(
                        dbc.Col(
                            dcc.Dropdown(["AND", "OR", "NOT"], "OR"), md=4, id="conj-1"
                        ),
                        style={"marginBottom": "5px"},
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(
                                    [
                                        "Primary Disease",
                                        "Age Category",
                                        "Sex",
                                        "GeneProb",
                                        "Essential Gene",
                                    ],
                                    "Primary Disease",
                                ),
                                md=4,
                                id="filter-2-a",
                            ),
                            dbc.Col(
                                dcc.Dropdown(["=", "!=", ">=", "<="], "="),
                                md=4,
                                id="filter-2-b",
                            ),
                            dbc.Col(
                                dbc.Input(placeholder="Type something..."),
                                md=4,
                                id="filter-2-c",
                            ),
                        ],
                        style={"marginBottom": "5px"},
                    ),
                    dbc.Row(
                        dbc.Col(
                            dcc.Dropdown(["AND", "OR", "NOT"], "OR"), md=4, id="conj-2"
                        ),
                        style={"marginBottom": "5px"},
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(
                                    [
                                        "Primary Disease",
                                        "Age Category",
                                        "Sex",
                                        "GeneProb",
                                        "Essential Gene",
                                    ],
                                    "Primary Disease",
                                ),
                                md=4,
                                id="filter-3-a",
                            ),
                            dbc.Col(
                                dcc.Dropdown(["=", "!=", ">=", "<="], "="),
                                md=4,
                                id="filter-3-b",
                            ),
                            dbc.Col(
                                dbc.Input(placeholder="Type something..."),
                                md=4,
                                id="filter-3-c",
                            ),
                        ],
                        style={"marginBottom": "5px"},
                    ),
                    dbc.Row(
                        dbc.Col(
                            dcc.Dropdown(["AND", "OR", "NOT"], "OR"), md=4, id="conj-3"
                        ),
                        style={"marginBottom": "5px"},
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(
                                    [
                                        "Primary Disease",
                                        "Age Category",
                                        "Sex",
                                        "GeneProb",
                                        "Essential Gene",
                                    ],
                                    "Primary Disease",
                                ),
                                md=4,
                                id="filter-4-a",
                            ),
                            dbc.Col(
                                dcc.Dropdown(["=", "!=", ">=", "<="], "="),
                                md=4,
                                id="filter-4-b",
                            ),
                            dbc.Col(
                                dbc.Input(placeholder="Type something..."),
                                md=4,
                                id="filter-4-c",
                            ),
                        ],
                        style={"marginBottom": "5px"},
                    ),
                    dbc.Row(
                        dbc.Col(
                            dcc.Dropdown(["AND", "OR", "NOT"], "OR"), md=4, id="conj-4"
                        ),
                        style={"marginBottom": "5px"},
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(
                                    [
                                        "Primary Disease",
                                        "Age Category",
                                        "Sex",
                                        "GeneProb",
                                        "Essential Gene",
                                    ],
                                    "Primary Disease",
                                ),
                                md=4,
                                id="filter-5-a",
                            ),
                            dbc.Col(
                                dcc.Dropdown(["=", "!=", ">=", "<="], "="),
                                md=4,
                                id="filter-5-b",
                            ),
                            dbc.Col(
                                dbc.Input(placeholder="Type something..."),
                                md=4,
                                id="filter-5-c",
                            ),
                        ],
                        style={"marginBottom": "5px"},
                    ),
                ]
            ),
            html.P(
                html.Span(
                    [
                        "Accept Filters? ",
                        dbc.Button(
                            "Yes",
                            color="success",
                            style={"marginRight": "5px"},
                            id="filters-yes-button",
                        ),
                        dbc.Button("No", color="danger"),
                    ]
                ),
                style={"marginTop": "15px", "marginBottom": "15px"},
            ),
            html.P(
                html.Span(
                    [
                        "Run Query? ",
                        dbc.Button("Yes", color="info", id="query-builder-button"),
                    ],
                )
            ),
            dbc.Spinner(html.Div(id="builder-query-table")),
        ]
    ),
    className="mt-3",
)

tab3_content = dbc.Card(
    dbc.CardBody(
        [
            html.P("Write your question in a box below", className="card-text"),
            dbc.Textarea(
                size="lg",
                className="mb-3",
                rows=15,
                debounce=True,
                style={
                    "backgroundColor": "#000",
                    "border": "1px solid #000",
                    "color": "#00ff00",
                    "padding": "8px",
                    "fontFamily": "courier new",
                },
            ),
            dbc.Button(
                "Run Query",
                color="success",
                style={"marginTop": "15px", "marginBottom": "15px"},
                id="run-nl-query-button",
            ),
            dbc.Spinner(html.Div(id="nl-query-table")),
        ]
    ),
    className="mt-3",
)

layout = dbc.Container(
    dbc.Row(
        dbc.Col(
            dbc.Tabs(
                [
                    dbc.Tab(tab1_content, label="SQL Interface"),
                    dbc.Tab(tab2_content, label="Advance Query Builder"),
                    dbc.Tab(
                        tab3_content,
                        label="Natural Language Query Interface",
                    ),
                ]
            ),
            md=12,
        ),
        style={"marginTop": "2em", "marginBottom": "2em"},
    )
)


@callback(
    Output("sql_query", "valid"),
    Output("sql_query", "invalid"),
    Output("validation_text", "children"),
    Input("sql_query", "value"),
)
def validate_query(query):
    if query is not None and "SELECT" in query:
        sleep(3)
        return (
            True,
            False,
            dbc.Badge(
                "Validation Results: This query will process 1.12 GB when run.",
                color="info",
            ),
        )
    elif query is not None and "SELECT" not in query:
        sleep(3)
        return (
            False,
            True,
            dbc.Badge("Validation Results: Invalid query!", color="danger"),
        )
    else:
        return False, False, ""


@callback(
    Output("run-sql-query-button", "disabled"),
    Output("sql-query-table", "children"),
    Input("run-sql-query-button", "n_clicks"),
)
def run_sql_query(n_clicks):
    if n_clicks is not None:
        sleep(5)
        return True, build_table()
    else:
        return False, None


@callback(
    Output("field-1", "children"),
    Output("field-2", "children"),
    Input("fields-yes-button", "n_clicks"),
    prevent_initial_call=True,
)
def build_fields(n_clicks):
    if n_clicks is not None:
        return dbc.Badge("GeneName", color="info", pill=True), dbc.Badge(
            "GeneProb", color="info", pill=True
        )


@callback(
    Output("filter-1-a", "children"),
    Output("filter-1-b", "children"),
    Output("filter-1-c", "children"),
    Output("conj-1", "children"),
    Output("filter-2-a", "children"),
    Output("filter-2-b", "children"),
    Output("filter-2-c", "children"),
    Output("conj-2", "children"),
    Output("filter-3-a", "children"),
    Output("filter-3-b", "children"),
    Output("filter-3-c", "children"),
    Output("conj-3", "children"),
    Output("filter-4-a", "children"),
    Output("filter-4-b", "children"),
    Output("filter-4-c", "children"),
    Output("conj-4", "children"),
    Output("filter-5-a", "children"),
    Output("filter-5-b", "children"),
    Output("filter-5-c", "children"),
    Input("filters-yes-button", "n_clicks"),
    prevent_initial_call=True,
)
def build_filters(n_clicks):
    if n_clicks is not None:
        return (
            dbc.Badge("Primary Disease", color="primary", pill=True),
            dbc.Badge("=", color="secondary", pill=True),
            dbc.Badge("Synovial Sarcoma", color="primary", pill=True),
            dbc.Badge("AND", color="success", pill=True),
            dbc.Badge("Age Category", color="primary", pill=True),
            dbc.Badge("=", color="secondary", pill=True),
            dbc.Badge("Pediatric", color="primary", pill=True),
            dbc.Badge("AND", color="success", pill=True),
            dbc.Badge("Sex", color="primary", pill=True),
            dbc.Badge("=", color="secondary", pill=True),
            dbc.Badge("Male", color="primary", pill=True),
            dbc.Badge("AND", color="success", pill=True),
            dbc.Badge("GeneProb", color="primary", pill=True),
            dbc.Badge(">=", color="secondary", pill=True),
            dbc.Badge("0.99", color="primary", pill=True),
            dbc.Badge("AND", color="success", pill=True),
            dbc.Badge("Essential Gene", color="primary", pill=True),
            dbc.Badge("=", color="secondary", pill=True),
            dbc.Badge("False", color="primary", pill=True),
        )


@callback(
    Output("builder-query-table", "children"),
    Input("query-builder-button", "n_clicks"),
)
def table_builder(n_clicks):
    sleep(2)
    if n_clicks is not None:
        return build_table()


@callback(
    Output("run-nl-query-button", "disabled"),
    Output("nl-query-table", "children"),
    Input("run-nl-query-button", "n_clicks"),
)
def run_nl_query(n_clicks):
    if n_clicks is not None:
        sleep(10)
        return True, build_table()
    else:
        return False, None
