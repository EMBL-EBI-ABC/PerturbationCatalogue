import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import requests


data_portal_layout = html.Div(
    [
        dbc.Row(
            [
                # Left column for filter cards
                dbc.Col(
                    [
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Sequence Type", className="card-title"),
                                    dbc.RadioItems(
                                        id="sequenceType",
                                        options=[],
                                        inline=True,
                                        style={"width": "100%"},
                                    ),
                                    dbc.Button(
                                        "Clear",
                                        id="clear-sequenceType",
                                        color="success",
                                        size="sm",
                                        style={"margin-top": "10px"},
                                    ),
                                ]
                            )
                        ),
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Gene Category", className="card-title"),
                                    dbc.RadioItems(
                                        id="geneCategory",
                                        options=[],
                                        inline=True,
                                        style={"width": "100%"},
                                    ),
                                    dbc.Button(
                                        "Clear",
                                        id="clear-geneCategory",
                                        color="success",
                                        size="sm",
                                        style={"margin-top": "10px"},
                                    ),
                                ]
                            )
                        ),
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Publication Year", className="card-title"),
                                    dbc.RadioItems(
                                        id="publicationYear",
                                        options=[],
                                        inline=True,
                                        style={"width": "100%"},
                                    ),
                                    dbc.Button(
                                        "Clear",
                                        id="clear-publicationYear",
                                        color="success",
                                        size="sm",
                                        style={"margin-top": "10px"},
                                    ),
                                ]
                            )
                        ),
                    ],
                    width=2,
                    style={
                        "padding": "30px 0px 0px 20px",
                        "display": "flex",
                        "flexDirection": "column",
                        "gap": "12px",
                    },
                ),
                # Right column for table and pagination controls
                dbc.Col(
                    [
                        dcc.Input(
                            id="search",
                            type="text",
                            placeholder="Search...",
                            debounce=True,
                            style={
                                "width": "100%",
                                "margin-bottom": "15px",
                            },
                        ),
                        html.Div(
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
                                    style={"width": "auto"},
                                ),
                                html.Span(id="pagination-info"),
                            ],
                            style={
                                "display": "flex",
                                "alignItems": "center",
                                "gap": "12px",
                            },
                        ),
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
                        html.Div(
                            id="data-table",
                            style={
                                "height": "100%",
                                "overflowY": "auto",
                                "width": "100%",
                            },
                        ),
                    ],
                    width=10,
                    style={"padding": "32px 25px 100px 25px"},
                ),
            ],
            className="g-0",
        ),
    ]
)


def details_layout(urn):
    api_url = f"https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search/{urn}"
    response = requests.get(api_url)

    if response.status_code != 200:
        return html.Div(
            "Error: Unable to fetch data from the API.", className="alert alert-danger"
        )

    data = response.json().get("results", [{}])[0]

    return html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H4(
                                        data.get("title", "N/A"), className="card-title"
                                    ),
                                    html.P(
                                        data.get("shortDescription", "N/A"),
                                        className="card-text",
                                    ),
                                    html.Div(
                                        [
                                            html.A(
                                                "View on MaveDB",
                                                href=f"https://www.mavedb.org/score-sets/{urn}/",
                                                className="btn btn-primary mb-3",
                                            )
                                        ]
                                    ),
                                    html.P(
                                        [html.Strong("URN: "), data.get("urn", "N/A")],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Sequence Type: "),
                                            data.get("sequenceType", "N/A"),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Gene Name: "),
                                            data.get("geneName", "N/A"),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Gene Category: "),
                                            data.get("geneCategory", "N/A"),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Publication: "),
                                            html.A(
                                                data.get("publicationUrl", "N/A"),
                                                href=data.get("publicationUrl", "#"),
                                                target="_blank",
                                            ),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Publication Year: "),
                                            str(data.get("publicationYear", "N/A")),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Number of Variants: "),
                                            str(data.get("numVariants", "N/A")),
                                        ],
                                        className="card-text",
                                    ),
                                ],
                                className="card-body",
                            )
                        ],
                        className="card",
                    )
                ],
                className="col-md-8 mx-auto",
            )
        ],
        className="container mt-4",
    )


def data_portal_callbacks(app):

    @app.callback(
        [
            Output("data-table", "children"),
            Output("sequenceType", "options"),
            Output("geneCategory", "options"),
            Output("publicationYear", "options"),
            Output("pagination", "max_value"),
            Output("pagination-info", "children"),
        ],
        [
            Input("search", "value"),
            Input("size", "value"),
            Input("pagination", "active_page"),
            Input("sequenceType", "value"),
            Input("geneCategory", "value"),
            Input("publicationYear", "value"),
        ],
    )
    def fetch_data(q, size, page, sequenceType, geneCategory, publicationYear):
        base_url = "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search"
        start = (page - 1) * size if page else 0

        params = {
            k: v
            for k, v in {
                "q": q,
                "size": size,
                "start": start,
                "sequenceType": sequenceType,
                "geneCategory": geneCategory,
                "publicationYear": publicationYear,
            }.items()
            if v not in [None, ""]
        }

        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            data = response.json()
            total = data.get("total", 0)
            max_pages = max(1, (total + size - 1) // size) if total > 0 else 1
            start_index = start + 1 if total > 0 else 0
            end_index = min(start + size, total)

            sequenceType_options = [
                {"label": b["key"], "value": b["key"]}
                for b in data.get("aggregations", {})
                .get("sequenceType", {})
                .get("buckets", [])
            ]
            geneCategory_options = [
                {"label": b["key"], "value": b["key"]}
                for b in data.get("aggregations", {})
                .get("geneCategory", {})
                .get("buckets", [])
            ]
            publicationYear_options = [
                {"label": str(b["key"]), "value": b["key"]}
                for b in data.get("aggregations", {})
                .get("publicationYear", {})
                .get("buckets", [])
            ]

            pagination_info = f"{start_index} - {end_index} of {total}"

            if data["results"]:
                columns = [
                    "URN",
                    "Sequence Type",
                    "Gene Name",
                    "Gene Category",
                    "Publication Year",
                    "Number of Variants",
                ]
                table_header = [html.Thead(html.Tr([html.Th(col) for col in columns]))]
                rows = []
                for row in data["results"]:
                    urn_link = html.A(
                        row["urn"],
                        href=f"/data-portal/{row['urn']}",
                        style={"whiteSpace": "nowrap", "textDecoration": "none"},
                    )
                    rows.append(
                        html.Tr(
                            [
                                html.Td(urn_link),
                                html.Td(row.get("sequenceType", "N/A")),
                                html.Td(row.get("geneName", "N/A")),
                                html.Td(row.get("geneCategory", "N/A")),
                                html.Td(row.get("publicationYear", "N/A")),
                                html.Td(row.get("numVariants", "N/A")),
                            ]
                        )
                    )
                table_body = [html.Tbody(rows)]
                table = dbc.Table(
                    table_header + table_body,
                    bordered=True,
                    hover=True,
                    responsive=True,
                )
            else:
                table = html.Div("No data found.")

            return (
                table,
                sequenceType_options,
                geneCategory_options,
                publicationYear_options,
                max_pages,
                pagination_info,
            )
        else:
            return html.Div("Error fetching data."), [], [], [], 1, ""

    @app.callback(
        [
            Output("sequenceType", "value"),
            Output("geneCategory", "value"),
            Output("publicationYear", "value"),
        ],
        [
            Input("clear-sequenceType", "n_clicks"),
            Input("clear-geneCategory", "n_clicks"),
            Input("clear-publicationYear", "n_clicks"),
        ],
    )
    def clear_filters(clear_seq, clear_gene, clear_pub):
        triggered_id = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

        if triggered_id == "clear-sequenceType":
            return [None, dash.no_update, dash.no_update]
        elif triggered_id == "clear-geneCategory":
            return [dash.no_update, None, dash.no_update]
        elif triggered_id == "clear-publicationYear":
            return [dash.no_update, dash.no_update, None]

        return [dash.no_update, dash.no_update, dash.no_update]


if __name__ == "__main__":
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    app.layout = data_portal_layout
    data_portal_callbacks(app)
    app.run_server(debug=True)
