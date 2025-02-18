import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable
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
                        "padding": "30px 0px 0px 20px",  # Increased outer padding for more space
                        "display": "flex",
                        "flexDirection": "column",
                        "gap": "12px",  # Reduced gap between cards
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
                            },  # Reduced margin bottom
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
                                "gap": "12px",  # Slightly reduced gap between controls
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
                                "margin-top": "15px",  # Consistent margin for pagination
                                "color": "white",
                            },
                        ),
                        html.Div(
                            DataTable(
                                id="data-table",
                                columns=[],  # Columns will be updated dynamically
                                data=[],  # Data will be updated dynamically
                                style_table={
                                    "height": "100%",
                                    "overflowY": "auto",
                                    "width": "100%",
                                },
                                style_cell={
                                    "whiteSpace": "normal",  # Wrap text inside cells
                                    "textAlign": "left",  # Left-align text inside cells
                                    "padding": "5px",  # Add padding for better readability
                                    # "fontSize": "13px",  # Smaller font size for the table
                                    "fontFamily": "system-ui",  # Set font to system UI
                                },
                                style_data_conditional=[
                                    {
                                        "if": {"column_id": col},
                                        "textAlign": "left",
                                    }
                                    for col in []  # Apply to all columns
                                ],
                            )
                        ),
                    ],
                    width=10,  # 80% width for table and pagination controls
                    style={
                        "padding": "32px 25px 100px 25px"
                    },  # Increased outer padding for more space
                ),
            ],
            className="g-0",  # Remove gutters
        ),
    ]
)


def details_layout(urn):
    # Make the API call to fetch data
    api_url = f"https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search/{urn}"
    response = requests.get(api_url)

    # Check if the API call was successful
    if response.status_code != 200:
        return html.Div(
            "Error: Unable to fetch data from the API.", className="alert alert-danger"
        )

    # Extract data from the API response
    data = response.json().get("results", [{}])[
        0
    ]  # Default to an empty dict if no results

    # Create the Bootstrap card layout
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

    # API call function
    @app.callback(
        [
            Output("data-table", "columns"),
            Output("data-table", "data"),
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

            # Dynamically generate columns based on keys in the first result
            if data["results"]:
                # Exclude the columns you want to hide
                columns = [
                    {"name": key, "id": key}
                    for key in data["results"][0].keys()
                    if key not in ["shortDescription", "publicationUrl", "title"]
                ]
                table_data = data["results"]
            else:
                columns = []
                table_data = []

            return (
                columns,
                table_data,
                sequenceType_options,
                geneCategory_options,
                publicationYear_options,
                max_pages,
                pagination_info,
            )
        else:
            return [], [], [], [], [], 1, ""

    # Clearing the filters

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
        # Reset values when respective "Clear" buttons are clicked
        triggered_id = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

        if triggered_id == "clear-sequenceType":
            return [None, dash.no_update, dash.no_update]
        elif triggered_id == "clear-geneCategory":
            return [dash.no_update, None, dash.no_update]
        elif triggered_id == "clear-publicationYear":
            return [dash.no_update, dash.no_update, None]

        # Return current values if no button was clicked
        return [dash.no_update, dash.no_update, dash.no_update]


# Run the app
if __name__ == "__main__":
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    app.layout = data_portal_layout
    data_portal_callbacks(app)
    app.run_server(debug=True)
