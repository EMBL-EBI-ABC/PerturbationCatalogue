import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable
import requests
import json

# Initialize the Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Define the layout
app.layout = html.Div(
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
                                ]
                            )
                        ),
                    ],
                    width=3,  # 20% width for filters
                    style={
                        "padding": "10px",
                        "display": "flex",
                        "flexDirection": "column",
                        "gap": "10px",
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
                            style={"width": "100%", "margin-bottom": "10px"},
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
                                "gap": "10px",
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
                                "margin-top": "10px",
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
                                    "fontSize": "12px",  # Smaller font size for the table
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
                    width=9,  # 80% width for table and pagination controls
                    style={"padding": "10px"},
                ),
            ],
            className="g-0",  # Remove gutters
        ),
    ]
)


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


# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)
