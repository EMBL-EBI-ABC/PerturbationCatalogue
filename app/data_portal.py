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
        dcc.Input(id="search", type="text", placeholder="Search...", debounce=True),
        dcc.Dropdown(id="sequenceType", placeholder="Select sequenceType"),
        dcc.Dropdown(id="geneCategory", placeholder="Select geneCategory"),
        dcc.Dropdown(id="publicationYear", placeholder="Select publicationYear"),
        html.Div(
            [
                html.Label("Items per page:"),
                dcc.Dropdown(
                    id="size",
                    options=[{"label": str(i), "value": i} for i in [10, 50, 100, 200]],
                    value=10,
                    clearable=False,
                ),
                html.Span(id="pagination-info"),
                dbc.Pagination(
                    id="pagination",
                    max_value=1,
                    active_page=1,
                    first_last=True,
                    fully_expanded=False,
                    previous_next=True,
                ),
            ],
            style={"display": "flex", "alignItems": "center", "gap": "10px"},
        ),
        DataTable(
            id="data-table",
            style_table={
                "width": "100%",
                "maxWidth": "100%",
                "overflowX": "auto",  # Allow for table overflow if needed, no horizontal scroll
            },
            style_cell={
                "whiteSpace": "normal",  # Ensures text wraps in cells
                "textOverflow": "ellipsis",  # Ensures ellipsis for overflow
                "maxWidth": "200px",  # Set a max width for each cell (optional)
            },
            style_header={"textAlign": "center"},
            style_data={"whiteSpace": "normal"},
            style_cell_conditional=[
                {
                    "if": {
                        "column_id": "column_name"
                    },  # Replace with actual column names if needed
                    "textAlign": "center",
                }
            ],
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

        # Dynamically generate columns based on keys in the first result (assuming all items have the same keys)
        if data["results"]:
            columns = [{"name": key, "id": key} for key in data["results"][0].keys()]
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
