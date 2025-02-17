import dash
from dash import dcc, html, Input, Output
import requests
import json

# Initialize the Dash app
app = dash.Dash(__name__)

# Define the layout
app.layout = html.Div(
    [
        dcc.Input(id="search", type="text", placeholder="Search...", debounce=True),
        dcc.Input(id="size", type="number", placeholder="Page size", value=10, min=1),
        dcc.Input(id="start", type="number", placeholder="Start", value=0, min=0),
        dcc.Dropdown(id="sequenceType", placeholder="Select sequenceType"),
        dcc.Dropdown(id="geneCategory", placeholder="Select geneCategory"),
        dcc.Dropdown(id="publicationYear", placeholder="Select publicationYear"),
        html.Pre(id="output", style={"white-space": "pre-wrap"}),
    ]
)


# API call function
@app.callback(
    [
        Output("output", "children"),
        Output("sequenceType", "options"),
        Output("geneCategory", "options"),
        Output("publicationYear", "options"),
    ],
    [
        Input("search", "value"),
        Input("size", "value"),
        Input("start", "value"),
        Input("sequenceType", "value"),
        Input("geneCategory", "value"),
        Input("publicationYear", "value"),
    ],
)
def fetch_data(q, size, start, sequenceType, geneCategory, publicationYear):
    base_url = "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search"
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
        aggregations = data.get("aggregations", {})
        sequenceType_options = [
            {"label": b["key"], "value": b["key"]}
            for b in aggregations.get("sequenceType", {}).get("buckets", [])
        ]
        geneCategory_options = [
            {"label": b["key"], "value": b["key"]}
            for b in aggregations.get("geneCategory", {}).get("buckets", [])
        ]
        publicationYear_options = [
            {"label": str(b["key"]), "value": b["key"]}
            for b in aggregations.get("publicationYear", {}).get("buckets", [])
        ]
        return (
            json.dumps(data, indent=2),
            sequenceType_options,
            geneCategory_options,
            publicationYear_options,
        )
    else:
        return f"Error: {response.status_code}\n{response.text}", [], [], []


# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)
