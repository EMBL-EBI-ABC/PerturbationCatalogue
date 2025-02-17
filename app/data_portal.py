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
        dcc.Input(id="sequenceType", type="text", placeholder="sequenceType"),
        dcc.Input(id="geneCategory", type="text", placeholder="geneCategory"),
        dcc.Input(id="publicationYear", type="text", placeholder="publicationYear"),
        html.Pre(id="output", style={"white-space": "pre-wrap"}),
    ]
)


# API call function
@app.callback(
    Output("output", "children"),
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
        "q": q or "",
        "size": size,
        "start": start,
        "sequenceType": sequenceType or "",
        "geneCategory": geneCategory or "",
        "publicationYear": publicationYear or "",
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return json.dumps(response.json(), indent=2)
    else:
        return f"Error: {response.status_code}\n{response.text}"


# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)
