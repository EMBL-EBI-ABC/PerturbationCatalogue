import dash
from dash import html, dcc, callback, Input, Output
import dash_bootstrap_components as dbc
import requests
import os

# BE endpoint
api_base_url = os.getenv("PERTURBATION_CATALOGUE_BE", "http://127.0.0.1:8000")

# Page registration
dash.register_page(
    __name__,
    path="/perturbations",
    relative_path="/perturbations",
    name="Perturbations",
    button="Open Perturbations",
    description="A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.",
    icon="bi-whirlpool",
)


def format_value(value):
    if isinstance(value, float):
        return f"{value:.2e}".replace("-", "−")
    if isinstance(value, str):
        return value.capitalize()
    return value


def layout():
    return dbc.Container(
        fluid=True,
        children=[
            # Title
            html.H1(
                "Perturbation Catalogue",
                className="text-center display-4 mt-5 mb-3",
                style={"font-size": "60px"},
            ),
            # Subtitle
            html.P(
                "A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.",
                className="text-center lead mb-5",
                style={"font-size": "25px"},
            ),
            # Main content
            html.Div(
                className="px-5",
                children=[
                    # Headers
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.H4(
                                        "According to this",
                                        className="fw-normal fst-italic mb-0",
                                    ),
                                    html.H3("Dataset"),
                                ],
                                width=4,
                            ),
                            dbc.Col(
                                [
                                    html.H4(
                                        "Introducing this",
                                        className="fw-normal fst-italic mb-0",
                                    ),
                                    html.H3("Perturbation"),
                                ],
                                width=2,
                            ),
                            dbc.Col(
                                [
                                    html.H4(
                                        "Leads to the following",
                                        className="fw-normal fst-italic mb-0",
                                    ),
                                    html.H3("Effect"),
                                ],
                                width=6,
                            ),
                        ]
                    ),
                    # Search fields
                    dbc.Row(
                        [
                            dbc.Col(
                                dbc.Input(
                                    id="dataset-filter",
                                    placeholder="Filter by dataset metadata",
                                ),
                                width=4,
                            ),
                            dbc.Col(
                                dbc.Input(
                                    id="perturbation-filter",
                                    placeholder="Filter by perturbed gene",
                                ),
                                width=2,
                            ),
                            dbc.Col(
                                dbc.Input(
                                    id="effect-filter",
                                    placeholder="Filter by affected gene",
                                ),
                                width=6,
                            ),
                        ],
                        className="mb-4",
                    ),
                    # Data grid
                    html.Div(id="data-grid"),
                ],
            ),
        ],
    )


def register_callbacks(app):
    @callback(
        Output("data-grid", "children"),
        [
            Input("dataset-filter", "value"),
            Input("perturbation-filter", "value"),
            Input("effect-filter", "value"),
        ],
    )
    def update_data_grid(dataset_metadata, perturbation_gene_name, effect_gene_name):
        params = {
            "dataset_metadata": dataset_metadata,
            "perturbation_gene_name": perturbation_gene_name,
            "effect_gene_name": effect_gene_name,
        }
        # Filter out None values
        params = {k: v for k, v in params.items() if v}

        try:
            response = requests.get(f"{api_base_url}/v1/search", params=params)
            response.raise_for_status()
            data = response.json()
        except requests.exceptions.RequestException as e:
            return dbc.Alert(f"Error fetching data from backend: {e}", color="danger")

        if not data:
            return dbc.Alert("No results found.", color="info")

        grid_layout = []
        for dataset in data:
            if not dataset.get("perturbations"):
                continue

            dataset_info = dbc.Col(
                [
                    html.B(dataset["dataset_id"]),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.I(
                                        f'{key.replace("_label", "").replace("_", " ").capitalize()}\t'
                                    ),
                                    html.Span(
                                        format_value(value), className="fw-semibold"
                                    ),
                                ]
                            )
                            for key, value in dataset.items()
                            if key != "dataset_id" and key != "perturbations" and value
                        ]
                    ),
                ],
                width=4,
            )

            perturbation_rows = []
            for i, p in enumerate(dataset["perturbations"]):
                log2fc = p.get("log2foldchange")
                arrow = ""
                if log2fc is not None:
                    if log2fc > 0:
                        arrow = "▲"
                    elif log2fc < 0:
                        arrow = "▼"

                perturbation_rows.append(
                    dbc.Row(
                        [
                            dbc.Col(
                                html.B(p["perturbation_gene_name"]),
                                width=3,
                            ),
                            dbc.Col(
                                html.Div(
                                    [
                                        html.B(f'{arrow} {p["effect_gene_name"]}'),
                                        html.Span(" ", className="me-3"),
                                        html.I("log2fc "),
                                        html.Span(
                                            format_value(p.get("log2foldchange")),
                                            className="fw-semibold me-3",
                                        ),
                                        html.I("padj "),
                                        html.Span(
                                            format_value(p.get("padj")),
                                            className="fw-semibold me-3",
                                        ),
                                        html.I("base mean "),
                                        html.Span(
                                            format_value(p.get("base_mean")),
                                            className="fw-semibold",
                                        ),
                                    ],
                                    className="d-flex align-items-start flex-wrap",
                                ),
                                width=9,
                            ),
                        ],
                        align="start",
                        className="border-top" if i > 0 else "",
                    )
                )

            grid_layout.append(
                dbc.Row(
                    [dataset_info, dbc.Col(perturbation_rows, width=8)],
                    align="start",
                    className="border-top py-3",
                )
            )

        return grid_layout
