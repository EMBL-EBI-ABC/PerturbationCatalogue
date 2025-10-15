import os
import dash
from dash import html, dcc, callback, Input, Output, State
import dash_bootstrap_components as dbc
import requests
import json

dash.register_page(
    __name__,
    path="/perturbations",
    relative_path="/perturbations",
    name="Perturbations",
    button="Explore perturbations",
    description="Search and filter CRISPR, MAVE and Perturb-Seq perturbation results.",
    icon="bi-tropical-storm",
)


# Helper function to format values
def format_value(value):
    if isinstance(value, float):
        return f"{value:.2e}".replace("-", "−")
    if isinstance(value, str):
        return value.capitalize() if value.lower() not in ["grch38"] else value
    return value


def get_arrow(log2fc):
    if log2fc > 0:
        return html.Span("▲", style={"color": "green"})
    elif log2fc < 0:
        return html.Span("▼", style={"color": "red"})
    return None


def create_label_value_pair(label, value, value_is_bold=True, spacing="me-3"):
    return html.Div(
        [
            html.I(f"{label} ", className="text-muted"),
            html.Span(value, className="fw-semibold" if value_is_bold else ""),
        ],
        className=f"d-inline-block {spacing}",
    )


def render_effect_row(effect, is_grouped=False):
    elements = [
        get_arrow(effect.get("log2foldchange", 0)),
        html.B(f" {effect['effect_gene_name']}"),
        html.Span(" ", className="me-3"),
    ]

    if "log2foldchange" in effect:
        elements.append(
            create_label_value_pair("log2fc", format_value(effect["log2foldchange"]))
        )
    if "padj" in effect:
        elements.append(create_label_value_pair("padj", format_value(effect["padj"])))
    if "base_mean" in effect:
        elements.append(
            create_label_value_pair("base mean", format_value(effect["base_mean"]))
        )

    return dbc.Row(
        dbc.Col(html.Div(elements, className="d-flex align-items-start flex-wrap"))
    )


def render_grouped_by_perturbation(dataset):
    rows = []
    for pert in dataset["data"]:
        pert_info = dbc.Col(
            [
                html.H5(pert["perturbation_gene_name"], className="fw-bold"),
                html.Div(["Affects ", html.B(pert["n_total"]), " genes"]),
                html.Div(
                    [html.Span("▲", style={"color": "green"}), f" {pert['n_up']} genes"]
                ),
                html.Div(
                    [html.Span("▼", style={"color": "red"}), f" {pert['n_down']} genes"]
                ),
            ],
            width=3,
        )

        effects_rows = [
            render_effect_row(eff, is_grouped=True) for eff in pert.get("effects", [])
        ]
        effects_col = dbc.Col(effects_rows, width=9)

        rows.append(dbc.Row([pert_info, effects_col], className="py-2"))

    return dbc.Col(rows, width=8)


def render_grouped_by_effect(dataset):
    rows = []
    for eff in dataset["data"]:
        pert_rows = []
        for pert in eff.get("perturbations", []):
            pert_col = dbc.Col(html.B(pert["perturbation_gene_name"]), width=3)

            effect_details = [
                get_arrow(pert.get("log2foldchange", 0)),
                html.Span(" ", className="me-3"),
            ]
            if "log2foldchange" in pert:
                effect_details.append(
                    create_label_value_pair(
                        "log2fc", format_value(pert["log2foldchange"])
                    )
                )
            if "padj" in pert:
                effect_details.append(
                    create_label_value_pair("padj", format_value(pert["padj"]))
                )

            effect_col = dbc.Col(
                html.Div(
                    effect_details, className="d-flex align-items-start flex-wrap"
                ),
                width=9,
            )
            pert_rows.append(
                dbc.Row([pert_col, effect_col], className="border-top py-2")
            )

        effect_info = dbc.Col(
            [
                html.H5(eff["effect_gene_name"], className="fw-bold"),
                html.Div(["Affected by ", html.B(eff["n_total"]), " genes"]),
                html.Div(
                    [
                        html.Span("▲", style={"color": "green"}),
                        " by ",
                        html.B(eff["n_up"]),
                        " genes",
                    ]
                ),
                html.Div(
                    [
                        html.Span("▼", style={"color": "red"}),
                        " by ",
                        html.B(eff["n_down"]),
                        " genes",
                    ]
                ),
                html.Div(
                    ["Base mean expression is ", html.B(format_value(eff["base_mean"]))]
                ),
            ],
            width=9,
        )

        # This is tricky. The effect info is in the right column, and the left is empty.
        # The perturbations are listed below.
        # The spec is a bit ambiguous here. Let's try to make it reasonable.
        # The effect info should be in the effect column. The perturbation column for the "header" is empty.

        # Let's restructure. The effect info is the header of the right column.
        # The left column is empty for the header.
        # Then the perturbation rows follow.

        # The spec says "display the perturbation rows". This implies they are separate.
        # "perturbed gene should be displayed in the perturbation column, and log2fc, padj must be displayed in the Effect column"

        # Let's create a container for the whole effect block
        effect_block = dbc.Row(
            [dbc.Col(width=3), effect_info]  # Empty perturbation column for the header
        )

        rows.append(effect_block)
        rows.extend(pert_rows)

    # This structure is different. The grouping info is on the right.
    # So the main split is 4 (dataset) and 8 (pert+effect).
    # Inside the 8, we have the grouped data.
    # The spec says "display the effect gene name in a larger font" in the *effect* column.
    # This means the perturbation column for that "header" row is empty.

    # Let's try to build the final structure for one dataset
    final_rows = []
    for eff in dataset["data"]:
        effect_info_col = dbc.Col(
            [
                html.H5(eff["effect_gene_name"], className="fw-bold"),
                html.Div(["Affected by ", html.B(eff["n_total"]), " genes"]),
                html.Div(
                    [
                        html.Span("▲", style={"color": "green"}),
                        " by ",
                        html.B(eff["n_up"]),
                        " genes",
                    ]
                ),
                html.Div(
                    [
                        html.Span("▼", style={"color": "red"}),
                        " by ",
                        html.B(eff["n_down"]),
                        " genes",
                    ]
                ),
                html.Div(
                    ["Base mean expression is ", html.B(format_value(eff["base_mean"]))]
                ),
            ],
            width=9,
        )

        header_row = dbc.Row([dbc.Col(width=3), effect_info_col], className="py-2")
        final_rows.append(header_row)

        for pert in eff.get("perturbations", []):
            pert_col = dbc.Col(html.B(pert["perturbation_gene_name"]), width=3)

            effect_details = [
                get_arrow(pert.get("log2foldchange", 0)),
                html.Span(" ", className="me-3"),
            ]
            if "log2foldchange" in pert:
                effect_details.append(
                    create_label_value_pair(
                        "log2fc", format_value(pert["log2foldchange"])
                    )
                )
            if "padj" in pert:
                effect_details.append(
                    create_label_value_pair("padj", format_value(pert["padj"]))
                )

            effect_col = dbc.Col(
                html.Div(
                    effect_details, className="d-flex align-items-start flex-wrap"
                ),
                width=9,
            )
            final_rows.append(dbc.Row([pert_col, effect_col], className="py-2"))

    return dbc.Col(final_rows, width=8)


def render_ungrouped(dataset):
    perturbation_rows = []
    for item in dataset["data"]:
        pert_col = dbc.Col(html.B(item["perturbation_gene_name"]), width=3)

        effect_elements = [
            get_arrow(item["log2foldchange"]),
            html.B(f" {item['effect_gene_name']}"),
            html.Span(" ", className="me-3"),
            create_label_value_pair("log2fc", format_value(item["log2foldchange"])),
            create_label_value_pair("padj", format_value(item["padj"])),
            create_label_value_pair("base mean", format_value(item["base_mean"])),
        ]
        effect_col = dbc.Col(
            html.Div(effect_elements, className="d-flex align-items-start flex-wrap"),
            width=9,
        )

        perturbation_rows.append(
            dbc.Row([pert_col, effect_col], align="start", className="py-2")
        )

    return dbc.Col(perturbation_rows, width=8)


def render_dataset_row(dataset, group_by):
    dataset_info = [html.B(dataset["dataset_id"])]
    for key, value in dataset.items():
        if key.endswith("_label") and value:
            label = key.replace("_label", "").replace("_", " ").capitalize()
            dataset_info.append(
                html.Div(
                    [
                        html.I(f"{label} ", className="text-muted"),
                        html.Span(format_value(value), className="fw-semibold"),
                    ]
                )
            )

    dataset_col = dbc.Col(dataset_info, width=4)

    if group_by == "perturbation_gene_name":
        data_col = render_grouped_by_perturbation(dataset)
    elif group_by == "effect_gene_name":
        data_col = render_grouped_by_effect(dataset)
    else:
        data_col = render_ungrouped(dataset)

    return dbc.Row([dataset_col, data_col], align="start", className="border-top py-3")


# Layout definition
def layout():
    return dbc.Container(
        fluid=True,
        className="px-5",
        children=[
            dcc.Store(id="group-by-store", data="perturbation_gene_name"),
            html.H1(
                "Perturbation Catalogue",
                className="text-center display-4 mt-5 mb-3",
                style={"font-size": "60px"},
            ),
            html.P(
                "A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.",
                className="text-center lead mb-5",
                style={"font-size": "25px"},
            ),
            # Headers
            dbc.Row(
                [
                    dbc.Col(
                        html.Div(
                            [
                                html.H4(
                                    "According to this",
                                    className="fw-normal fst-italic mb-0",
                                ),
                                html.H3("Dataset"),
                            ]
                        ),
                        width=4,
                    ),
                    dbc.Col(
                        html.Div(
                            [
                                html.H4(
                                    "Introducing this",
                                    className="fw-normal fst-italic mb-0",
                                ),
                                html.H3(
                                    "Perturbation",
                                    id="group-by-perturbation",
                                    n_clicks=0,
                                    className="group-by-selector",
                                ),
                            ]
                        ),
                        width=2,
                    ),
                    dbc.Col(
                        html.Div(
                            [
                                html.H4(
                                    "Leads to the following",
                                    className="fw-normal fst-italic mb-0",
                                ),
                                html.H3(
                                    "Effect",
                                    id="group-by-effect",
                                    n_clicks=0,
                                    className="group-by-selector",
                                ),
                            ]
                        ),
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
                            id="effect-filter", placeholder="Filter by affected gene"
                        ),
                        width=6,
                    ),
                ],
                className="mb-4",
            ),
            # Data grid
            html.Div(
                id="data-grid",
                children=[
                    dbc.Spinner(color="primary"),
                ],
            ),
        ],
    )


@callback(
    Output("group-by-store", "data"),
    Output("group-by-perturbation", "style"),
    Output("group-by-effect", "style"),
    Input("group-by-perturbation", "n_clicks"),
    Input("group-by-effect", "n_clicks"),
    State("group-by-store", "data"),
)
def update_grouping(pert_clicks, effect_clicks, current_grouping):
    style = {
        "display": "inline-block",
        "border": "1px dotted black",
        "border-radius": "8px",
        "padding": "2px",
        "cursor": "pointer",
    }
    selected_style = {**style, "background-color": "#FFF3CD"}

    ctx = dash.callback_context

    # This handles the user clicking a selector
    if ctx.triggered:
        triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if triggered_id == "group-by-perturbation":
            if current_grouping == "perturbation_gene_name":
                # User clicked the selected one, so unselect
                return None, style, style
            else:
                # User clicked the unselected one, so select
                return "perturbation_gene_name", selected_style, style
        elif triggered_id == "group-by-effect":
            if current_grouping == "effect_gene_name":
                # User clicked the selected one, so unselect
                return None, style, style
            else:
                # User clicked the unselected one, so select
                return "effect_gene_name", style, selected_style

    # This handles the initial page load, applying styles based on the default state
    pert_style = (
        selected_style if current_grouping == "perturbation_gene_name" else style
    )
    eff_style = selected_style if current_grouping == "effect_gene_name" else style
    return dash.no_update, pert_style, eff_style


@callback(
    Output("data-grid", "children"),
    Input("dataset-filter", "value"),
    Input("perturbation-filter", "value"),
    Input("effect-filter", "value"),
    Input("group-by-store", "data"),
)
def update_data_grid(dataset_filter, perturbation_filter, effect_filter, group_by):
    be_url = os.environ.get("PERTURBATION_CATALOGUE_BE", "http://127.0.0.1:8000")

    params = {
        "dataset": dataset_filter,
        "perturbation_gene_name": perturbation_filter,
        "effect_gene_name": effect_filter,
        "group_by": group_by,
    }
    # Filter out None values
    params = {k: v for k, v in params.items() if v}

    try:
        response = requests.get(f"{be_url}/v1/search", params=params)
        response.raise_for_status()
        data = response.json()

        if not data:
            return dbc.Alert("No results found.", color="info")

        return html.Div([render_dataset_row(ds, group_by) for ds in data])

    except requests.exceptions.RequestException as e:
        return dbc.Alert(f"Error connecting to backend: {e}", color="danger")
    except json.JSONDecodeError:
        return dbc.Alert("Error decoding backend response.", color="danger")
