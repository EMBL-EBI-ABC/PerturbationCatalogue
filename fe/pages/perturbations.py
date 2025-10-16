import dash
from dash import html, dcc, callback, Input, Output, State
import dash_bootstrap_components as dbc
import requests
import os
import json

dash.register_page(
    __name__,
    path="/perturbations",
    relative_path="/perturbations",
    name="Perturbations",
    button="Explore perturbations",
    description="A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.",
    icon="bi-tropical-storm",
)


# --- Reusable component for rendering labels and values ---
def render_label_value(label, value, value_class="fw-semibold"):
    return html.Div(
        [
            html.Span(f"{label} ", className="fw-light"),
            html.Span(value, className=value_class),
        ]
    )


# --- Data Rendering Functions ---


def render_dataset_cell(item):
    dataset = item.get("dataset", {})
    dataset_id = dataset.get("dataset_id", "Unknown Dataset")
    formatted_dataset_id = dataset_id.replace("_", " ").capitalize()
    metadata = [
        ("Tissue", dataset.get("tissue", "N/A")),
        ("Cell type", dataset.get("cell_type", "N/A")),
        ("Cell line", dataset.get("cell_line", "N/A")),
        ("Sex", dataset.get("sex", "N/A")),
        ("Developmental stage", dataset.get("developmental_stage", "N/A")),
        ("Disease", dataset.get("disease", "N/A")),
        ("Library perturbation type", dataset.get("library_perturbation_type", "N/A")),
    ]
    return html.Div(
        [
            html.H3(formatted_dataset_id, className="fw-bold"),
            *[
                render_label_value(label, str(value).capitalize())
                for label, value in metadata
                if value != "N/A"
            ],
        ],
        style={"align-self": "start"},
    )


def render_perturbation_cell(p, is_grouped=False):
    if not is_grouped:
        return html.Div(
            html.Span(p.get("perturbation_gene_name", "N/A"), className="fw-semibold"),
            style={"align-self": "start"},
        )

    # Grouped rendering
    return html.Div(
        [
            html.H3(p.get("perturbation_gene_name", "N/A"), className="fw-bold"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Span("Affects ", className="fw-light"),
                            html.Span(
                                f"{p.get('n_total', 0)}", className="fw-semibold"
                            ),
                            html.Span(" phenotypes", className="fw-light"),
                        ]
                    ),
                    html.Div(
                        [
                            html.Span("△ ", style={"color": "#2acc06"}),
                            html.Span(
                                f"{p.get('n_up', 0)}",
                                style={"color": "#2acc06"},
                                className="fw-semibold",
                            ),
                            html.Span(" ▽ ", style={"color": "#ff4824"}),
                            html.Span(
                                f"{p.get('n_down', 0)}",
                                style={"color": "#ff4824"},
                                className="fw-semibold",
                            ),
                        ]
                    ),
                ]
            ),
        ],
        style={"align-self": "start"},
    )


def render_change_cell(c):
    log2fc_val = c.get("log2fc", 0)
    padj_val = c.get("padj", 0)

    log2fc_str = f"{log2fc_val:.2f}".replace("-", "−")
    padj_str = f"{padj_val:.2e}".replace("e-", "e−")

    direction = c.get("direction", "increased")
    arrow_char, color = (
        ("△", "#2acc06") if direction == "increased" else ("▽", "#ff4824")
    )

    return html.Div(
        [
            html.Div(
                [
                    html.Span("padj ", className="fw-light"),
                    html.Span(padj_str, className="fw-semibold"),
                ],
                style={"min-width": "110px"},
            ),
            html.Div(
                [
                    html.Span("log2fc ", className="fw-light"),
                    html.Span(log2fc_str, className="fw-semibold"),
                ],
                style={"min-width": "90px"},
            ),
            html.Div(
                html.Span(arrow_char, style={"color": color, "font-weight": "bold"}),
                style={"min-width": "30px"},
            ),
        ],
        style={
            "align-self": "start",
            "display": "flex",
            "justify-content": "space-between",
            "width": "100%",
        },
    )


def render_phenotype_cell(ph, is_grouped=False):
    if not is_grouped:
        return html.Div(
            html.Span(ph.get("phenotype_gene_name", "N/A"), className="fw-semibold"),
            style={"align-self": "start"},
        )

    # Grouped rendering
    base_mean = ph.get("base_mean", 0)
    try:
        base_mean_formatted = f"{int(base_mean):,}"
    except (ValueError, TypeError):
        base_mean_formatted = "N/A"

    return html.Div(
        [
            html.H3(ph.get("phenotype_gene_name", "N/A"), className="fw-bold"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Span("Base expression ", className="fw-light"),
                            html.Span(base_mean_formatted, className="fw-semibold"),
                        ]
                    ),
                    html.Div(
                        [
                            html.Span("Affected by ", className="fw-light"),
                            html.Span(
                                f"{ph.get('n_total', 0)}", className="fw-semibold"
                            ),
                            html.Span(" perturbations", className="fw-light"),
                        ]
                    ),
                    html.Div(
                        [
                            html.Span("△ ", style={"color": "#2acc06"}),
                            html.Span(
                                f"{ph.get('n_up', 0)}",
                                style={"color": "#2acc06"},
                                className="fw-semibold",
                            ),
                            html.Span(" ▽ ", style={"color": "#ff4824"}),
                            html.Span(
                                f"{ph.get('n_down', 0)}",
                                style={"color": "#ff4824"},
                                className="fw-semibold",
                            ),
                        ]
                    ),
                ]
            ),
        ],
        style={"align-self": "start"},
    )


# --- Page Layout ---


def layout():
    return dbc.Container(
        fluid=True,
        className="px-5",
        children=[
            html.H1(
                "Perturbation Catalogue",
                className="text-center display-4 mb-3",
                style={"font-size": "60px", "margin-top": "4.5rem"},
            ),
            html.P(
                "A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.",
                className="text-center lead",
                style={"font-size": "25px", "margin-bottom": "4.5rem"},
            ),
            html.Div(
                className="d-flex justify-content-center mb-5",
                children=[
                    html.Span("Group by ", className="align-self-center me-2"),
                    dbc.ButtonGroup(
                        [
                            dbc.Button(
                                "Perturbation",
                                id="group-by-perturbation-btn",
                                color="primary",
                            ),
                            dbc.Button(
                                "Phenotype", id="group-by-phenotype-btn", color="light"
                            ),
                        ]
                    ),
                ],
            ),
            # --- CSS Grid for Headers, Filters, and Data ---
            html.Div(
                className="mb-5",
                style={
                    "display": "grid",
                    "grid-template-columns": "4fr 3fr 3fr 3fr",
                    "row-gap": "0.5rem",
                    "column-gap": "1rem",
                },
                children=[
                    # --- Headers ---
                    html.Div(
                        [
                            html.H4(
                                "According to", className="fw-normal fst-italic mb-0"
                            ),
                            html.H3("Dataset"),
                        ]
                    ),
                    html.Div(
                        [
                            html.H4(
                                "Introducing", className="fw-normal fst-italic mb-0"
                            ),
                            html.H3("Perturbation"),
                        ]
                    ),
                    html.Div(
                        [
                            html.H4("Leads to", className="fw-normal fst-italic mb-0"),
                            html.H3("Change"),
                        ]
                    ),
                    html.Div(
                        [
                            html.H4("Affecting", className="fw-normal fst-italic mb-0"),
                            html.H3("Phenotype"),
                        ]
                    ),
                    # --- Filter Fields ---
                    dbc.Input(
                        id="filter-dataset", placeholder="Filter by dataset metadata"
                    ),
                    dbc.Input(
                        id="filter-perturbation", placeholder="Filter by perturbed gene"
                    ),
                    dbc.ButtonGroup(
                        [
                            dbc.Button(
                                "△ Up", id="filter-change-up-btn", color="light"
                            ),
                            dbc.Button(
                                "▽ Down", id="filter-change-down-btn", color="light"
                            ),
                            dbc.Button(
                                "△▽ Both", id="filter-change-both-btn", color="primary"
                            ),
                        ],
                        id="filter-change-group",
                    ),
                    dbc.Input(
                        id="filter-phenotype", placeholder="Filter by phenotype gene"
                    ),
                    # --- Data Grid Container ---
                    html.Div(
                        style={
                            "grid-column": "1 / -1",
                            "border-top": "2px solid #dee2e6",
                            "margin-top": "2rem",
                            "padding-top": "2rem",
                        }
                    ),
                    html.Div(
                        id="data-grid-container",
                        style={"grid-column": "1 / -1", "display": "contents"},
                    ),
                ],
            ),
        ],
    )


# --- Callbacks ---


@callback(
    [
        Output("group-by-perturbation-btn", "color"),
        Output("group-by-phenotype-btn", "color"),
    ],
    [
        Input("group-by-perturbation-btn", "n_clicks"),
        Input("group-by-phenotype-btn", "n_clicks"),
    ],
)
def update_group_by_selection(pert_clicks, pheno_clicks):
    ctx = dash.callback_context
    if not ctx.triggered:
        return "primary", "light"
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "group-by-perturbation-btn":
        return "primary", "light"
    else:
        return "light", "primary"


@callback(
    [
        Output("filter-change-up-btn", "color"),
        Output("filter-change-down-btn", "color"),
        Output("filter-change-both-btn", "color"),
    ],
    [
        Input("filter-change-up-btn", "n_clicks"),
        Input("filter-change-down-btn", "n_clicks"),
        Input("filter-change-both-btn", "n_clicks"),
    ],
)
def update_change_filter_selection(up_clicks, down_clicks, both_clicks):
    ctx = dash.callback_context
    if not ctx.triggered:
        return "light", "light", "primary"
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "filter-change-up-btn":
        return "primary", "light", "light"
    elif button_id == "filter-change-down-btn":
        return "light", "primary", "light"
    else:  # both
        return "light", "light", "primary"


@callback(
    Output("data-grid-container", "children"),
    [
        Input("filter-dataset", "value"),
        Input("filter-perturbation", "value"),
        Input("filter-phenotype", "value"),
        Input("group-by-perturbation-btn", "color"),
        Input("filter-change-up-btn", "color"),
        Input("filter-change-down-btn", "color"),
    ],
)
def update_data_grid(
    dataset_filter,
    pert_filter,
    pheno_filter,
    group_by_pert_color,
    change_up_color,
    change_down_color,
):
    backend_url = os.getenv("PERTURBATION_CATALOGUE_BE", "http://127.0.0.1:8000")
    api_endpoint = f"{backend_url}/v1/search"

    group_by = (
        "perturbation_gene_name"
        if group_by_pert_color == "primary"
        else "phenotype_gene_name"
    )

    params = {"group_by": group_by}
    if dataset_filter:
        params["dataset_metadata"] = dataset_filter
    if pert_filter:
        params["perturbation_gene_name"] = pert_filter
    if pheno_filter:
        params["phenotype_gene_name"] = pheno_filter

    if change_up_color == "primary":
        params["change_direction"] = "increased"
    elif change_down_color == "primary":
        params["change_direction"] = "decreased"

    try:
        response = requests.get(api_endpoint, params=params)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        return [
            html.Div(
                f"Error fetching data: {e}",
                style={
                    "grid-column": "1 / -1",
                    "text-align": "center",
                    "padding": "2rem",
                },
            )
        ]
    except json.JSONDecodeError:
        return [
            html.Div(
                "Error: Invalid JSON response from server.",
                style={
                    "grid-column": "1 / -1",
                    "text-align": "center",
                    "padding": "2rem",
                },
            )
        ]

    if not data:
        return [
            html.Div(
                "No results found.",
                style={
                    "grid-column": "1 / -1",
                    "text-align": "center",
                    "padding": "2rem",
                },
            )
        ]

    grid_cells = []
    for item_index, item in enumerate(data):
        if item_index > 0:
            grid_cells.append(
                html.Div(
                    style={
                        "grid-column": "1 / -1",
                        "border-top": "2px solid #dee2e6",
                        "margin-top": "2rem",
                        "padding-top": "2rem",
                    }
                )
            )

        if group_by == "perturbation_gene_name":
            group_items = item.get("by_perturbation", [])

            group_rowspans = []
            for g in group_items:
                sub_items_len = len(g.get("change_phenotype", []))
                row_span = max(1, sub_items_len)
                if sub_items_len == 10:
                    row_span += 1
                group_rowspans.append(row_span)

            dataset_rowspan = sum(group_rowspans) + max(0, len(group_items) - 1)
            if len(group_items) == 3:
                dataset_rowspan += 1

            dataset_cell = render_dataset_cell(item)
            if dataset_rowspan > 0:
                dataset_cell.style["grid-row"] = f"span {dataset_rowspan}"
            grid_cells.append(dataset_cell)

            for i, group_item in enumerate(group_items):
                if i > 0:
                    grid_cells.append(
                        html.Div(
                            style={
                                "grid-column": "2 / 5",
                                "border-top": "1px solid #dee2e6",
                                "margin-top": "1rem",
                                "margin-bottom": "1rem",
                            }
                        )
                    )

                pert_cell = render_perturbation_cell(
                    group_item.get("perturbation", {}), is_grouped=True
                )
                sub_items = group_item.get("change_phenotype", [])
                pert_rowspan = group_rowspans[i]
                if pert_rowspan > 0:
                    pert_cell.style["grid-row"] = f"span {pert_rowspan}"

                grid_cells.append(pert_cell)

                if not sub_items:
                    grid_cells.extend(
                        [
                            html.Div(),
                            html.Div(),
                        ]
                    )
                else:
                    grid_cells.append(
                        render_change_cell(sub_items[0].get("change", {}))
                    )
                    grid_cells.append(
                        render_phenotype_cell(
                            sub_items[0].get("phenotype", {}), is_grouped=False
                        )
                    )

                    for sub_item in sub_items[1:]:
                        grid_cells.append(
                            render_change_cell(sub_item.get("change", {}))
                        )
                        grid_cells.append(
                            render_phenotype_cell(
                                sub_item.get("phenotype", {}), is_grouped=False
                            )
                        )

                if len(sub_items) == 10:
                    grid_cells.append(
                        html.Div(
                            [
                                html.Span("Displaying first 10 entries "),
                                dcc.Link("View all", href="#"),
                            ],
                            style={
                                "grid-column": "3 / 5",
                                "text-align": "left",
                                "background-color": "#f8f9fa",
                                "padding": "0.25rem",
                                "border-radius": "0.25rem",
                            },
                        )
                    )

            if len(group_items) == 3:
                grid_cells.append(
                    html.Div(
                        [
                            html.Span("Displaying first 3 groups "),
                            dcc.Link("View all", href="#"),
                        ],
                        style={
                            "grid-column": "2 / 5",
                            "text-align": "left",
                            "background-color": "#f8f9fa",
                            "padding": "0.25rem",
                            "border-radius": "0.25rem",
                            "margin-top": "0.5rem",
                        },
                    )
                )

        else:  # group by phenotype
            group_items = item.get("by_phenotype", [])

            group_rowspans = []
            for g in group_items:
                sub_items_len = len(g.get("perturbation_change", []))
                row_span = max(1, sub_items_len)
                if sub_items_len == 10:
                    row_span += 1
                group_rowspans.append(row_span)

            dataset_rowspan = sum(group_rowspans) + max(0, len(group_items) - 1)
            if len(group_items) == 3:
                dataset_rowspan += 1

            dataset_cell = render_dataset_cell(item)
            if dataset_rowspan > 0:
                dataset_cell.style["grid-row"] = f"span {dataset_rowspan}"
            grid_cells.append(dataset_cell)

            for i, group_item in enumerate(group_items):
                if i > 0:
                    grid_cells.append(
                        html.Div(
                            style={
                                "grid-column": "2 / 5",
                                "border-top": "1px solid #dee2e6",
                                "margin-top": "1rem",
                                "margin-bottom": "1rem",
                            }
                        )
                    )

                pheno_cell = render_phenotype_cell(
                    group_item.get("phenotype", {}), is_grouped=True
                )
                sub_items = group_item.get("perturbation_change", [])
                pheno_rowspan = group_rowspans[i]
                if pheno_rowspan > 0:
                    pheno_cell.style["grid-row"] = f"span {pheno_rowspan}"

                if not sub_items:
                    grid_cells.extend([html.Div(), html.Div(), pheno_cell])
                else:
                    grid_cells.append(
                        render_perturbation_cell(
                            sub_items[0].get("perturbation", {}), is_grouped=False
                        )
                    )
                    grid_cells.append(
                        render_change_cell(sub_items[0].get("change", {}))
                    )
                    grid_cells.append(pheno_cell)

                    for sub_item in sub_items[1:]:
                        grid_cells.append(
                            render_perturbation_cell(
                                sub_item.get("perturbation", {}), is_grouped=False
                            )
                        )
                        grid_cells.append(
                            render_change_cell(sub_item.get("change", {}))
                        )

                if len(sub_items) == 10:
                    grid_cells.append(html.Div())
                    grid_cells.append(
                        html.Div(
                            [
                                html.Span("Displaying first 10 entries "),
                                dcc.Link("View all", href="#"),
                            ],
                            style={
                                "grid-column": "2 / 4",
                                "text-align": "left",
                                "background-color": "#f8f9fa",
                                "padding": "0.25rem",
                                "border-radius": "0.25rem",
                            },
                        )
                    )

            if len(group_items) == 3:
                grid_cells.append(
                    html.Div(
                        [
                            html.Span("Displaying first 3 groups "),
                            dcc.Link("View all", href="#"),
                        ],
                        style={
                            "grid-column": "2 / 5",
                            "text-align": "left",
                            "background-color": "#f8f9fa",
                            "padding": "0.25rem",
                            "border-radius": "0.25rem",
                            "margin-top": "0.5rem",
                        },
                    )
                )

    return grid_cells
