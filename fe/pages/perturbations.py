import dash
from dash import html, dcc, callback, Input, Output, State, ALL
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

# --- Constants ---
FACET_FIELDS = [
    "tissue",
    "cell_type",
    "cell_line",
    "sex",
    "developmental_stage",
    "disease",
    "library_perturbation_type",
]


# --- Reusable component for rendering labels and values ---
def render_label_value(label, value, value_class="fw-semibold"):
    return html.Div(
        [
            html.Span(f"{label} ", className="fw-light"),
            html.Span(value, className=value_class),
        ]
    )


# --- Data Rendering Functions for Perturb-Seq ---


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
                if value and value != "N/A"
            ],
        ],
        style={"align-self": "start"},
    )


def render_perturbation_cell_perturb_seq(p, is_grouped=False):
    if not is_grouped:
        return html.Div(
            html.Span(p.get("perturbation_gene_name", "N/A"), className="fw-semibold"),
            style={"align-self": "start"},
        )

    gene_name = p.get("perturbation_gene_name", "N/A")
    knockout_text = (
        "Multiple gene knockout" if "|" in gene_name else "Single gene knockout"
    )

    return html.Div(
        [
            html.H3(gene_name, className="fw-bold mb-1"),
            html.Div(knockout_text, className="fw-light"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Span("Affects ", className="fw-light"),
                            html.Span(
                                f'{p.get("n_total", 0)}', className="fw-semibold"
                            ),
                            html.Span(" phenotypes", className="fw-light"),
                        ]
                    ),
                    html.Div(
                        [
                            html.Span("△ ", style={"color": "#2acc06"}),
                            html.Span(
                                f'{p.get("n_up", 0)}',
                                style={"color": "#2acc06"},
                                className="fw-semibold",
                            ),
                            html.Span(" ▽ ", style={"color": "#ff4824"}),
                            html.Span(
                                f'{p.get("n_down", 0)}',
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


def render_change_cell_perturb_seq(c):
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


def render_phenotype_cell_perturb_seq(ph, is_grouped=False):
    if not is_grouped:
        return html.Div(
            html.Span(ph.get("phenotype_gene_name", "N/A"), className="fw-semibold"),
            style={"align-self": "start"},
        )

    base_mean = ph.get("base_mean", 0)
    try:
        base_mean_formatted = f"{int(base_mean):,}"
    except (ValueError, TypeError):
        base_mean_formatted = "N/A"

    return html.Div(
        [
            html.H3(ph.get("phenotype_gene_name", "N/A"), className="fw-bold mb-1"),
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
                                f'{ph.get("n_total", 0)}', className="fw-semibold"
                            ),
                            html.Span(" perturbations", className="fw-light"),
                        ]
                    ),
                    html.Div(
                        [
                            html.Span("△ ", style={"color": "#2acc06"}),
                            html.Span(
                                f'{ph.get("n_up", 0)}',
                                style={"color": "#2acc06"},
                                className="fw-semibold",
                            ),
                            html.Span(" ▽ ", style={"color": "#ff4824"}),
                            html.Span(
                                f'{ph.get("n_down", 0)}',
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


# --- Data Rendering Functions for CRISPR/MAVE ---


def render_perturbation_cell_crispr_mave(p):
    gene_name = p.get("perturbation_gene_name", "N/A")
    pert_name = p.get("perturbation_name")  # MAVE specific
    return html.Div(
        [
            html.Span(gene_name, className="fw-semibold"),
            html.Br(),
            html.Span(pert_name, className="fw-light") if pert_name else None,
        ],
        style={"align-self": "start"},
    )


def render_change_cell_crispr_mave(c):
    score_val = c.get("score_value", 0)
    score_str = f"{score_val:.2f}".replace("-", "−")
    return html.Div(
        [
            html.Span("Score ", className="fw-light"),
            html.Span(score_str, className="fw-semibold"),
        ],
        style={"align-self": "start"},
    )


def render_phenotype_cell_crispr_mave(ph):
    return html.Div(
        html.Span(ph.get("score_name", "N/A"), className="fw-semibold"),
        style={"align-self": "start"},
    )


# --- Truncation Notice ---
def render_truncation_notice(text, grid_column):
    return html.Div(
        html.I(text),
        style={
            "grid-column": grid_column,
            "text-align": "center",
            "background-color": "#f8f9fa",
            "padding": "0.25rem",
            "border-radius": "0.25rem",
            "margin-top": "0.5rem",
        },
    )


# --- Page Layout ---
def layout():
    return html.Div(
        [
            dcc.Store(id="api-response-store"),
            dcc.Store(id="group-by-store", data="perturbation_gene_name"),
            dcc.Store(id="change-direction-store"),
            html.Div(
                className="text-center",
                style={
                    "background-image": 'url("https://new-pc-fe-959149465821.europe-west2.run.app/assets/banner.png")',
                    "background-size": "cover",
                    "background-position": "center",
                    "padding": "4.5rem 2rem",
                    "color": "white",
                },
                children=[
                    html.H1(
                        "Perturbation Catalogue",
                        className="display-4 mb-3 fw-bold",
                        style={
                            "font-size": "60px",
                            "text-shadow": "2px 2px 4px rgba(0, 0, 0, 0.6)",
                        },
                    ),
                    html.P(
                        "A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.",
                        className="lead",
                        style={
                            "font-size": "25px",
                            "text-shadow": "1px 1px 3px rgba(0, 0, 0, 0.6)",
                        },
                    ),
                ],
            ),
            dbc.Container(
                fluid=True,
                className="px-5 pt-5 mb-5",
                children=[
                    dbc.Row(
                        [
                            # --- Left Column: Facet Filters ---
                            dbc.Col(
                                width=3,
                                children=[
                                    html.H4("Filter by", className="mb-3"),
                                    html.Div(
                                        id="facet-filters-container",
                                        children=[
                                            html.Div(
                                                [
                                                    html.H5(
                                                        field.replace("_", " ").title(),
                                                        className="mt-3",
                                                    ),
                                                    dbc.Checklist(
                                                        options=[],
                                                        id={
                                                            "type": "facet-filter",
                                                            "index": field,
                                                        },
                                                    ),
                                                ],
                                                id=f"facet-wrapper-{field}",
                                                style={"display": "none"},
                                            )
                                            for field in FACET_FIELDS
                                        ],
                                    ),
                                ],
                            ),
                            # --- Right Column: Controls and Data ---
                            dbc.Col(
                                width=9,
                                children=[
                                    # --- Controls ---
                                    html.Div(
                                        className="d-flex justify-content-center align-items-center mb-4",
                                        children=[
                                            html.Span(
                                                "Group by ",
                                                className="align-self-center me-2",
                                                style={"font-size": "1.25rem"},
                                            ),
                                            dbc.ButtonGroup(
                                                [
                                                    dbc.Button(
                                                        "Perturbation",
                                                        id="group-by-perturbation-btn",
                                                        color="primary",
                                                        n_clicks=0,
                                                    ),
                                                    dbc.Button(
                                                        "Phenotype",
                                                        id="group-by-phenotype-btn",
                                                        color="light",
                                                        n_clicks=0,
                                                    ),
                                                ],
                                                size="lg",
                                            ),
                                            html.Div(className="flex-grow-1"),  # Spacer
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupText("Max Datasets"),
                                                    dbc.Input(
                                                        id="limit-datasets",
                                                        type="number",
                                                        value=10,
                                                        min=1,
                                                    ),
                                                ],
                                                className="ms-3",
                                                style={"width": "200px"},
                                            ),
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupText("Max Groups"),
                                                    dbc.Input(
                                                        id="limit-groups",
                                                        type="number",
                                                        value=10,
                                                        min=1,
                                                    ),
                                                ],
                                                className="ms-3",
                                                style={"width": "200px"},
                                            ),
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupText("Max Rows"),
                                                    dbc.Input(
                                                        id="limit-rows",
                                                        type="number",
                                                        value=10,
                                                        min=1,
                                                    ),
                                                ],
                                                className="ms-3",
                                                style={"width": "200px"},
                                            ),
                                        ],
                                    ),
                                    # --- Results Summary ---
                                    html.Div(
                                        id="results-summary-container",
                                        className="text-center mb-4",
                                    ),
                                    # --- CSS Grid for Headers, Filters, and Data ---
                                    html.Div(
                                        style={
                                            "display": "grid",
                                            "grid-template-columns": "4fr 3fr 3fr 3fr",
                                            "row-gap": "0.25rem",
                                            "column-gap": "1.5rem",
                                            "grid-auto-rows": "min-content",
                                        },
                                        children=[
                                            # --- Headers ---
                                            html.Div(
                                                [
                                                    html.H4(
                                                        "According to",
                                                        className="fw-normal fst-italic mb-0",
                                                    ),
                                                    html.H3("Dataset"),
                                                ]
                                            ),
                                            html.Div(
                                                [
                                                    html.H4(
                                                        "Introducing",
                                                        className="fw-normal fst-italic mb-0",
                                                    ),
                                                    html.H3("Perturbation"),
                                                ]
                                            ),
                                            html.Div(
                                                [
                                                    html.H4(
                                                        "Leads to",
                                                        className="fw-normal fst-italic mb-0",
                                                    ),
                                                    html.H3("Change"),
                                                ]
                                            ),
                                            html.Div(
                                                [
                                                    html.H4(
                                                        "Affecting",
                                                        className="fw-normal fst-italic mb-0",
                                                    ),
                                                    html.H3("Phenotype"),
                                                ]
                                            ),
                                            # --- Filter Fields ---
                                            dbc.Input(
                                                id="filter-dataset",
                                                placeholder="Filter by dataset metadata",
                                            ),
                                            dbc.Input(
                                                id="filter-perturbation",
                                                placeholder="Filter by perturbed gene",
                                            ),
                                            dbc.ButtonGroup(
                                                [
                                                    dbc.Button(
                                                        "△ Up",
                                                        id="filter-change-up-btn",
                                                        color="light",
                                                        n_clicks=0,
                                                    ),
                                                    dbc.Button(
                                                        "▽ Down",
                                                        id="filter-change-down-btn",
                                                        color="light",
                                                        n_clicks=0,
                                                    ),
                                                    dbc.Button(
                                                        "△▽ Both",
                                                        id="filter-change-both-btn",
                                                        color="secondary",
                                                        n_clicks=0,
                                                    ),
                                                ],
                                                id="filter-change-group",
                                            ),
                                            dbc.Input(
                                                id="filter-phenotype",
                                                placeholder="Filter by phenotype gene",
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
                                                style={
                                                    "grid-column": "1 / -1",
                                                    "display": "contents",
                                                },
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                        ]
                    )
                ],
            ),
        ]
    )


# --- Callbacks ---


@callback(
    Output("group-by-store", "data"),
    Output("group-by-perturbation-btn", "color"),
    Output("group-by-phenotype-btn", "color"),
    Input("group-by-perturbation-btn", "n_clicks"),
    Input("group-by-phenotype-btn", "n_clicks"),
)
def update_group_by_selection(pert_clicks, pheno_clicks):
    ctx = dash.callback_context
    if not ctx.triggered or (pert_clicks == 0 and pheno_clicks == 0):
        return "perturbation_gene_name", "primary", "light"
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if button_id == "group-by-perturbation-btn":
        return "perturbation_gene_name", "primary", "light"
    else:
        return "phenotype_gene_name", "light", "primary"


@callback(
    Output("change-direction-store", "data"),
    Output("filter-change-up-btn", "color"),
    Output("filter-change-down-btn", "color"),
    Output("filter-change-both-btn", "color"),
    Input("filter-change-up-btn", "n_clicks"),
    Input("filter-change-down-btn", "n_clicks"),
    Input("filter-change-both-btn", "n_clicks"),
)
def update_change_filter_selection(up_clicks, down_clicks, both_clicks):
    ctx = dash.callback_context
    if not ctx.triggered or (up_clicks == 0 and down_clicks == 0 and both_clicks == 0):
        return None, "light", "light", "secondary"
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if button_id == "filter-change-up-btn":
        return "increased", "secondary", "light", "light"
    elif button_id == "filter-change-down-btn":
        return "decreased", "light", "secondary", "light"
    else:  # both
        return None, "light", "light", "secondary"


@callback(
    Output("api-response-store", "data"),
    Input("filter-dataset", "value"),
    Input("filter-perturbation", "value"),
    Input("filter-phenotype", "value"),
    Input("group-by-store", "data"),
    Input("change-direction-store", "data"),
    Input("limit-datasets", "value"),
    Input("limit-groups", "value"),
    Input("limit-rows", "value"),
    [
        Input({"type": "facet-filter", "index": field}, "value")
        for field in FACET_FIELDS
    ],
)
def fetch_data_from_be(
    dataset_filter,
    pert_filter,
    pheno_filter,
    group_by,
    change_direction,
    limit_datasets,
    limit_groups,
    limit_rows,
    *facet_values,
):
    backend_url = os.getenv("PERTURBATION_CATALOGUE_BE", "http://127.0.0.1:8000")
    api_endpoint = f"{backend_url}/v1/search"

    params = {
        "group_by": group_by,
        "max_datasets_per_modality": limit_datasets,
        "max_top_level": limit_groups,
        "max_rows": limit_rows,
    }
    if dataset_filter:
        params["dataset_metadata"] = dataset_filter
    if pert_filter:
        params["perturbation_gene_name"] = pert_filter
    if pheno_filter:
        params["phenotype_gene_name"] = pheno_filter
    if change_direction:
        params["change_direction"] = change_direction

    for i, field in enumerate(FACET_FIELDS):
        if facet_values[i]:
            params[field] = ",".join(facet_values[i])

    try:
        response = requests.get(api_endpoint, params=params)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        return {"error": f"Error fetching data: {e}"}
    except json.JSONDecodeError:
        return {"error": "Invalid JSON response from server."}


@callback(
    [Output(f"facet-wrapper-{field}", "style") for field in FACET_FIELDS]
    + [
        Output({"type": "facet-filter", "index": field}, "options")
        for field in FACET_FIELDS
    ],
    Input("api-response-store", "data"),
)
def update_facet_filters(data):
    if not data or "facet_counts" not in data:
        return [{"display": "none"}] * len(FACET_FIELDS) + [[]] * len(FACET_FIELDS)

    facet_counts = data["facet_counts"]
    styles = []
    options_list = []
    for field in FACET_FIELDS:
        if field in facet_counts and facet_counts[field]:
            styles.append({"display": "block"})
            options = [
                {"label": f'{item["value"]} ({item["count"]})', "value": item["value"]}
                for item in facet_counts[field]
            ]
            options_list.append(options)
        else:
            styles.append({"display": "none"})
            options_list.append([])
    return styles + options_list


@callback(
    Output("results-summary-container", "children"), Input("api-response-store", "data")
)
def update_results_summary(data):
    if not data or "modalities" not in data:
        return []

    summary_parts = []
    for modality in data["modalities"]:
        count = modality.get("total_datasets_count", 0)
        name = modality.get("modality", "N/A")
        if count > 0:
            summary_parts.append(f"{count} {name} datasets")

    return (
        html.Span(f"Found: {', '.join(summary_parts)}")
        if summary_parts
        else "No results found."
    )


@callback(
    Output("data-grid-container", "children"),
    Input("api-response-store", "data"),
    State("group-by-store", "data"),
    State("limit-datasets", "value"),
    State("limit-groups", "value"),
    State("limit-rows", "value"),
)
def update_data_grid(data, group_by, limit_datasets, limit_groups, limit_rows):
    if not data:
        return []
    if "error" in data:
        return [
            html.Div(
                data["error"],
                style={
                    "grid-column": "1 / -1",
                    "text-align": "center",
                    "padding": "2rem",
                },
            )
        ]
    if not data.get("modalities"):
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

    for modality in data["modalities"]:
        modality_name = modality["modality"]
        total_datasets = modality["total_datasets_count"]
        datasets = modality.get("datasets", [])

        if not datasets:
            continue

        grid_cells.append(
            html.Div(
                html.H2(modality_name.replace("-", " ").title()),
                style={"grid-column": "1 / -1", "margin-top": "2rem"},
            )
        )

        for item_index, item in enumerate(datasets):
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

            if modality_name == "perturb-seq":
                render_perturb_seq_dataset(
                    grid_cells, item, group_by, limit_groups, limit_rows
                )
            else:  # crispr-screen or mave
                render_crispr_mave_dataset(grid_cells, item, modality_name, limit_rows)

        if len(datasets) < total_datasets:
            grid_cells.append(
                render_truncation_notice(
                    f"Displaying top {len(datasets)} of {total_datasets} datasets",
                    "1 / -1",
                )
            )

    return (
        grid_cells
        if grid_cells
        else [
            html.Div(
                "No results found for the selected modalities.",
                style={
                    "grid-column": "1 / -1",
                    "text-align": "center",
                    "padding": "2rem",
                },
            )
        ]
    )


def render_perturb_seq_dataset(grid_cells, item, group_by, limit_groups, limit_rows):
    group_key, sub_item_key = (
        ("by_perturbation", "change_phenotype")
        if group_by == "perturbation_gene_name"
        else ("by_phenotype", "perturbation_change")
    )

    group_items = item.get(group_key, [])
    total_groups = 0  # BE doesn't provide this yet

    row_spans = []
    for g in group_items:
        sub_items_len = len(g.get(sub_item_key, []))
        row_span = max(1, sub_items_len)
        if (
            sub_items_len == limit_rows
            and g.get("perturbation", {}).get("n_total", 0) > limit_rows
        ):
            row_span += 1  # For truncation notice
        row_spans.append(row_span)

    dataset_rowspan = sum(row_spans) + max(0, len(group_items) - 1)
    if (
        len(group_items) == limit_groups
    ):  # and total_groups > limit_groups: # BE doesn't provide total_groups
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

        sub_items = group_item.get(sub_item_key, [])

        if group_by == "perturbation_gene_name":
            pert = group_item.get("perturbation", {})
            cell = render_perturbation_cell_perturb_seq(pert, is_grouped=True)
            total_sub_items = pert.get("n_total", 0)
            cell.style["grid-row"] = f"span {row_spans[i]}"
            grid_cells.append(cell)

            if not sub_items:
                grid_cells.extend([html.Div("-"), html.Div("-")])
            else:
                for j, sub_item in enumerate(sub_items):
                    if j > 0:
                        grid_cells.append(
                            html.Div()
                        )  # Empty cell for the grouped perturbation
                    grid_cells.append(
                        render_change_cell_perturb_seq(sub_item.get("change", {}))
                    )
                    grid_cells.append(
                        render_phenotype_cell_perturb_seq(sub_item.get("phenotype", {}))
                    )
        else:  # group by phenotype
            pheno = group_item.get("phenotype", {})
            cell = render_phenotype_cell_perturb_seq(pheno, is_grouped=True)
            total_sub_items = pheno.get("n_total", 0)
            cell.style["grid-row"] = f"span {row_spans[i]}"

            if not sub_items:
                grid_cells.extend([html.Div("-"), html.Div("-"), cell])
            else:
                for j, sub_item in enumerate(sub_items):
                    if j > 0:
                        grid_cells.extend(
                            [html.Div(), html.Div()]
                        )  # Empty cells for change and grouped phenotype
                    grid_cells.append(
                        render_perturbation_cell_perturb_seq(
                            sub_item.get("perturbation", {})
                        )
                    )
                    grid_cells.append(
                        render_change_cell_perturb_seq(sub_item.get("change", {}))
                    )
                    if j == 0:
                        grid_cells.append(cell)

        if len(sub_items) == limit_rows and total_sub_items > limit_rows:
            grid_cells.append(
                render_truncation_notice(
                    f"Displaying top {limit_rows} of {total_sub_items} rows",
                    "3 / 5" if group_by == "perturbation_gene_name" else "2 / 4",
                )
            )

    if len(group_items) == limit_groups:  # and total_groups > limit_groups:
        grid_cells.append(
            render_truncation_notice(f"Displaying top {limit_groups} groups", "2 / 5")
        )


def render_crispr_mave_dataset(grid_cells, item, modality_name, limit_rows):
    data_rows = item.get("data", [])
    total_rows = 0  # BE doesn't provide this

    dataset_rowspan = max(1, len(data_rows))
    if len(data_rows) == limit_rows:  # and total_rows > limit_rows:
        dataset_rowspan += 1

    dataset_cell = render_dataset_cell(item)
    dataset_cell.style["grid-row"] = f"span {dataset_rowspan}"
    grid_cells.append(dataset_cell)

    if not data_rows:
        grid_cells.extend([html.Div("-"), html.Div("-"), html.Div("-")])
    else:
        for i, row in enumerate(data_rows):
            if i > 0:
                grid_cells.append(html.Div())  # Empty cell for dataset
            grid_cells.append(
                render_perturbation_cell_crispr_mave(row.get("perturbation", {}))
            )
            grid_cells.append(render_change_cell_crispr_mave(row.get("change", {})))
            grid_cells.append(
                render_phenotype_cell_crispr_mave(row.get("phenotype", {}))
            )

    if len(data_rows) == limit_rows:  # and total_rows > limit_rows:
        grid_cells.append(
            render_truncation_notice(f"Displaying top {limit_rows} rows", "2 / 5")
        )
