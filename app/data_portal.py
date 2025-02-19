from collections import namedtuple

import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import requests

FilterField = namedtuple("FilterField", ["id", "title"])
filter_fields = [
    FilterField(id="sequenceType", title="Sequence Type"),
    FilterField(id="geneCategory", title="Gene Category"),
    FilterField(id="publicationYear", title="Publication Year"),
]

sortable_columns = {
    "publicationYear": "Publication Year",
    "numVariants": "Number of Variants",
}

main_layout = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5(field.title, className="card-title"),
                                    dbc.RadioItems(
                                        id=field.id,
                                        options=[],
                                        inline=True,
                                        style={"width": "100%"},
                                    ),
                                    dbc.Button(
                                        "Clear",
                                        id=f"clear-{field.id}",
                                        color="success",
                                        size="sm",
                                        style={"margin-top": "10px"},
                                    ),
                                ]
                            )
                        )
                        for field in filter_fields
                    ],
                    width=2,
                    lg=2,
                    md=4,
                    sm=12,
                    xs=12,
                    style={
                        "padding": "30px 0px 0px 20px",
                        "display": "flex",
                        "flexDirection": "column",
                        "gap": "12px",
                    },
                ),
                dbc.Col(
                    [
                        dcc.Input(
                            id="search",
                            type="text",
                            placeholder="Search...",
                            debounce=True,
                            style={
                                "width": "100%",
                                "margin-bottom": "20px",
                            },
                        ),
                        dcc.Store(
                            id="sort-store",
                            data={"field": "publicationYear", "order": "desc"},
                        ),
                        html.Div(
                            id="data-table",
                            style={
                                "overflowY": "auto",
                                "width": "100%",
                            },
                        ),
                        html.Div(
                            [
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
                                html.Label("Items per page:"),
                                dcc.Dropdown(
                                    id="size",
                                    options=[
                                        {"label": str(i), "value": i}
                                        for i in [10, 50, 100, 200]
                                    ],
                                    value=10,
                                    clearable=False,
                                    style={"width": "70px"},
                                ),
                                html.Span(id="pagination-info"),
                            ],
                            style={
                                "display": "flex",
                                "alignItems": "center",
                                "gap": "12px",
                            },
                        ),
                    ],
                    width=10,
                    lg=10,
                    md=8,
                    sm=12,
                    xs=12,
                    style={"padding": "32px 25px 25px 25px"},
                ),
            ],
            className="g-0",
        ),
    ]
)


def create_sort_header(column_name, field_name, current_sort):
    if field_name not in sortable_columns and field_name is not None:
        return html.Th(column_name)

    if field_name is None:
        return html.Th(column_name)

    is_sorted = current_sort["field"] == field_name
    order = current_sort["order"] if is_sorted else None
    arrow = "🡫" if order == "desc" else "🡩" if order == "asc" else ""
    return html.Th(
        html.Div(
            [
                column_name + " " + arrow,
            ],
            style={
                "display": "flex",
                "alignItems": "center",
                "cursor": "pointer",
                "color": "blue" if is_sorted else "inherit",
            },
        ),
        id={"type": "sort-header-container", "field": field_name},
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


def register_callbacks(app):
    @app.callback(
        Output("sort-store", "data"),
        [
            Input({"type": "sort-header-container", "field": dash.ALL}, "n_clicks"),
        ],
        State("sort-store", "data"),
    )
    def update_sort(container_clicks, current_sort):
        if not dash.callback_context.triggered:
            return dash.no_update

        triggered_id = dash.callback_context.triggered[0]["prop_id"]
        if not triggered_id:
            return dash.no_update

        field = eval(triggered_id.split(".")[0])["field"]

        if current_sort["field"] != field:
            return {"field": field, "order": "asc"}
        else:
            if current_sort["order"] == "asc":
                return {"field": field, "order": "desc"}
            else:
                return {"field": field, "order": "asc"}

    @app.callback(
        [
            Output("data-table", "children"),
            *[Output(f.id, "options") for f in filter_fields],
            Output("pagination", "max_value"),
            Output("pagination-info", "children"),
        ],
        [
            Input("search", "value"),
            Input("size", "value"),
            Input("pagination", "active_page"),
            Input("sort-store", "data"),
            *[Input(f.id, "value") for f in filter_fields],
        ],
    )
    def fetch_data(
        q, size, page, sort_data, sequenceType, geneCategory, publicationYear
    ):
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
                "sort_field": sort_data.get("field"),
                "sort_order": sort_data.get("order"),
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

            filter_options = [
                [
                    {"label": b["key"], "value": b["key"]}
                    for b in data.get("aggregations", {})
                    .get(filter_field.id, {})
                    .get("buckets", [])
                ]
                for filter_field in filter_fields
            ]

            pagination_info = f"{start_index} – {end_index} of {total}"

            if data["results"]:
                columns = [
                    ("URN", None),
                    ("Sequence Type", None),
                    ("Gene Name", None),
                    ("Gene Category", None),
                    ("Publication Year", "publicationYear"),
                    ("Number of Variants", "numVariants"),
                ]

                table_header = [
                    html.Thead(
                        html.Tr(
                            [
                                create_sort_header(col[0], col[1], sort_data)
                                for col in columns
                            ]
                        )
                    )
                ]

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
                *filter_options,
                max_pages,
                pagination_info,
            )
        else:
            return html.Div("Error fetching data."), [], [], [], 1, ""

    @app.callback(
        [Output(field.id, "value") for field in filter_fields],
        [Input(f"clear-{field.id}", "n_clicks") for field in filter_fields],
    )
    def clear_filters(*args):
        # Get the ID of the triggered input
        triggered_id = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
        # Initialize the result with no_update for all outputs
        result = [dash.no_update] * len(filter_fields)
        # Find the index of the triggered filter and set its value to None
        for i, field in enumerate(filter_fields):
            if triggered_id == f"clear-{field.id}":
                result[i] = None
                break
        return result


def resolver(url_parts):
    # If there's something after the main URL, return details page.
    if len(url_parts) == 1:
        return details_layout(url_parts[0])
    # Otherwise, return the main data portal page.
    return main_layout
