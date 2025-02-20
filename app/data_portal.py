from collections import namedtuple

import dash
from dash import html, Input, Output, State
import dash_bootstrap_components as dbc
import requests

from elastic_table import details, table

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


def create_table_header(column_name, field_name, current_sort):
    if field_name not in sortable_columns and field_name is not None:
        return html.Th(column_name)

    is_sorted = current_sort.get("field") == field_name
    order = current_sort.get("order") if is_sorted else None
    arrow = "🡫" if order == "desc" else "🡩" if order == "asc" else ""

    return (
        html.Th(
            html.Div(
                [f"{column_name} {arrow}"],
                className="d-flex align-items-center",
                style={
                    "cursor": "pointer",
                    "color": "blue" if is_sorted else "inherit",
                },
            ),
            id={"type": "sort-header-container", "field": field_name},
        )
        if field_name is not None
        else html.Th(column_name)
    )


def create_table(data, sort_data):
    if not data:
        return html.Div("No data found.")

    columns = [
        ("URN", None),
        ("Sequence Type", None),
        ("Gene Name", None),
        ("Gene Category", None),
        ("Publication Year", "publicationYear"),
        ("Number of Variants", "numVariants"),
    ]

    rows = [
        html.Tr(
            [
                html.Td(
                    html.A(
                        row["urn"],
                        href=f"/data-portal/{row['urn']}",
                        className="text-decoration-none text-nowrap",
                    )
                ),
                html.Td(row.get("sequenceType", "N/A")),
                html.Td(row.get("geneName", "N/A")),
                html.Td(row.get("geneCategory", "N/A")),
                html.Td(row.get("publicationYear", "N/A")),
                html.Td(row.get("numVariants", "N/A")),
            ]
        )
        for row in data
    ]

    return dbc.Table(
        [
            html.Thead(
                html.Tr(
                    [create_table_header(col[0], col[1], sort_data) for col in columns]
                )
            ),
            html.Tbody(rows),
        ],
        bordered=True,
        hover=True,
        responsive=True,
    )


def register_callbacks(app):
    @app.callback(
        Output("sort-store", "data"),
        Input({"type": "sort-header-container", "field": dash.ALL}, "n_clicks"),
        State("sort-store", "data"),
    )
    def update_sort(_, current_sort):
        if not dash.callback_context.triggered:
            return dash.no_update

        field = eval(dash.callback_context.triggered[0]["prop_id"].split(".")[0])[
            "field"
        ]
        if current_sort["field"] != field:
            return {"field": field, "order": "asc"}
        return {
            "field": field,
            "order": "desc" if current_sort["order"] == "asc" else "asc",
        }

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
    def fetch_data(q, size, page, sort_data, *filter_values):
        params = {
            "q": q,
            "size": size,
            "start": (page - 1) * size if page else 0,
            **{f.id: v for f, v in zip(filter_fields, filter_values) if v},
            "sort_field": sort_data.get("field"),
            "sort_order": sort_data.get("order"),
        }

        response = requests.get(
            "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search",
            params={k: v for k, v in params.items() if v is not None},
        )

        if response.status_code != 200:
            return html.Div("Error fetching data."), [], [], [], 1, ""

        data = response.json()
        total = data.get("total", 0)
        start = params["start"]

        filter_options = [
            [
                {"label": f"{b['key']} ({b['doc_count']})", "value": b["key"]}
                for b in data.get("aggregations", {}).get(f.id, {}).get("buckets", [])
            ]
            for f in filter_fields
        ]

        return (
            create_table(data.get("results", []), sort_data),
            *filter_options,
            max(1, (total + size - 1) // size) if total > 0 else 1,
            f"{start + 1} – {min(start + size, total)} of {total}" if total > 0 else "",
        )

    @app.callback(
        [Output(f.id, "value") for f in filter_fields],
        [Input(f"clear-{f.id}", "n_clicks") for f in filter_fields],
    )
    def clear_filters(*_):
        triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
        return [
            None if f"clear-{f.id}" == triggered else dash.no_update
            for f in filter_fields
        ]


def resolver(url_parts):
    return (
        details.layout(url_parts[0])
        if len(url_parts) == 1
        else table.layout(filter_fields)
    )
