import functools
import json
import os

import dash
from dash import html, Output, Input, callback
from dash.dependencies import State

from .elastic_table import ElasticTable, Column


# Common parameters.

api_base_url = os.getenv("PERTURBATION_CATALOGUE_BE")


# DepMap.


def high_dependency_genes(
    data, display_links=True, max_other_genes=None, component_id=None
):
    """
    Dynamic layout for the list of high dependency genes with Dash callbacks for expand/collapse.

    Parameters:
    - data: List of gene dictionaries
    - display_links: Boolean to determine if genes should be clickable
    - max_other_genes: Maximum number of "Other genes" to display initially.
                       If None, all genes are displayed without expand/collapse functionality.
    - component_id: Unique identifier for this instance of the component
                    (required when using multiple instances)
    """
    # Generate a random ID if none provided
    if component_id is None:
        import uuid

        component_id = str(uuid.uuid4())[:8]

    # Create elements for MaveDB genes
    mavedb_genes = sorted(
        (g for g in data if g.get("xref") == "MaveDB"), key=lambda x: x["name"]
    )

    mavedb_elements = [html.I("Genes in MaveDB: ")] + [
        html.Span(
            [
                html.B(
                    html.A(
                        g["name"],
                        href="#",
                        id={"type": "gene-link", "index": g["name"]},
                        style={"textDecoration": "none"},
                    )
                    if display_links
                    else g["name"]
                ),
                html.Span(" "),
            ]
        )
        for g in mavedb_genes
    ]

    # Get sorted other genes
    other_genes = sorted(
        (g for g in data if g.get("xref") != "MaveDB"), key=lambda x: x["name"]
    )

    # Create elements for other genes
    if max_other_genes is None or len(other_genes) <= max_other_genes:
        # Show all genes if max_other_genes is None or there are fewer genes than the limit
        other_elements = [html.I("Other genes: ")] + [
            html.Span(g["name"] + " ") for g in other_genes
        ]

        return html.Div(
            [
                html.P(mavedb_elements, style={"marginBottom": "5px"}),
                html.P(
                    other_elements,
                    style={"marginBottom": "0px"},
                    className="text-muted",
                ),
            ],
            id=f"gene-list-container-{component_id}",
        )
    else:
        # Create limited view with expand/collapse functionality
        visible_genes = [
            html.Span(g["name"] + " ") for g in other_genes[:max_other_genes]
        ]

        # Create full list (initially hidden)
        all_genes = [html.Span(g["name"] + " ") for g in other_genes]

        # Create the expand/collapse elements with pattern-matching IDs
        expand_link = html.A(
            "[...]",
            href="#",
            id={"type": "expand-genes", "id": component_id},
            style={"marginLeft": "5px", "cursor": "pointer"},
        )

        collapse_link = html.A(
            "[Collapse]",
            href="#",
            id={"type": "collapse-genes", "id": component_id},
            style={"marginLeft": "10px", "cursor": "pointer"},
        )

        return html.Div(
            [
                html.P(mavedb_elements, style={"marginBottom": "5px"}),
                html.P(
                    [html.I("Other genes: ")] + visible_genes + [expand_link],
                    id={"type": "limited-genes-view", "id": component_id},
                    style={"marginBottom": "0px", "display": "block"},
                    className="text-muted",
                ),
                html.P(
                    [html.I("Other genes: ")] + all_genes + [collapse_link],
                    id={"type": "all-genes-view", "id": component_id},
                    style={"marginBottom": "0px", "display": "none"},
                    className="text-muted",
                ),
            ],
            id=f"gene-list-container-{component_id}",
        )


depmap_table = ElasticTable(
    id="depmap",
    api_endpoint=f"{api_base_url}/depmap/search",
    columns=[
        # Special field, not displayed anywhere and only used for filtering.
        Column(
            field_name="xref",
            display_name="Cross-references",
            display_table=False,
            display_details=False,
            filterable=True,
        ),
        # Subtitle.
        Column(
            field_name="OncotreePrimaryDisease",
            display_name="Oncotree Primary Disease",
            sortable=True,
            default_sort="asc",
            display_details="subtitle",
        ),
        # Data columns.
        Column(
            field_name="CellLineName",
            display_name="Cell Line Name",
            display_table=False,
            display_details="text",
        ),
        Column(
            field_name="OncotreeLineage",
            display_name="Oncotree Lineage",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="OncotreeSubtype",
            display_name="Oncotree Subtype",
            display_table=False,
            display_details="text",
        ),
        Column(
            field_name="Age",
            display_name="Age",
            sortable=True,
            display_details="text",
        ),
        Column(
            field_name="AgeCategory",
            display_name="Age Category",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="Sex",
            display_name="Sex",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="PrimaryOrMetastasis",
            display_name="Primary or Metastasis",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="SampleCollectionSite",
            display_name="Sample Collection Site",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="CatalogNumber",
            display_name="Catalog Number",
            display_table=False,
            display_details="text",
        ),
        Column(
            field_name="high_dependency_genes",
            display_name="High Dependency Genes",
            display_table=functools.partial(high_dependency_genes, max_other_genes=30),
            display_details=functools.partial(
                high_dependency_genes, display_links=False
            ),
        ),
        # Title; displayed last in the table view.
        Column(
            field_name="ModelID",
            display_name="Model ID",
            display_details="title",
            display_table=lambda model_id: html.A(
                model_id,
                href=f"/data-portal/depmap/{model_id}",
                className="text-decoration-none text-nowrap",
            ),
        ),
    ],
    details_button_name="View on DepMap",
    details_button_link=lambda ModelID: f"https://depmap.org/portal/cell_line/{ModelID}",
    title="DepMap",
    description=(
        [
            html.P(
                "The table below lists genes which have >95% probability of killing a particular cancer cell line if they are knocked out.",
                className="mb-2",
            ),
            html.P(
                "Genes which are broadly essential to all cells have been excluded.",
                className="mb-2",
            ),
            html.P(
                [
                    "Details on how the data was processed are available in a ",
                    html.A(
                        "Jupyter Notebook",
                        href="https://github.com/EMBL-EBI-ABC/PerturbationCatalogue/blob/main/data_exploration/DepMap/DepMap.ipynb",
                    ),
                    ".",
                ],
                className="mb-2",
            ),
            html.P(
                "Click on a highlighted gene to view relevant functional data in MaveDB.",
                className="mb-0",
            ),
        ]
    ),
)

dash.register_page(
    "data-portal-details-depmap",
    path_template="/data-portal/depmap/<record_id>",
    layout=depmap_table.details_layout,
)


# MaveDB.

mavedb_table = ElasticTable(
    id="mavedb",
    api_endpoint=f"{api_base_url}/mavedb/search",
    columns=[
        Column(
            field_name="normalisedGeneName",
            display_name="Gene",
            display_details="text",
        ),
        Column(
            field_name="geneName",
            display_name="Domain",
            display_details="text",
        ),
        Column(
            field_name="geneCategory",
            display_name="Gene Category",
            display_details="text",
            filterable=True,
        ),
        # Special columns: title and subtitle, for details view only.
        Column(
            field_name="title",
            display_table=False,
            display_details="title",
        ),
        Column(
            field_name="shortDescription",
            display_table=False,
            display_details="subtitle",
        ),
        # Data columns.
        Column(
            field_name="sequenceType",
            display_name="Sequence Type",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="publicationUrl",
            display_name="Publication URL",
            display_table=False,
            display_details="link",
        ),
        Column(
            field_name="publicationYear",
            display_name="Publication Year",
            display_details="text",
            filterable=True,
            sortable=True,
            default_sort="desc",
        ),
        Column(
            field_name="numVariants",
            display_name="Number of Variants",
            display_details="text",
            sortable=True,
        ),
        # Special column: URN, to serve as the primary key and link construction.
        Column(
            field_name="urn",
            display_name="MaveDB ID",
            display_details="text",
            display_table=lambda urn: html.A(
                urn,
                href=f"/data-portal/mavedb/{urn}",
                className="text-decoration-none text-nowrap",
            ),
        ),
    ],
    details_button_name="View on MaveDB",
    details_button_link=lambda urn: f"https://www.mavedb.org/score-sets/{urn}/",
    title="MaveDB",
    description="MaveDB catalogs multiplexed assays of variant effects (MAVEs), which measure the functional impact of genetic variants on specific gene domains.",
)

dash.register_page(
    "data-portal-details-mavedb",
    path_template="/data-portal/mavedb/<record_id>",
    layout=mavedb_table.details_layout,
)


# Main data portal page.


def complete_layout():
    return html.Div(
        [datasource.table_layout() for datasource in (depmap_table, mavedb_table)]
    )


dash.register_page(
    __name__,
    path="/data-portal",
    name="Data Portal",
    button="Open Data Portal",
    description="The Data Portal allows users to sort and filter metadata using a set of predefined filters, and it also has free-text search capabilities.",
    icon="bi-table",
    layout=complete_layout,
)


# Callbacks.


def register_callbacks(app):
    mavedb_table.register_callbacks(app)
    depmap_table.register_callbacks(app)

    @callback(
        Output("elastic-table-mavedb-search", "value"),
        Input({"type": "gene-link", "index": dash.dependencies.ALL}, "n_clicks"),
        State({"type": "gene-link", "index": dash.dependencies.ALL}, "id"),
        prevent_initial_call=True,
    )
    def gene_link_clicked(n_clicks, ids):
        if not any(n_clicks):
            return dash.no_update

        # Get the context of the callback to determine which input triggered it
        ctx = dash.callback_context
        if not ctx.triggered:
            return dash.no_update

        # Get the id of the component that triggered the callback
        triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

        # Extract the index from the triggered component's id
        try:
            triggered_dict = json.loads(triggered_id)
            gene_name = triggered_dict["index"]
            return gene_name
        except:
            return dash.no_update

    # Simple clientside callback for scrolling
    app.clientside_callback(
        """
        function(value) {
            if (value) {
                document.getElementById('elastic-table-mavedb-search').scrollIntoView({behavior: 'smooth'});
            }
            return null;
        }
        """,
        Output("elastic-table-mavedb-search", "_", allow_duplicate=True),
        Input("elastic-table-mavedb-search", "value"),
        prevent_initial_call=True,
    )

    @app.callback(
        [
            Output({"type": "limited-genes-view", "id": dash.MATCH}, "style"),
            Output({"type": "all-genes-view", "id": dash.MATCH}, "style"),
        ],
        [Input({"type": "expand-genes", "id": dash.MATCH}, "n_clicks")],
        [
            State({"type": "limited-genes-view", "id": dash.MATCH}, "style"),
            State({"type": "all-genes-view", "id": dash.MATCH}, "style"),
        ],
    )
    def expand_gene_list(n_clicks, limited_style, all_style):
        if n_clicks is None:
            raise dash.exceptions.PreventUpdate

        limited_style = limited_style or {}
        all_style = all_style or {}

        limited_style.update({"display": "none"})
        all_style.update({"display": "block"})

        return limited_style, all_style

    @app.callback(
        [
            Output(
                {"type": "limited-genes-view", "id": dash.MATCH},
                "style",
                allow_duplicate=True,
            ),
            Output(
                {"type": "all-genes-view", "id": dash.MATCH},
                "style",
                allow_duplicate=True,
            ),
        ],
        [Input({"type": "collapse-genes", "id": dash.MATCH}, "n_clicks")],
        [
            State({"type": "limited-genes-view", "id": dash.MATCH}, "style"),
            State({"type": "all-genes-view", "id": dash.MATCH}, "style"),
        ],
        prevent_initial_call=True,
    )
    def collapse_gene_list(n_clicks, limited_style, all_style):
        if n_clicks is None:
            raise dash.exceptions.PreventUpdate

        limited_style = limited_style or {}
        all_style = all_style or {}

        limited_style.update({"display": "block"})
        all_style.update({"display": "none"})

        return limited_style, all_style
