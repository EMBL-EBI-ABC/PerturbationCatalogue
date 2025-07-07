import base64
import functools
import json
import os

import brotli
import dash
import requests
from dash import html, Output, Input, callback, MATCH
from dash.dependencies import State
import msgpack

from .elastic_table import ElasticTable, Column


# Common parameters.
api_base_url = os.getenv("PERTURBATION_CATALOGUE_BE")

# Prefetch gene lists for cross-references.
genes_sets = {}
for source_id in ["depmap", "mavedb", "perturb-seq"]:
    response = requests.get(f"{api_base_url}/{source_id}/genes")
    response.raise_for_status()
    genes_sets[source_id] = set(response.json())


# Cross-referencing helpers.


def create_cross_reference_link(gene_name):
    """Creates a styled link that pre-fills the search for a gene across all tables."""
    # Create a state where the search is pre-filled for all tables.
    depmap_state = depmap_table.default_state.copy()
    depmap_state["search"] = gene_name
    depmap_query = serialise_state(depmap_state, depmap_table.default_state)

    mavedb_state = mavedb_table.default_state.copy()
    mavedb_state["search"] = gene_name
    mavedb_query = serialise_state(mavedb_state, mavedb_table.default_state)

    perturb_seq_state = perturb_seq_table.default_state.copy()
    perturb_seq_state["search"] = gene_name
    perturb_seq_query = serialise_state(
        perturb_seq_state, perturb_seq_table.default_state
    )

    href = f"/data-portal?depmap={depmap_query}&mavedb={mavedb_query}&perturb_seq={perturb_seq_query}"

    return html.A(
        html.B(gene_name),
        href=href,
        style={"textDecoration": "none", "color": "#198754"},  # Bootstrap success green
        target="_blank",
    )


def format_gene_with_cross_reference(gene_name, current_source_id):
    """Formats a gene name as a cross-reference link if it exists in other data sources.
    Otherwise, returns the plain gene name."""
    other_source_ids = set(genes_sets.keys()) - {current_source_id}
    is_cross_referenced = any(
        gene_name in genes_sets[other_id] for other_id in other_source_ids
    )

    if is_cross_referenced:
        return create_cross_reference_link(gene_name)
    else:
        return gene_name


# DepMap.


def high_dependency_genes(data, display_links=True, max_other_genes=None):
    """Dynamic layout for the list of high dependency genes with toggle for other genes."""

    cross_ref_genes = []
    other_genes_list = []
    for g in data:
        # For DepMap, a cross-reference is to MaveDB or Perturb-Seq.
        if any(
            g["name"] in genes_sets[other_id] for other_id in ["mavedb", "perturb-seq"]
        ):
            cross_ref_genes.append(g)
        else:
            other_genes_list.append(g)

    # Create elements for genes with cross-references.
    cross_ref_elements = [html.I("Genes in other data sources: ")]
    for g in sorted(cross_ref_genes, key=lambda x: x["name"]):
        if display_links:
            link = create_cross_reference_link(g["name"])
        else:
            link = html.B(g["name"])
        cross_ref_elements.append(html.Span([link, html.Span(" ")]))

    # Split other genes into displayed and toggled.
    other_genes = sorted(other_genes_list, key=lambda x: x["name"])
    other_genes_always_displayed = (
        other_genes[:max_other_genes] if max_other_genes else other_genes
    )
    other_genes_toggled = other_genes[max_other_genes:] if max_other_genes else []

    # Create layout.
    return html.Div(
        [
            html.P(cross_ref_elements, style={"marginBottom": "5px"}),
            html.Div(
                [
                    html.I("Other genes: "),
                    html.Span(
                        " ".join([g["name"] for g in other_genes_always_displayed])
                        + " ",
                        className="text-muted",
                        style={"display": "inline"},
                    ),
                    html.Span(
                        " ".join([g["name"] for g in other_genes_toggled]) + " ",
                        className="text-muted",
                        style={"display": "none"},
                        id={"type": "other-genes-list", "index": str(id(data))},
                    ),
                    *(
                        [
                            html.Button(
                                "Expand",
                                id={
                                    "type": "toggle-other-genes",
                                    "index": str(id(data)),
                                },
                                className="btn btn-sm",
                                style={
                                    "marginLeft": "5px",
                                    "backgroundColor": "#f8f9fa",
                                    "color": "#6c757d",
                                    "border": "1px solid #dee2e6",
                                    "padding": "2px 8px",
                                    "fontSize": "12px",
                                    "verticalAlign": "top",
                                },
                            )
                        ]
                        if other_genes_toggled
                        else []
                    ),
                ],
                style={"marginBottom": "0px"},
            ),
        ]
    )


depmap_table = ElasticTable(
    id="depmap",
    api_endpoint=f"{api_base_url}/depmap/search",
    columns=[
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
            display_table=functools.partial(high_dependency_genes, max_other_genes=25),
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
                "Click on a highlighted gene to view relevant functional data in other sources.",
                className="mb-0",
            ),
        ]
    ),
    default_page_size=10,
)

dash.register_page(
    "data-portal-details-depmap",
    path_template="/data-portal/depmap/<record_id>",
    layout=depmap_table.details_layout,
)


# MaveDB.

display_gene_mavedb = functools.partial(
    format_gene_with_cross_reference, current_source_id="mavedb"
)

mavedb_table = ElasticTable(
    id="mavedb",
    api_endpoint=f"{api_base_url}/mavedb/search",
    columns=[
        Column(
            field_name="normalisedGeneName",
            display_name="Gene",
            display_table=display_gene_mavedb,
            display_details=display_gene_mavedb,
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


# Perturb-Seq.

display_gene_perturb_seq = functools.partial(
    format_gene_with_cross_reference, current_source_id="perturb-seq"
)

perturb_seq_table = ElasticTable(
    id="perturb_seq",
    api_endpoint=f"{api_base_url}/perturb-seq/search",
    columns=[
        Column(
            field_name="study_id",
            display_name="Study ID",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="perturbation",
            display_name="Perturbation",
            display_details=display_gene_perturb_seq,
            display_table=display_gene_perturb_seq,
            filterable=True,
        ),
        Column(
            field_name="gene",
            display_name="Gene",
            display_details=display_gene_perturb_seq,
            display_table=display_gene_perturb_seq,
            filterable=True,
        ),
        Column(
            field_name="log2fc",
            display_name="Log2 Fold Change",
            display_details="text",
            sortable=True,
        ),
        Column(
            field_name="pvalue",
            display_name="P-value",
            display_details="text",
            sortable=True,
        ),
        Column(
            field_name="padj",
            display_name="Adjusted P-value",
            display_details="text",
            sortable=True,
            default_sort="asc",
        ),
        Column(
            field_name="mean_control",
            display_name="Mean Control",
            display_details="text",
            sortable=True,
        ),
        Column(
            field_name="mean_perturbed",
            display_name="Mean Perturbed",
            display_details="text",
            sortable=True,
        ),
        # Special column: record_id, to serve as the primary key and link construction.
        Column(
            field_name="record_id",
            display_name="Record ID",
            # display_details="title",
            # display_table=lambda record_id: html.A(
            #     record_id,
            #     href=f"/data-portal/perturb-seq/{record_id}",
            #     className="text-decoration-none text-nowrap",
            # ),
            display_table=False,
            display_details=False,
        ),
    ],
    title="Perturb-Seq",
    description="Perturb-Seq studies are curated and pseudobulk differential expression is computed against controls. Only genes with a significant change are displayed: adjusted p-value <= 0.05 and log2fc is either <= -1 or >= 1",
)

dash.register_page(
    "data-portal-details-perturb-seq",
    path_template="/data-portal/perturb-seq/<record_id>",
    layout=perturb_seq_table.details_layout,
)


# State serialisation and deserialisation.


def serialise_state(state, default_state):
    # Only record the values which are different from the defaults.
    state_diff = {}
    for k, v in state.items():
        if v != default_state[k]:
            state_diff[k] = v
    # Make sure the state is treated as initial load on deserialisation.
    state_diff["initial_load"] = True
    packed_data = msgpack.packb(state_diff)
    return base64.urlsafe_b64encode(brotli.compress(packed_data)).decode().rstrip("=")


def deserialise_state(state):
    if state:
        padded_state = state + "=" * (-len(state) % 4)
        decoded_data = base64.urlsafe_b64decode(padded_state)
        return msgpack.unpackb(brotli.decompress(decoded_data), raw=False)
    else:
        return None


# Main data portal page.


def complete_layout(**kwargs):
    return html.Div(
        [
            depmap_table.table_layout(deserialise_state(kwargs.get("depmap"))),
            mavedb_table.table_layout(deserialise_state(kwargs.get("mavedb"))),
            perturb_seq_table.table_layout(
                deserialise_state(kwargs.get("perturb_seq"))
            ),
        ]
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
    perturb_seq_table.register_callbacks(app)

    @callback(
        [
            Output({"type": "other-genes-list", "index": MATCH}, "style"),
            Output({"type": "toggle-other-genes", "index": MATCH}, "children"),
        ],
        [Input({"type": "toggle-other-genes", "index": MATCH}, "n_clicks")],
        [State({"type": "toggle-other-genes", "index": MATCH}, "children")],
    )
    def toggle_genes(n_clicks, current_label):
        """Toggle the visibility of hidden genes and update button label."""
        if n_clicks is None:
            return {"display": "none"}, "Expand"

        if current_label == "Expand":
            return {"display": "inline"}, "Collapse"
        else:
            return {"display": "none"}, "Expand"

    @callback(
        Output("url", "search"),
        Input("elastic-table-depmap-state", "data"),
        Input("elastic-table-mavedb-state", "data"),
        Input("elastic-table-perturb_seq-state", "data"),
        prevent_initial_call=False,
    )
    def update_url_with_state(depmap_data, mavedb_data, perturb_seq_data):
        base_path = "/data-portal"
        query_string = (
            f"?depmap={serialise_state(depmap_data, depmap_table.default_state)}"
            f"&mavedb={serialise_state(mavedb_data, mavedb_table.default_state)}"
            f"&perturb_seq={serialise_state(perturb_seq_data, perturb_seq_table.default_state)}"
        )
        return query_string
