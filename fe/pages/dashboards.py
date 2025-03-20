import os
import pandas as pd
import dash
import dash_bootstrap_components as dbc
import plotly.express as px
from dash import html, dcc, callback, Output, Input

dash.register_page(
    __name__,
    path="/dashboards",
    name="Dashboards",
    button="Explore Dashboards",
    description="The Dashboards tab provides a visual overview of the existing data.",
    icon="bi-bar-chart-line",
)

DEPMMAP_DATA = pd.read_parquet(
    os.path.abspath("./") + "/parquet_data_sources" + "/depmap_metadata.parquet")
ONCTOTREE_LINEAGE_VALUES = DEPMMAP_DATA["OncotreeLineage"].unique()

MAVEDB_DATA = pd.read_parquet(
    os.path.abspath("./") + "/parquet_data_sources" + "/mavedb_metadata.parquet")
GENE_CATEGORY_VALUES = MAVEDB_DATA["geneCategory"].unique()
SEQUENCE_TYPE_VALUES = MAVEDB_DATA["sequenceType"].unique()

layout = dbc.Container(
    [
        dbc.Row(
            dbc.Col(html.H1("Patients information"), md=4),
            style={"marginTop": "2em"},
        ),
        dbc.Row([
            dbc.Col([
                html.Label("Choose Y-axis category"),
                dcc.Dropdown(["Sex", "Patient Race", "Primary Or Metastasis"], "Sex",
                             id="y-axis-dropdown")
            ], md=3),
            dbc.Col([
                html.Label("Choose X-axis category"),
                dcc.Dropdown(["Age", "Age Category"], "Age",
                             id="x-axis-dropdown")
            ], md=3)
        ]
        ),
        dbc.Row(
            dbc.Col(dbc.Spinner(dcc.Graph(id="age-barchart")), md=12)
        ),
        dbc.Row(
            dbc.Col(html.H1("Samples information"), md=4),
        ),
        dbc.Row(
            dbc.Col([
                html.Label("Choose Oncotree Lineage"),
                dcc.Dropdown(ONCTOTREE_LINEAGE_VALUES, ONCTOTREE_LINEAGE_VALUES[0],
                             id="oncotree-lineage-dropdown")
            ], md=3)
        ),
        dbc.Row(
            [
                dbc.Col(dbc.Spinner(dcc.Graph(id="primary-disease-piechart")), md=6),
                dbc.Col(dbc.Spinner(dcc.Graph(id="collection-site-barchart")), md=6)
            ]
        ),
        dbc.Row(
            dbc.Col(html.H1("Genes/Variants information"), md=8),
        ),
        dbc.Row([
            dbc.Col([
                html.Label("Choose Gene Category"),
                dcc.Dropdown(GENE_CATEGORY_VALUES, GENE_CATEGORY_VALUES[0],
                             id="gene-category-dropdown")
            ], md=3),
            dbc.Col([
                html.Label("Choose Sequence Type"),
                dcc.Dropdown(SEQUENCE_TYPE_VALUES, SEQUENCE_TYPE_VALUES[0],
                             id="sequence-type-dropdown")
            ], md=3)]
        ),
        dbc.Row(
            dbc.Col(dbc.Spinner(dcc.Graph(id="num-variants-histogram")), md=12)
        ),
    ]
)


@callback(
    Output("age-barchart", "figure"),
    Input("y-axis-dropdown", "value"),
    Input("x-axis-dropdown", "value"),
)
def build_patients_dashboards(y_axis_value, x_axis_value):
    formatted_y_axis_value = y_axis_value.replace(" ", "")
    formatted_x_axis_value = x_axis_value.replace(" ", "")
    return px.bar(
        DEPMMAP_DATA.groupby(by=[formatted_y_axis_value, formatted_x_axis_value],
                             as_index=False).size(),
        x=formatted_x_axis_value,
        y="size", color=formatted_y_axis_value,
        labels={formatted_x_axis_value: x_axis_value,
                formatted_y_axis_value: y_axis_value})


@callback(
    Output("primary-disease-piechart", "figure"),
    Output("collection-site-barchart", "figure"),
    Input("oncotree-lineage-dropdown", "value")
)
def build_samples_dashboards(oncotree_lineage_value):
    primary_disease_piechart = px.pie(
        DEPMMAP_DATA[DEPMMAP_DATA["OncotreeLineage"] == oncotree_lineage_value].groupby(
            by="OncotreePrimaryDisease", as_index=False).size(), values="size",
        names="OncotreePrimaryDisease", title="Oncotree Primary Disease",
        hole=.3)
    collection_site_barchart = px.bar(
        DEPMMAP_DATA[DEPMMAP_DATA["OncotreeLineage"] == oncotree_lineage_value].groupby(
            by="SampleCollectionSite",
            as_index=False).size(),
        x="size", y="SampleCollectionSite", title="Samples Collection Site",
        orientation="h", labels={"SampleCollectionSite": "Samples Collectoin Site"})
    return primary_disease_piechart, collection_site_barchart


@callback(
    Output("num-variants-histogram", "figure"),
    Input("gene-category-dropdown", "value"),
    Input("sequence-type-dropdown", "value")
)
def build_num_variants_histogram(gene_category_value, sequence_type_value):
    return px.histogram(
        MAVEDB_DATA[(MAVEDB_DATA["geneCategory"] == gene_category_value) & (
                MAVEDB_DATA["sequenceType"] == sequence_type_value)],
        x="numVariants",
        labels={"numVariants": "Number of Variants"}, log_y=True,
        title="Log(count) Number of Variants", )
