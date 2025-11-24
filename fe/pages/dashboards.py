# import pandas as pd
# import dash
# import dash_bootstrap_components as dbc
# import plotly.express as px
# from dash import html, dcc, callback, Output, Input
# from pathlib import Path

# dash.register_page(
#     __name__,
#     path="/dashboards",
#     relative_path="/dashboards",
#     name="Dashboards",
#     button="Explore Dashboards",
#     description="The Dashboards tab provides a visual overview of the existing data.",
#     icon="bi-bar-chart-line",
# )

# DEPMAP_DATA = pd.read_parquet(
#     Path(".") / "parquet_data_sources" / "depmap_metadata.parquet"
# )
# ONCTOTREE_LINEAGE_VALUES = DEPMAP_DATA["OncotreeLineage"].unique()

# MAVEDB_DATA = pd.read_parquet(
#     Path(".") / "parquet_data_sources" / "mavedb_metadata.parquet"
# )
# GENE_CATEGORY_VALUES = MAVEDB_DATA["geneCategory"].unique()
# SEQUENCE_TYPE_VALUES = [
#     value.upper() if value == "dna" else value.capitalize()
#     for value in MAVEDB_DATA["sequenceType"].unique()
# ]

# PERTURB_SEQ_DATA = pd.read_parquet(
#     Path(".") / "parquet_data_sources" / "perturb-seq.parquet"
# )
# STUDIES_VALUES = list(PERTURB_SEQ_DATA["study_id"].unique())
# STUDIES_VALUES.append("All")

# layout = dbc.Container(
#     [
#         dbc.Row(
#             dbc.Col(html.H1("Individuals information"), md=12),
#             style={"marginTop": "2em"},
#         ),
#         dbc.Row(
#             [
#                 dbc.Col(
#                     [
#                         html.Label("Choose Y-axis category"),
#                         dcc.Dropdown(
#                             ["Sex", "Genetic Ancestry Group", "Primary Or Metastasis"],
#                             "Sex",
#                             id="y-axis-dropdown",
#                         ),
#                     ],
#                     md=3,
#                 ),
#                 dbc.Col(
#                     [
#                         html.Label("Choose X-axis category"),
#                         dcc.Dropdown(
#                             ["Age", "Age Category"], "Age", id="x-axis-dropdown"
#                         ),
#                     ],
#                     md=3,
#                 ),
#             ]
#         ),
#         dbc.Row(dbc.Col(dbc.Spinner(dcc.Graph(id="age-barchart")), md=12)),
#         dbc.Row(
#             dbc.Col(html.H1("Samples information"), md=12),
#         ),
#         dbc.Row(
#             dbc.Col(
#                 [
#                     html.Label("Choose Oncotree Lineage"),
#                     dcc.Dropdown(
#                         ONCTOTREE_LINEAGE_VALUES,
#                         ONCTOTREE_LINEAGE_VALUES[0],
#                         id="oncotree-lineage-dropdown",
#                     ),
#                 ],
#                 md=3,
#             )
#         ),
#         dbc.Row(
#             [
#                 dbc.Col(dbc.Spinner(dcc.Graph(id="primary-disease-piechart")), md=6),
#                 dbc.Col(dbc.Spinner(dcc.Graph(id="collection-site-barchart")), md=6),
#             ]
#         ),
#         dbc.Row(
#             dbc.Col(html.H1("Genes/Variants information"), md=12),
#         ),
#         dbc.Row(
#             [
#                 dbc.Col(
#                     [
#                         html.Label("Choose Gene Category"),
#                         dcc.Dropdown(
#                             GENE_CATEGORY_VALUES,
#                             GENE_CATEGORY_VALUES[0],
#                             id="gene-category-dropdown",
#                         ),
#                     ],
#                     md=3,
#                 ),
#                 dbc.Col(
#                     [
#                         html.Label("Choose Sequence Type"),
#                         dcc.Dropdown(
#                             SEQUENCE_TYPE_VALUES,
#                             SEQUENCE_TYPE_VALUES[0],
#                             id="sequence-type-dropdown",
#                         ),
#                     ],
#                     md=3,
#                 ),
#             ]
#         ),
#         dbc.Row(dbc.Col(dbc.Spinner(dcc.Graph(id="num-variants-histogram")), md=12)),
#         dbc.Row(
#             dbc.Col(
#                 [
#                     html.Label("Choose x-axis range"),
#                     dcc.Slider(
#                         0,
#                         700000,
#                         value=50000,
#                         marks={
#                             0: {"label": "0"},
#                             50000: {"label": "50k"},
#                             100000: {"label": "100k"},
#                             200000: {"label": "200k"},
#                             300000: {"label": "300k"},
#                             400000: {"label": "400k"},
#                             500000: {"label": "500k"},
#                             600000: {"label": "600k"},
#                             700000: {"label": "700k"},
#                         },
#                         id="x-axis-range-slider",
#                     ),
#                 ],
#                 md=12,
#             )
#         ),
#         dbc.Row(
#             dbc.Col(html.H1("Perturb-seq information"), md=12),
#         ),
#         dbc.Row(
#             [
#                 dbc.Col(
#                     [
#                         html.Label("Choose Study"),
#                         dcc.Dropdown(
#                             STUDIES_VALUES,
#                             STUDIES_VALUES[-1],
#                             id="study-dropdown",
#                         ),
#                     ],
#                     md=3,
#                 ),
#                 dbc.Col(
#                     [
#                         html.Label("Choose Perturbed Gene"),
#                         dcc.Dropdown(
#                             id="perturbed-gene-dropdown",
#                         ),
#                     ],
#                     md=3,
#                 ),
#             ]
#         ),
#         dbc.Row(dbc.Col(dbc.Spinner(dcc.Graph(id="volcano-plot")), md=12)),
#     ]
# )


# @callback(
#     Output("age-barchart", "figure"),
#     Input("y-axis-dropdown", "value"),
#     Input("x-axis-dropdown", "value"),
# )
# def build_patients_dashboards(y_axis_value, x_axis_value):
#     formatted_y_axis_value = y_axis_value.replace(" ", "")
#     if formatted_y_axis_value == "GeneticAncestryGroup":
#         formatted_y_axis_value = "PatientRace"
#     formatted_x_axis_value = x_axis_value.replace(" ", "")
#     return px.bar(
#         DEPMAP_DATA.groupby(
#             by=[formatted_y_axis_value, formatted_x_axis_value], as_index=False
#         ).size(),
#         x=formatted_x_axis_value,
#         y="size",
#         color=formatted_y_axis_value,
#         labels={
#             formatted_x_axis_value: x_axis_value,
#             formatted_y_axis_value: y_axis_value,
#             "size": "Count",
#         },
#     )


# @callback(
#     Output("primary-disease-piechart", "figure"),
#     Output("collection-site-barchart", "figure"),
#     Input("oncotree-lineage-dropdown", "value"),
# )
# def build_samples_dashboards(oncotree_lineage_value):
#     primary_disease_piechart = px.pie(
#         DEPMAP_DATA[DEPMAP_DATA["OncotreeLineage"] == oncotree_lineage_value]
#         .groupby(by="OncotreePrimaryDisease", as_index=False)
#         .size(),
#         values="size",
#         names="OncotreePrimaryDisease",
#         title="Oncotree Primary Disease",
#         hole=0.3,
#     )
#     collection_site_barchart = px.bar(
#         DEPMAP_DATA[DEPMAP_DATA["OncotreeLineage"] == oncotree_lineage_value]
#         .groupby(by="SampleCollectionSite", as_index=False)
#         .size(),
#         x="size",
#         y="SampleCollectionSite",
#         title="Samples Collection Site",
#         orientation="h",
#         labels={"SampleCollectionSite": "Samples Collectoin Site", "size": "Count"},
#     )
#     return primary_disease_piechart, collection_site_barchart


# @callback(
#     Output("num-variants-histogram", "figure"),
#     Input("gene-category-dropdown", "value"),
#     Input("sequence-type-dropdown", "value"),
#     Input("x-axis-range-slider", "value"),
# )
# def build_num_variants_histogram(
#     gene_category_value, sequence_type_value, x_axis_range_value
# ):
#     sequence_type_value = sequence_type_value.lower()
#     return px.histogram(
#         MAVEDB_DATA[
#             (MAVEDB_DATA["geneCategory"] == gene_category_value)
#             & (MAVEDB_DATA["sequenceType"] == sequence_type_value)
#         ],
#         x="numVariants",
#         nbins=1000,
#         labels={"numVariants": "Number of Variants", "count": "Count"},
#         range_x=[0, x_axis_range_value],
#     )


# @callback(
#     Output("perturbed-gene-dropdown", "options"),
#     Output("perturbed-gene-dropdown", "value"),
#     Output("volcano-plot", "figure"),
#     Input("study-dropdown", "value"),
#     Input("perturbed-gene-dropdown", "value"),
# )
# def build_volcano_plot(study_value, perturbed_gene_value):
#     if study_value == "All":
#         perturbed_genes = sorted(list(PERTURB_SEQ_DATA["perturbation"].unique()))
#         if perturbed_gene_value not in perturbed_genes:
#             perturbed_gene = perturbed_genes[0]
#         else:
#             perturbed_gene = perturbed_gene_value
#         fig = px.scatter(
#             PERTURB_SEQ_DATA[PERTURB_SEQ_DATA["perturbation"] == perturbed_gene],
#             x="log2fc",
#             y="-log_padj",
#             hover_data=["gene", "perturbation", "log2fc", "-log_padj", "study_id"],
#             labels={"log2fc": "Log(FC)", "-log_padj": "-Log(p-value adjusted)"},
#         )
#     else:
#         perturbed_genes = sorted(
#             list(
#                 PERTURB_SEQ_DATA[PERTURB_SEQ_DATA["study_id"] == study_value][
#                     "perturbation"
#                 ].unique()
#             )
#         )
#         if perturbed_gene_value not in perturbed_genes:
#             perturbed_gene = perturbed_genes[0]
#         else:
#             perturbed_gene = perturbed_gene_value
#         fig = px.scatter(
#             PERTURB_SEQ_DATA[
#                 (PERTURB_SEQ_DATA["study_id"] == study_value)
#                 & (PERTURB_SEQ_DATA["perturbation"] == perturbed_gene)
#             ],
#             x="log2fc",
#             y="-log_padj",
#             hover_data=["gene", "perturbation", "log2fc", "-log_padj", "study_id"],
#             labels={"log2fc": "Log(FC)", "-log_padj": "-Log(p-value adjusted)"},
#         )
#     return perturbed_genes, perturbed_gene, fig
