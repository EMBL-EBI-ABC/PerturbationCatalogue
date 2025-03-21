import dash
from dash import html
import dash_bootstrap_components as dbc

dash.register_page(
    __name__,
    path="/about",
    name="About",
    button="Contact Us",
    description="Here you can find more information about how to get help or propose new features.",
    icon="bi-question-circle",
)

layout = dbc.Container(
    dbc.Row(
        dbc.Col(
            dbc.Card(
                [
                    dbc.CardBody(
                        [
                            html.H3("About", className="card-title"),
                            html.P(
                                "Perturbation Catalogue is an integrated and curated "
                                "resource for human gene perturbation (e.g., CRISPR), "
                                "variant analysis (e.g., MAVE), and expression data "
                                "(e.g., Perturb-seq). This catalogue will provide "
                                "programmatic access and a user-friendly portal for "
                                "comprehensive data exploration, visualisation and "
                                "retrieval, facilitating computational advancements in "
                                "functional genetics and drug target discovery. "
                                "Incorporating data from multiple sources, the catalogue "
                                "will enhance statistical power for protein-coding and "
                                "regulatory target prioritisation. Metadata standard "
                                "development and thorough data curation and harmonisation "
                                "will ensure the catalogue's utility for novel algorithm "
                                "development and machine learning applications."
                            ),
                            html.P(
                                [
                                    "If you have any questions and feature requests please "
                                    "contact us using this email ",
                                    html.A(
                                        "perturbation-catalogue-help@ebi.ac.uk",
                                        href="mailto:perturbation-catalogue-help@ebi.ac.uk",
                                        target="_blank",
                                        style={"textDecoration": "none"},
                                    ),
                                ]
                            ),
                            html.P(
                                [
                                    "Or create issues directly in our GitHub repository ",
                                    html.A(
                                        "https://github.com/EMBL-EBI-ABC/"
                                        "PerturbationCatalogue/issues",
                                        href="https://github.com/EMBL-EBI-ABC/"
                                        "PerturbationCatalogue/issues",
                                        target="_blank",
                                        style={"textDecoration": "none"},
                                    ),
                                ]
                            ),
                        ]
                    ),
                    dbc.CardImg(src="/assets/pc_architecture.png", bottom=True),
                ]
            ),
            md=12,
        ),
        style={"marginTop": "2em", "marginBottom": "2em"},
    )
)
