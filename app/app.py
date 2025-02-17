import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

# URLs for iframes
DASHBOARDS_URL = "https://lookerstudio.google.com/embed/reporting/86ab32f9-151b-4f91-87eb-060a22f2f890/page/QJFZE"
DATA_ANALYTICS_URL = "https://lookerstudio.google.com/embed/reporting/8e98079d-144b-4583-a108-844b2ed3adf7/page/6eWZE"
ABOUT_URL = "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/redoc"

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Top navigation bar
navbar = dbc.Navbar(
    dbc.Container(
        [
            dbc.Nav(
                [
                    dbc.NavItem(
                        dbc.NavLink("Home", href="/", style={"color": "white"})
                    ),
                    dbc.NavItem(
                        dbc.NavLink(
                            "Data Portal", href="/data-portal", style={"color": "white"}
                        )
                    ),
                    dbc.NavItem(
                        dbc.NavLink(
                            "Dashboards", href="/dashboards", style={"color": "white"}
                        )
                    ),
                    dbc.NavItem(
                        dbc.NavLink(
                            "Data Analytics",
                            href="/data-analytics",
                            style={"color": "white"},
                        )
                    ),
                    dbc.NavItem(
                        dbc.NavLink("About", href="/about", style={"color": "white"})
                    ),
                ],
            ),
        ],
        className="d-flex justify-content-center",
    ),
    color="success",
    dark=True,
    className="fixed-top",
)

# Footer
footer = html.Footer(
    dbc.Container(
        [
            html.Div("Powered by ", style={"display": "inline"}),
            html.Img(
                src="/assets/embl-ebi-logo.png",
                height="20px",
                style={"display": "inline"},
            ),
        ],
        className="text-center py-2",
    ),
    className="fixed-bottom bg-dark text-white",
)

# Home page layout
home_layout = html.Div(
    [
        html.Div(
            html.Img(src="/assets/banner.png", style={"width": "100%"}),
            className="mb-4",
        ),
        html.H1("Perturbation Catalogue", className="text-center my-4"),
        dbc.Container(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Card(
                                [
                                    dbc.CardBody(
                                        [
                                            html.H5("Data Portal"),
                                            html.P(
                                                "The Data Portal allows users to sort and filter metadata using a set of predefined filters, and it also has free-text search capabilities."
                                            ),
                                            dbc.Button(
                                                "Open Data Portal",
                                                href="/data-portal",
                                                color="success",
                                            ),
                                        ]
                                    )
                                ]
                            ),
                            md=4,
                        ),
                        dbc.Col(
                            dbc.Card(
                                [
                                    dbc.CardBody(
                                        [
                                            html.H5("Dashboards"),
                                            html.P(
                                                "The Dashboards tab provides a visual overview of the existing data."
                                            ),
                                            dbc.Button(
                                                "Explore Dashboards",
                                                href="/dashboards",
                                                color="success",
                                            ),
                                        ]
                                    )
                                ]
                            ),
                            md=4,
                        ),
                        dbc.Col(
                            dbc.Card(
                                [
                                    dbc.CardBody(
                                        [
                                            html.H5("Data Analytics"),
                                            html.P(
                                                "The Data Analytics tab allows users to search the data using the Data Warehouse."
                                            ),
                                            dbc.Button(
                                                "Data Analytics",
                                                href="/data-analytics",
                                                color="success",
                                            ),
                                        ]
                                    )
                                ]
                            ),
                            md=4,
                        ),
                    ],
                    className="mb-4",
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Card(
                                [
                                    dbc.CardBody(
                                        [
                                            html.H5("About"),
                                            html.P(
                                                "Here you can find more information about how to get help or propose new features."
                                            ),
                                            dbc.Button(
                                                "Contact us",
                                                href="/about",
                                                color="success",
                                            ),
                                        ]
                                    )
                                ]
                            ),
                            md=4,
                        ),
                    ],
                    className="mb-4 justify-content-center",
                ),
            ]
        ),
    ]
)


# Layout for iframe-based pages
def iframe_layout(url):
    return html.Div(
        [
            html.Iframe(
                src=url,
                style={
                    "width": "100%",
                    "height": "calc(100vh - 100px)",
                    "border": "none",
                },
            )
        ],
        className="mt-5",
    )


# Define the layout
app.layout = html.Div(
    [
        navbar,
        dcc.Location(id="url", refresh=False),
        html.Div(id="page-content", className="mt-5"),
        footer,
    ]
)


# Callback for page routing
@app.callback(dash.Output("page-content", "children"), [dash.Input("url", "pathname")])
def display_page(pathname):
    if pathname == "/dashboards":
        return iframe_layout(DASHBOARDS_URL)
    elif pathname == "/data-analytics":
        return iframe_layout(DATA_ANALYTICS_URL)
    elif pathname == "/about":
        return iframe_layout(ABOUT_URL)
    elif pathname == "/data-portal":
        return html.Div(
            [html.H2("Data Portal - Placeholder")]
        )  # To be implemented later
    return home_layout


if __name__ == "__main__":
    app.run_server(debug=True)
