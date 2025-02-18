import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

from data_portal import data_portal_layout, details_layout, data_portal_callbacks

# URLs for iframes
DASHBOARDS_URL = "https://lookerstudio.google.com/embed/reporting/86ab32f9-151b-4f91-87eb-060a22f2f890/page/QJFZE"
DATA_ANALYTICS_URL = "https://lookerstudio.google.com/embed/reporting/8e98079d-144b-4583-a108-844b2ed3adf7/page/6eWZE"
ABOUT_URL = "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/redoc"

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True,
)

app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <!-- Google tag (gtag.js) -->
        <script async src="https://www.googletagmanager.com/gtag/js?id=G-2TM7RP1SB5"></script>
        <script>
            window.dataLayer = window.dataLayer || [];
            function gtag(){dataLayer.push(arguments);}
            gtag('js', new Date());

            gtag('config', 'G-2TM7RP1SB5');
        </script>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""

# Top navigation bar
navbar = dbc.Navbar(
    dbc.Container(
        [
            dbc.Nav(
                [
                    dbc.NavItem(
                        dbc.NavLink(
                            "Home",
                            href="/",
                            id="nav-home",
                            style={"color": "white"},
                            className="nav-link-custom",  # Add a custom class for styling
                        )
                    ),
                    dbc.NavItem(
                        dbc.NavLink(
                            "Data Portal",
                            href="/data-portal",
                            id="nav-data-portal",
                            style={"color": "white"},
                            className="nav-link-custom",
                        )
                    ),
                    dbc.NavItem(
                        dbc.NavLink(
                            "Dashboards",
                            href="/dashboards",
                            id="nav-dashboards",
                            style={"color": "white"},
                            className="nav-link-custom",
                        )
                    ),
                    dbc.NavItem(
                        dbc.NavLink(
                            "Data Analytics",
                            href="/data-analytics",
                            id="nav-data-analytics",
                            style={"color": "white"},
                            className="nav-link-custom",
                        )
                    ),
                    dbc.NavItem(
                        dbc.NavLink(
                            "About",
                            href="/about",
                            id="nav-about",
                            style={"color": "white"},
                            className="nav-link-custom",
                        )
                    ),
                ],
                className="ml-auto",
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
                height="29px",
                style={"display": "inline"},
            ),
        ],
        className="text-center py-3",
    ),
    className="fixed-bottom bg-light text-dark",
)

# Home page layout
home_layout = html.Div(
    [
        html.Div(
            style={"position": "relative", "text-align": "center"},
            children=[
                html.Img(
                    src="https://acxngcvroo.cloudimg.io/v7/https://www.embl.org/files/wp-content/uploads/roundels.png",
                    style={"width": "100%", "background-color": "rgb(0, 112, 73, 0.8)"},
                ),
                html.H1(
                    "Perturbation Catalogue",
                    style={
                        "font-size": "300%",
                        "position": "absolute",
                        "top": "50%",
                        "left": "50%",
                        "transform": "translate(-50%, -50%)",
                        "color": "white",
                    },
                ),
            ],
        ),
        html.H1("", className="text-center my-4"),
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
        dcc.Store(id="cookie-consent-store", storage_type="local"),
        # Cookie Banner
        html.Div(
            id="cookie-banner",
            children=[
                html.Div(
                    [
                        html.P(
                            "This website uses cookies to improve your experience.",
                            className="mb-0 me-3",
                        ),
                        html.Button(
                            "Accept",
                            id="accept-cookies",
                            className="btn btn-success me-2",
                        ),
                        html.Button(
                            "Reject", id="reject-cookies", className="btn btn-danger"
                        ),
                    ],
                    className="alert alert-dark d-flex justify-content-center align-items-center gap-2 px-3 py-2",
                    style={
                        "position": "fixed",
                        "bottom": "0",
                        "width": "100%",
                        "zIndex": "1050",
                        "borderRadius": "0",
                    },
                )
            ],
            hidden=True,
        ),
        navbar,
        dcc.Location(id="url", refresh=False),
        html.Div(id="page-content", className="mt-5"),
        footer,
    ]
)


@app.callback(
    [
        dash.Output("cookie-banner", "hidden"),
        dash.Output("cookie-consent-store", "data"),
    ],
    [
        dash.Input("accept-cookies", "n_clicks"),
        dash.Input("reject-cookies", "n_clicks"),
    ],
    [dash.State("cookie-consent-store", "data")],
    prevent_initial_call=False,
)
def handle_cookie_consent(accept_clicks, reject_clicks, store_data):
    ctx = dash.callback_context
    if store_data:
        return True, dash.no_update  # Hide banner if choice was already made

    if ctx.triggered:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if button_id in ["accept-cookies", "reject-cookies"]:
            return True, {"cookie_consent": button_id}

    return False, dash.no_update  # Show banner if no choice stored


# Initialise callbacks for the data portal component.
data_portal_callbacks(app)


# Callback to update page content and navbar active state
@app.callback(
    [
        dash.Output("page-content", "children"),
        dash.Output("nav-home", "style"),
        dash.Output("nav-data-portal", "style"),
        dash.Output("nav-dashboards", "style"),
        dash.Output("nav-data-analytics", "style"),
        dash.Output("nav-about", "style"),
    ],
    [dash.Input("url", "pathname")],
)
def display_page(pathname):
    # Define the default style (non-active)
    default_style = {"color": "white"}

    # Define the active style (bold text)
    active_style = {"color": "white", "fontWeight": "bold"}

    # Determine which page is active
    home_style = active_style if pathname == "/" else default_style
    data_portal_style = (
        active_style if pathname.startswith("/data-portal") else default_style
    )
    dashboards_style = active_style if pathname == "/dashboards" else default_style
    data_analytics_style = (
        active_style if pathname == "/data-analytics" else default_style
    )
    about_style = active_style if pathname == "/about" else default_style

    # Return the corresponding layout and styles
    if pathname == "/dashboards":
        return (
            iframe_layout(DASHBOARDS_URL),
            home_style,
            data_portal_style,
            dashboards_style,
            data_analytics_style,
            about_style,
        )
    elif pathname == "/data-analytics":
        return (
            iframe_layout(DATA_ANALYTICS_URL),
            home_style,
            data_portal_style,
            dashboards_style,
            data_analytics_style,
            about_style,
        )
    elif pathname == "/about":
        return (
            iframe_layout(ABOUT_URL),
            home_style,
            data_portal_style,
            dashboards_style,
            data_analytics_style,
            about_style,
        )
    elif pathname == "/data-portal":
        return (
            data_portal_layout,
            home_style,
            data_portal_style,
            dashboards_style,
            data_analytics_style,
            about_style,
        )
    elif pathname and pathname.startswith("/data-portal/urn:"):
        urn = pathname.split("/")[-1]
        return (
            details_layout(urn),
            home_style,
            data_portal_style,
            dashboards_style,
            data_analytics_style,
            about_style,
        )
    return (
        home_layout,
        home_style,
        data_portal_style,
        dashboards_style,
        data_analytics_style,
        about_style,
    )


if __name__ == "__main__":
    app.run_server(debug=True)
