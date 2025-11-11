import dash
from dash import html
import dash_bootstrap_components as dbc

import cookie_banner
import google_analytics

# Initialise the app.
app = dash.Dash(
    __name__,
    external_stylesheets=[
        dbc.themes.BOOTSTRAP,
        "https://cdn.jsdelivr.net/npm/bootstrap-icons/font/bootstrap-icons.css",
    ],
    suppress_callback_exceptions=True,
    use_pages=True,
    # url_base_pathname="/perturbation-catalogue/",
)

# Initialise callbacks for external components. This ensures that interactivity defined
# in the sub-modules can work across the app. Note that these have to be imported
# *after* app initialisation.
# from pages import data_portal
# import navbar

# cookie_banner.register_callbacks(app)
# navbar.register_callbacks(app)
# data_portal.register_callbacks(app)

# Inject Google Analytics scripts.
app.index_string = google_analytics.inject

app.title = "Perturbation Catalogue"

# Footer.
footer = html.Footer(
    html.Div(
        [
            html.Div(
                "Perturbation Catalogue is funded by ",
                style={"display": "inline", "margin-right": "3px", "color": "black"},
            ),
            html.A(
                href="https://www.ebi.ac.uk/",
                target="_blank",
                children=[
                    html.Img(
                        src="/perturbation-catalogue/assets/embl-ebi-logo.png",
                        height="30px",
                        style={
                            "display": "inline",
                            "vertical-align": "top",
                            "margin-right": "10px",
                            "margin-top": "-2px",
                        },
                    ),
                ],
            ),
            html.Div(
                "and ",
                style={"display": "inline", "margin-right": "5px", "color": "black"},
            ),
            html.A(
                href="https://www.opentargets.org/",
                target="_blank",
                children=[
                    html.Img(
                        src="/perturbation-catalogue/assets/open-targets-logo.png",
                        height="32px",
                        style={
                            "display": "inline",
                            "vertical-align": "top",
                            "margin-top": "-3px",
                        },
                    ),
                ],
            ),
        ],
        className="text-left py-3",
        style={"padding-left": "20px", "color": "black"},
    ),
    style={
        "backgroundColor": "#f0f0f0",
        "position": "relative",
    },
)

# Overall app layout.
app.layout = html.Div(
    [
        # Cookie banner.
        cookie_banner.store,
        cookie_banner.layout,
        html.Header(
            dbc.Container(
                [
                    html.A(
                        [
                            html.Img(
                                src="/assets/logo.png",
                                alt="Perturbation Catalogue logo",
                                className="header-logo-img",
                            )
                        ],
                        href="/",
                    ),
                    html.Nav(
                        [
                            html.A(
                                "API Documentation", href="#", className="header-link"
                            ),
                            html.A("About", href="#", className="header-link"),
                        ],
                        className="header-links",
                    ),
                ],
                fluid=True,
                className="header-content",
            ),
            className="app-header",
        ),
        html.Main(dash.page_container, className="app-main"),
        html.Footer(
            dbc.Container(
                html.Div(
                    [
                        html.Span(
                            "Perturbation Catalogue is funded by",
                            className="footer-text",
                        ),
                        html.Div(
                            [
                                html.Img(
                                    src="/assets/embl-ebi-logo.png", alt="EMBL-EBI logo"
                                ),
                                html.Span("and"),
                                html.Img(
                                    src="/assets/open-targets-logo.png",
                                    alt="Open Targets logo",
                                ),
                            ],
                            className="footer-logos",
                        ),
                    ],
                    className="footer-content",
                ),
                fluid=True,
            ),
            className="app-footer",
        ),
    ],
    className="app-shell",
)

# Expose the server variable for Gunicorn.
server = app.server


if __name__ == "__main__":
    app.run(debug=True)
