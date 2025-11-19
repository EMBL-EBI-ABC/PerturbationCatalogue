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
    url_base_pathname="/perturbation-catalogue/",
)

# Initialise callbacks for external components. This ensures that interactivity defined
# in the sub-modules can work across the app. Note that these have to be imported
# *after* app initialisation.
# from pages import data_portal
# import navbar

# Import pages to ensure they are registered
from pages import api, about

# cookie_banner.register_callbacks(app)
# navbar.register_callbacks(app)
# data_portal.register_callbacks(app)

# Inject Google Analytics scripts.
app.index_string = google_analytics.inject

app.title = "Perturbation Catalogue"

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
                                src="/perturbation-catalogue/assets/logo.png",
                                alt="Perturbation Catalogue logo",
                                className="header-logo-img",
                            )
                        ],
                        href="/perturbation-catalogue/",
                    ),
                    html.Nav(
                        [
                            html.A(
                                "API Documentation",
                                href="/perturbation-catalogue/api",
                                className="header-link",
                            ),
                            html.A(
                                "About",
                                href="/perturbation-catalogue/about",
                                className="header-link",
                            ),
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
                                html.A(
                                    href="https://www.ebi.ac.uk/",
                                    target="_blank",
                                    children=[
                                        html.Img(
                                            src="/perturbation-catalogue/assets/embl-ebi-logo.png",
                                            alt="EMBL-EBI logo",
                                        )
                                    ],
                                ),
                                html.Span("and"),
                                html.A(
                                    href="https://www.opentargets.org/",
                                    target="_blank",
                                    children=[
                                        html.Img(
                                            src="/perturbation-catalogue/assets/open-targets-logo.png",
                                            alt="Open Targets logo",
                                        )
                                    ],
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
