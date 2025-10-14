import dash
from dash import html, dcc
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
from pages import data_portal, perturbations
import navbar

cookie_banner.register_callbacks(app)
navbar.register_callbacks(app)
data_portal.register_callbacks(app)
perturbations.register_callbacks(app)

# Inject Google Analytics scripts.
app.index_string = google_analytics.inject

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
        # URL tracker.
        dcc.Location(id="url", refresh=False),
        # Cookie banner.
        cookie_banner.store,
        cookie_banner.layout,
        # Navigation bar.
        navbar.layout,
        # Main content.
        dbc.Container(
            dash.page_container,
            fluid=True,
            className="p-0 m-0",
        ),
        # Footer.
        footer,
    ],
    className="d-flex flex-column vh-100",
)

# Expose the server variable for Gunicorn.
server = app.server


if __name__ == "__main__":
    app.run(debug=True)
