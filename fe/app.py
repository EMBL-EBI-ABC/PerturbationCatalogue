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
)

# Initialise callbacks for external components. This ensures that interactivity defined
# in the sub-modules can work across the app. Note that these have to be imported
# *after* app initialisation.
from pages import data_portal
import navbar

cookie_banner.register_callbacks(app)
navbar.register_callbacks(app)
data_portal.register_callbacks(app)

# Inject Google Analytics scripts.
app.index_string = google_analytics.inject

# Footer.
footer = html.Footer(
    dbc.Container(
        [
            html.Div("Powered by ", style={"display": "inline", "margin-right": "5px"}),
            html.Img(
                src="/assets/embl-ebi-logo.png",
                height="30px",
                style={"display": "inline"},
            ),
        ],
        className="text-left py-3 text-light",
        fluid=True,
        style={"padding-left": "20px"},
    ),
    style={"backgroundColor": "rgb(74, 78, 80)", "position": "relative"},
)

# Overall app layout.
app.layout = html.Div(
    [
        # Cookie banner.
        cookie_banner.store,
        cookie_banner.layout,
        # Navigation bar.
        navbar.layout,
        # Main content.
        dash.page_container,
        # Footer.
        footer,
    ],
    className="d-flex flex-column vh-100",
)

# Expose the server variable for Gunicorn.
server = app.server


if __name__ == "__main__":
    app.run_server(debug=True)
