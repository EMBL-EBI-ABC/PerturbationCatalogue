from collections import namedtuple

import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

import cookie_banner
import google_analytics
import data_portal
import home_page

# URLs for iframes.
DASHBOARDS_URL = "https://lookerstudio.google.com/embed/reporting/86ab32f9-151b-4f91-87eb-060a22f2f890/page/QJFZE"
DATA_ANALYTICS_URL = "https://lookerstudio.google.com/embed/reporting/8e98079d-144b-4583-a108-844b2ed3adf7/page/6eWZE"
ABOUT_URL = "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/redoc"

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True,
)

# Initialise callbacks for external components. This ensures that interactivity defined
# in the sub-modules can work across the app.
cookie_banner.register_callbacks(app)
data_portal.register_callbacks(app)

# Inject Google Analytics scripts.
app.index_string = google_analytics.inject

# Define app pages.
Page = namedtuple("Page", ["name", "selector", "button", "description", "resolver"])
pages = [
    Page(
        name="Home",
        selector="home",
        button=None,
        description=None,
        resolver=lambda url: home_page.layout(pages),
    ),
    Page(
        name="Data Portal",
        selector="data-portal",
        button="Open Data Portal",
        description="The Data Portal allows users to sort and filter metadata using a set of predefined filters, and it also has free-text search capabilities.",
        resolver=lambda url: data_portal.resolver(url.split("/")[2:]),
    ),
    Page(
        name="Dashboards",
        selector="dashboards",
        button="Explore Dashboards",
        description="The Dashboards tab provides a visual overview of the existing data.",
        resolver=lambda url: iframe_layout(DASHBOARDS_URL),
    ),
    Page(
        name="Data Analytics",
        selector="data-analytics",
        button="Data Analytics",
        description="The Data Analytics tab allows users to search the data using the Data Warehouse.",
        resolver=lambda url: iframe_layout(DATA_ANALYTICS_URL),
    ),
    Page(
        name="About",
        selector="about",
        button="Contact Us",
        description="Here you can find more information about how to get help or propose new features.",
        resolver=lambda url: iframe_layout(ABOUT_URL),
    ),
]

# Top navigation bar.
navbar = dbc.Navbar(
    dbc.Container(
        [
            dbc.Nav(
                [
                    dbc.NavItem(
                        dbc.NavLink(
                            page.name,
                            href=f"/{page.selector}",
                            id=f"nav-{page.selector}",
                            style={"color": "white"},
                            className="nav-link-custom",
                        )
                    )
                    for page in pages
                ],
            ),
        ],
        className="justify-content-center",
    ),
    color="success",
    className="sticky-top",
)

# Footer.
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
    className="bg-light text-dark",
)


# Layout for iframe-based pages.
def iframe_layout(url):
    return html.Div(
        [
            html.Iframe(
                src=url,
                style={
                    "width": "100%",
                    "height": "calc(100vh - 130px)",
                    "border": "none",
                },
            )
        ],
        className="mt-0",
    )


# Define the layout.
app.layout = html.Div(
    [
        # Cookie banner.
        cookie_banner.store,
        cookie_banner.layout,
        # Navigation bar.
        navbar,
        dcc.Location(id="url", refresh=False),
        # Main content.
        html.Div(
            id="page-content",
            className="flex-grow-1",
        ),
        # Footer.
        footer,
    ],
    className="d-flex flex-column min-vh-100",
)


# Callback to update page content and navbar active state.
@app.callback(
    [
        dash.Output("page-content", "children"),
        *[
            dash.Output(f"nav-{page.selector}", "style") for page in pages
        ],  # Current navbar element styles.
    ],
    [dash.Input("url", "pathname")],
)
def display_page(pathname):
    # Available navbar element styles.
    default_style = {"color": "white"}
    active_style = {"color": "white", "fontWeight": "bold"}

    # Determine which page is active.
    styles = {
        page.selector: (
            active_style
            if pathname == f"/{page.selector}"
            or (page.selector == "home" and pathname == "/")
            else default_style
        )
        for page in pages
    }

    # Determine the layout to return based on the pathname.
    layout = home_page.layout(pages)
    for page in pages:
        if pathname.startswith("/" + page.selector):
            layout = page.resolver(pathname)

    # Return the corresponding layout and styles.
    return [layout] + [styles[page.selector] for page in pages]


server = app.server

if __name__ == "__main__":
    app.run_server(debug=True)
