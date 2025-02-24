import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback, html, dcc

from pages._order import get_pages, get_home_page

home_page = get_home_page()

layout = dbc.Navbar(
    [
        dcc.Location(id="navbar_url", refresh=False),
        dbc.Container(
            [
                dbc.Nav(
                    [
                        dbc.NavItem(
                            dbc.NavLink(
                                [
                                    html.I(className=f"bi {home_page['icon']} me-2"),
                                    home_page["name"],
                                ],
                                href=home_page["supplied_path"],
                                id=f"nav-{home_page['supplied_path'][1:]}",
                                className="nav-link-custom text-white",
                            )
                        )
                    ],
                    className="me-auto",
                ),
                dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
                dbc.Collapse(
                    dbc.Nav(
                        [
                            dbc.NavItem(
                                dbc.NavLink(
                                    [
                                        html.I(className=f"bi {page['icon']} me-2"),
                                        page["name"],
                                    ],
                                    href=page["supplied_path"],
                                    id=f"nav-{page['supplied_path'][1:]}",
                                    className="nav-link-custom text-white",
                                )
                            )
                            for page in get_pages(include_home=False, require_icon=True)
                        ],
                        className="ms-auto",
                        navbar=True,
                    ),
                    id="navbar-collapse",
                    is_open=False,
                    navbar=True,
                ),
            ],
            fluid=True,
            className="pe-4",
        ),
    ],
    color=None,
    dark=True,
    className="sticky-top",
    style={"backgroundColor": "rgb(23, 140, 67)"},
)


def register_callbacks(app):

    @callback(
        Output("navbar-collapse", "is_open"),
        Input("navbar-toggler", "n_clicks"),
        State("navbar-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_navbar(n_clicks, is_open):
        """Collapse and expand the navbar."""
        return not is_open

    @callback(
        [
            # Current navbar element styles.
            dash.Output(f"nav-{page['supplied_path'][1:]}", "style")
            for page in get_pages(require_icon=True)
        ],
        [dash.Input("navbar_url", "pathname")],
    )
    def highlight_bold(pathname):
        """Highlight the active navbar element in bold."""

        # Available navbar element styles.
        default_style = {"color": "white"}
        active_style = {"color": "white", "fontWeight": "bold"}

        # Determine which page is active.
        styles = {
            page["supplied_path"][1:]: (
                active_style
                if pathname.startswith(f"/{page['supplied_path'][1:]}")
                or (page["supplied_path"][1:] == "home" and pathname == "/")
                else default_style
            )
            for page in get_pages(require_icon=True)
        }

        # Return the corresponding layout and styles.
        return [
            styles[page["supplied_path"][1:]] for page in get_pages(require_icon=True)
        ]
