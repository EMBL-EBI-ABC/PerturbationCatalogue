import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback, html


def layout(pages):
    return dbc.Navbar(
        dbc.Container(
            [
                dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
                dbc.Collapse(
                    dbc.Nav(
                        [
                            dbc.NavItem(
                                dbc.NavLink(
                                    page.name,
                                    href=f"/{page.selector}",
                                    id=f"nav-{page.selector}",
                                    className="nav-link-custom text-white",
                                )
                            )
                            for page in pages
                        ],
                        className="ms-auto",  # Push menu items to the right
                        navbar=True,
                    ),
                    id="navbar-collapse",
                    is_open=False,
                    navbar=True,
                ),
            ],
            className="justify-content-between",
        ),
        color=None,  # Remove default Bootstrap color
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
        return not is_open
