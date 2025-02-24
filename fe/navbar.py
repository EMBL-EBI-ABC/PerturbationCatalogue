import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback, html


def layout(pages):
    home_page = next(page for page in pages if page.selector == "home")
    other_pages = [page for page in pages if page.selector != "home"]

    return dbc.Navbar(
        dbc.Container(
            [
                dbc.Nav(
                    [
                        dbc.NavItem(
                            dbc.NavLink(
                                [
                                    html.I(className=f"bi {home_page.icon} me-2"),
                                    home_page.name,
                                ],
                                href=f"/{home_page.selector}",
                                id=f"nav-{home_page.selector}",
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
                                        html.I(className=f"bi {page.icon} me-2"),
                                        page.name,
                                    ],
                                    href=f"/{page.selector}",
                                    id=f"nav-{page.selector}",
                                    className="nav-link-custom text-white",
                                )
                            )
                            for page in other_pages
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
        return not is_open
