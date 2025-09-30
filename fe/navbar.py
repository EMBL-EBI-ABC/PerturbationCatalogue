import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback, html

from pages._order import get_pages, get_home_page, REQUEST_DATASET_PAGE


def nav_item(page):
    return dbc.NavItem(
        dbc.NavLink(
            [
                html.I(className=f"bi {page['icon']} me-2"),
                page["name"],
            ],
            href=page["relative_path"],
            id=f"nav-{page['relative_path'][1:]}",
            className="nav-link-custom text-white",
        )
    )


layout = dbc.NavbarSimple(
    children=[
        nav_item(page) for page in get_pages(include_home=False, require_icon=True)
    ]
    + [
        dbc.NavItem(
            dbc.NavLink(
                [
                    html.I(className=f"bi {REQUEST_DATASET_PAGE['icon']} me-2"),
                    REQUEST_DATASET_PAGE["name"],
                ],
                href=REQUEST_DATASET_PAGE["href"],
                className="nav-link-custom text-white",
                external_link=True,
                target="_blank",
            )
        )
    ],
    brand=nav_item(get_home_page()),
    brand_href=get_home_page()["relative_path"],
    color=None,
    className="sticky-top",
    style={"backgroundColor": "var(--custom-color)"},
    fluid=True,
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
            dash.Output(f"nav-{page['relative_path'][1:]}", "style")
            for page in get_pages(require_icon=True)
        ],
        [dash.Input("_pages_location", "pathname")],
    )
    def highlight_bold(pathname):
        """Highlight the active navbar element in bold."""

        # Available navbar element styles.
        default_style = {"color": "white"}
        active_style = {"color": "white", "fontWeight": "bold"}

        # Determine which page is active.
        styles = [
            (
                active_style
                if (pathname == page["relative_path"] == "/")
                or (
                    page["relative_path"] != "/"
                    and pathname.startswith(page["relative_path"])
                )
                else default_style
            )
            for page in get_pages(require_icon=True)
        ]

        # Return the corresponding navigation item styles.
        return styles
