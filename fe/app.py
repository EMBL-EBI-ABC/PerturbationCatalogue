import dash
from dash import html, Input, Output, State
import dash_bootstrap_components as dbc
import os

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

# Import pages to ensure they are registered
from pages import api, about, dataset, targets, datasets, chat, login

cookie_banner.register_callbacks(app)

# Login form: POST credentials to backend, store token in sessionStorage
app.clientside_callback(
    """
    function(n_clicks, email, password) {
        if (!n_clicks) return [window.dash_clientside.no_update, ""];
        if (!email || !password) return [window.dash_clientside.no_update, "Please enter email and password"];
        var urlEl = document.getElementById("login-backend-url");
        var baseUrl = urlEl ? urlEl.getAttribute("data-url") : "";
        return fetch(baseUrl + "/v1/auth/login", {
            method: "POST",
            headers: {"Content-Type": "application/json"},
            body: JSON.stringify({email: email, password: password})
        })
        .then(function(r) {
            if (r.ok) return r.json();
            return r.json().then(function(err) { return Promise.reject(err); });
        })
        .then(function(data) {
            sessionStorage.setItem("auth_token", data.token);
            sessionStorage.setItem("auth_user", JSON.stringify(data.user));
            window.location.href = "/perturbation-catalogue/chat";
            return [window.dash_clientside.no_update, ""];
        })
        .catch(function(err) {
            var msg = (err && err.detail) ? err.detail : "Invalid email or password";
            return [window.dash_clientside.no_update, msg];
        });
    }
    """,
    [Output("login-redirect", "data"), Output("login-error", "children")],
    Input("login-btn", "n_clicks"),
    [State("login-email", "value"), State("login-password", "value")],
    prevent_initial_call=True,
)

# Logout: clear sessionStorage and redirect to login
app.clientside_callback(
    """
    function(n_clicks) {
        if (!n_clicks) return window.dash_clientside.no_update;
        sessionStorage.removeItem("auth_token");
        sessionStorage.removeItem("auth_user");
        window.location.href = "/perturbation-catalogue/login";
        return window.dash_clientside.no_update;
    }
    """,
    Output("logout-link", "style"),
    Input("logout-link", "n_clicks"),
    prevent_initial_call=True,
)

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
                                src="/perturbation-catalogue/assets/2026-Final-Perturbation-Catalogue-Logo.png",
                                alt="Perturbation Catalogue logo",
                                className="header-logo-img",
                            )
                        ],
                        href="/perturbation-catalogue/",
                    ),
                    html.Nav(
                        [
                            html.A(
                                "Targets",
                                href="/perturbation-catalogue/targets",
                                className="header-link",
                            ),
                            html.A(
                                "Datasets",
                                href="/perturbation-catalogue/datasets",
                                className="header-link",
                            ),
                            html.A(
                                "AI Explorer (Beta)",
                                href="/perturbation-catalogue/chat",
                                className="header-link header-link-ai",
                            ),
                            html.A(
                                "API Documentation",
                                href=os.getenv(
                                    "PERTURBATION_CATALOGUE_BE",
                                    "https://perturbation-catalogue-be-november-prototype-959149465821.europe-west2.run.app",
                                )
                                + "/docs",
                                target="_blank",
                                className="header-link",
                            ),
                            html.A(
                                "About",
                                href="/perturbation-catalogue/about",
                                className="header-link",
                            ),
                            html.A(
                                "Request dataset",
                                href="https://pollunit.com/polls/zl6y1cje-smx6jikfvfwaw",
                                target="_blank",
                                className="header-link",
                            ),
                            html.A(
                                "Logout",
                                id="logout-link",
                                href="#",
                                className="header-link",
                                style={"display": "none"},
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
