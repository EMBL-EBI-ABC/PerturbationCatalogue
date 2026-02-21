import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
import os

from utils import BACKEND_URL

dash.register_page(
    __name__,
    path="/chat",
    relative_path="/chat",
    name="AI Explorer (Beta)",
    description="Ask questions about perturbation data using natural language.",
)

layout = html.Div(
    [
        # Hidden Graph to force Dash to load Plotly.js
        dcc.Graph(style={"display": "none"}, id="chat-plotly-loader"),
        # Hidden div to pass backend URL to JavaScript
        html.Div(
            id="chat-backend-url",
            **{"data-url": BACKEND_URL},
            style={"display": "none"},
        ),
        # ── Chat sidebar ──
        html.Div(
            [
                # Sidebar header
                html.Div(
                    [
                        html.Div(
                            [
                                html.I(className="bi bi-chat-dots-fill", style={"color": "#007B53"}),
                                html.Span("AI Explorer", className="chat-sidebar-title"),
                            ],
                            className="chat-sidebar-brand",
                        ),
                        html.Button(
                            html.I(className="bi bi-chevron-left"),
                            id="chat-collapse-btn",
                            className="chat-collapse-btn",
                            title="Collapse chat",
                        ),
                    ],
                    className="chat-sidebar-header",
                ),
                # Messages area
                html.Div(id="chat-messages", className="chat-messages"),
                # Spinner
                html.Div(
                    [
                        html.I(className="bi bi-flower1 spinner-icon"),
                        html.Span("", className="spinner-status"),
                    ],
                    id="chat-spinner",
                    className="chat-spinner",
                    style={"display": "none"},
                ),
                # Input area (pinned to bottom)
                html.Div(
                    [
                        html.Div(
                            [
                                dcc.Input(
                                    id="chat-input",
                                    type="text",
                                    placeholder="Ask about perturbation data...",
                                    className="chat-input",
                                    autoComplete="off",
                                ),
                                html.Button(
                                    [html.I(className="bi bi-send-fill")],
                                    id="chat-send-btn",
                                    className="chat-send-btn",
                                ),
                            ],
                            className="chat-input-area",
                        ),
                        # Suggestion buttons
                        html.Div(
                            [
                                html.Button(
                                    "TP53 perturbation data",
                                    className="chat-suggestion-btn",
                                    **{
                                        "data-query": "What perturbation data is available for TP53 across all modalities?"
                                    },
                                ),
                                html.Button(
                                    "BRCA2 protein info",
                                    className="chat-suggestion-btn",
                                    **{
                                        "data-query": "Look up the BRCA2 protein and show me a gene card with its function, domains, and disease associations"
                                    },
                                ),
                                html.Button(
                                    "KRAS 3D structure",
                                    className="chat-suggestion-btn",
                                    **{
                                        "data-query": "Show me the predicted 3D protein structure for KRAS"
                                    },
                                ),
                                html.Button(
                                    "Is EGFR druggable?",
                                    className="chat-suggestion-btn",
                                    **{
                                        "data-query": "Use Pharos to check if EGFR is druggable. What is its target development level?"
                                    },
                                ),
                                html.Button(
                                    "CRISPR cancer papers",
                                    className="chat-suggestion-btn",
                                    **{
                                        "data-query": "Find recent papers about CRISPR screens in cancer"
                                    },
                                ),
                                html.Button(
                                    "BRCA1 variants",
                                    className="chat-suggestion-btn",
                                    **{
                                        "data-query": "What are the known protein variants and mutagenesis data for BRCA1?"
                                    },
                                ),
                            ],
                            id="chat-suggestions",
                            className="chat-suggestions",
                        ),
                    ],
                    className="chat-sidebar-bottom",
                ),
            ],
            className="chat-sidebar",
            id="chat-sidebar",
        ),
        # ── Dashboard canvas ──
        html.Div(
            [
                # Dashboard header bar
                html.Div(
                    [
                        html.Div(
                            [
                                html.Button(
                                    html.I(className="bi bi-chat-dots-fill"),
                                    id="chat-expand-btn",
                                    className="chat-expand-btn",
                                    title="Open chat",
                                    style={"display": "none"},
                                ),
                                html.H5(
                                    "Dashboard",
                                    className="mb-0",
                                    style={"fontWeight": "600"},
                                ),
                            ],
                            className="d-flex align-items-center",
                            style={"gap": "0.625rem"},
                        ),
                        html.Div(
                            [
                                html.Span(
                                    "",
                                    id="chat-panel-count",
                                    className="dashboard-panel-count",
                                ),
                                html.Button(
                                    [
                                        html.I(className="bi bi-trash3 me-1"),
                                        "Clear all",
                                    ],
                                    id="chat-clear-portal-btn",
                                    className="dashboard-clear-btn",
                                    style={"display": "none"},
                                ),
                            ],
                            className="d-flex align-items-center",
                            style={"gap": "0.75rem"},
                        ),
                    ],
                    className="dashboard-header",
                ),
                # Dashboard content area
                html.Div(
                    [
                        html.Div(
                            id="chat-data-portal", className="data-portal-content"
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.I(
                                            className="bi bi-bar-chart-line",
                                            style={
                                                "fontSize": "2.5rem",
                                                "color": "#c4cad3",
                                            },
                                        ),
                                        html.P(
                                            "Your dashboard will build here",
                                            style={
                                                "fontWeight": "600",
                                                "color": "#6b7280",
                                                "marginBottom": "0.25rem",
                                                "marginTop": "0.75rem",
                                            },
                                        ),
                                        html.P(
                                            "Ask a question in the chat to start exploring. "
                                            "Charts, tables, networks, and protein structures will appear as tiles.",
                                            style={
                                                "color": "#9ca3af",
                                                "fontSize": "0.85rem",
                                                "maxWidth": "380px",
                                            },
                                        ),
                                    ],
                                    className="text-center",
                                ),
                            ],
                            id="chat-portal-placeholder",
                            className="dashboard-empty-state",
                        ),
                    ],
                    className="dashboard-scroll-area",
                ),
            ],
            className="dashboard-canvas",
        ),
    ],
    className="ai-explorer-layout",
    id="ai-explorer-layout",
)
