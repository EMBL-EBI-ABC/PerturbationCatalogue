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
                # Collapsed strip (visible only when sidebar is collapsed)
                html.Div(
                    [
                        html.Button(
                            html.I(className="bi bi-chat-dots-fill"),
                            id="chat-strip-expand-btn",
                            className="chat-strip-btn",
                            title="Open chat",
                        ),
                        html.Span(
                            "",
                            id="chat-unread-badge",
                            className="chat-unread-badge",
                        ),
                    ],
                    className="chat-collapsed-strip",
                ),
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
                html.Div(
                    [
                        # Welcome section (shown when chat is empty)
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.I(
                                            className="bi bi-stars",
                                            style={
                                                "fontSize": "1.25rem",
                                                "color": "#007B53",
                                            },
                                        ),
                                        html.Span(
                                            "How can I help you explore?",
                                            style={
                                                "fontWeight": "600",
                                                "fontSize": "0.9rem",
                                                "color": "#1a1d21",
                                            },
                                        ),
                                    ],
                                    className="chat-welcome-header",
                                ),
                                html.P(
                                    "Ask about genes, proteins, and perturbation experiments.",
                                    className="chat-welcome-subtitle",
                                ),
                                html.Div(
                                    [
                                        html.Button(
                                            [
                                                html.I(
                                                    className="bi bi-search me-2"
                                                ),
                                                "TP53 perturbation data",
                                            ],
                                            className="chat-suggestion-card",
                                            **{
                                                "data-query": "What perturbation data is available for TP53 across all modalities?"
                                            },
                                        ),
                                        html.Button(
                                            [
                                                html.I(
                                                    className="bi bi-card-text me-2"
                                                ),
                                                "BRCA2 protein info",
                                            ],
                                            className="chat-suggestion-card",
                                            **{
                                                "data-query": "Look up the BRCA2 protein and show me a gene card with its function, domains, and disease associations"
                                            },
                                        ),
                                        html.Button(
                                            [
                                                html.I(
                                                    className="bi bi-box me-2"
                                                ),
                                                "KRAS 3D structure",
                                            ],
                                            className="chat-suggestion-card",
                                            **{
                                                "data-query": "Show me the predicted 3D protein structure for KRAS"
                                            },
                                        ),
                                        html.Button(
                                            [
                                                html.I(
                                                    className="bi bi-capsule me-2"
                                                ),
                                                "Is EGFR druggable?",
                                            ],
                                            className="chat-suggestion-card",
                                            **{
                                                "data-query": "Use Pharos to check if EGFR is druggable. What is its target development level?"
                                            },
                                        ),
                                    ],
                                    className="chat-suggestion-cards",
                                ),
                            ],
                            id="chat-welcome",
                            className="chat-welcome",
                        ),
                    ],
                    id="chat-messages",
                    className="chat-messages",
                ),
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
                                dcc.Textarea(
                                    id="chat-input",
                                    placeholder="Ask about perturbation data...",
                                    className="chat-input",
                                ),
                                html.Button(
                                    [html.I(className="bi bi-send-fill")],
                                    id="chat-send-btn",
                                    className="chat-send-btn",
                                ),
                            ],
                            className="chat-input-area",
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
                                        html.Div(
                                            [
                                                html.I(className="bi bi-bar-chart-line"),
                                                html.I(
                                                    className="bi bi-table",
                                                    style={
                                                        "position": "absolute",
                                                        "bottom": "12px",
                                                        "right": "14px",
                                                        "fontSize": "1rem",
                                                        "opacity": "0.5",
                                                    },
                                                ),
                                                html.I(
                                                    className="bi bi-diagram-3",
                                                    style={
                                                        "position": "absolute",
                                                        "top": "12px",
                                                        "left": "14px",
                                                        "fontSize": "1rem",
                                                        "opacity": "0.5",
                                                    },
                                                ),
                                            ],
                                            className="empty-state-illustration",
                                        ),
                                        html.P(
                                            "Your dashboard will build here",
                                            className="empty-state-title",
                                        ),
                                        html.P(
                                            "Ask a question in the chat to start exploring. "
                                            "Charts, tables, networks, and protein structures will appear as tiles.",
                                            className="empty-state-description",
                                        ),
                                        html.Div(
                                            [
                                                html.Span(
                                                    "Try:",
                                                    className="empty-state-try-label",
                                                ),
                                                html.Button(
                                                    "Show TP53 data",
                                                    className="empty-state-try-btn",
                                                    **{
                                                        "data-query": "What perturbation data is available for TP53 across all modalities?"
                                                    },
                                                ),
                                                html.Button(
                                                    "BRCA2 protein",
                                                    className="empty-state-try-btn",
                                                    **{
                                                        "data-query": "Look up the BRCA2 protein and show me a gene card"
                                                    },
                                                ),
                                                html.Button(
                                                    "KRAS structure",
                                                    className="empty-state-try-btn",
                                                    **{
                                                        "data-query": "Show me the predicted 3D protein structure for KRAS"
                                                    },
                                                ),
                                            ],
                                            className="empty-state-try-buttons",
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
