# =============================================================================
# IMPORTS
# =============================================================================
import os
import re
import numpy as np
import pandas as pd
import dash
import dash_cytoscape as cyto
from dash import dcc, html, dash_table, Input, Output, State, ctx
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx
import gseapy as gp


# =============================================================================
# DATA LOADING
# Resolve the versioned network CSV from depmap_version.txt so the app
# never needs a code change when DepMap releases new data.
# Load the FDR 10% file (superset of FDR 5%).  All callbacks filter this
# DataFrame at query time based on the user's FDR selection, so we never
# need to re-read from disk.
# =============================================================================
_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "required_data")

def _resolve_data_path():
    version_file = os.path.join(_DATA_DIR, "depmap_version.txt")
    if not os.path.isfile(version_file):
        raise FileNotFoundError(
            f"depmap_version.txt not found at {version_file}.\n"
            f"Run the pipeline (step1 → step2 → step3) to generate the data files."
        )
    release = open(version_file).read().strip()
    match = re.search(r"(\d+Q\d+)", release, re.IGNORECASE)
    if not match:
        raise ValueError(f"Could not parse version from depmap_version.txt: {release!r}")
    version = match.group(1)
    return os.path.join(_DATA_DIR, f"depmap_{version}_gls_whole_coessential_network_FDR_10.csv")

DATA_PATH = _resolve_data_path()

df_all = pd.read_csv(DATA_PATH)

# All genes that appear in the FDR 10% network; used to populate both dropdowns.
all_genes = sorted(set(df_all["source"]).union(set(df_all["target"])))


def _load_dataset_stats():
    """Return (n_cell_lines, n_genes_profiled).

    n_cell_lines:     row count from CRISPRGeneEffect CSV.
    n_genes_profiled: line count from the GLS genes.txt (all genes that passed
                      NA filtering and entered the GLS analysis).
    """
    version_raw = open(os.path.join(_DATA_DIR, "depmap_version.txt")).read().strip()
    m = re.search(r"(\d+Q\d+)", version_raw, re.IGNORECASE)
    version = m.group(1) if m else None

    n_cell_lines = None
    if version:
        crispr_path = os.path.join(_DATA_DIR, f"CRISPRGeneEffect_{version}.csv")
        if os.path.isfile(crispr_path):
            with open(crispr_path) as fh:
                n_cell_lines = sum(1 for _ in fh) - 1  # subtract header

    n_genes_profiled = None
    if version:
        genes_path = os.path.join(_DATA_DIR, f"depmap_{version}_genes.txt")
        if os.path.isfile(genes_path):
            with open(genes_path) as fh:
                n_genes_profiled = sum(1 for _ in fh)

    return n_cell_lines, n_genes_profiled


_N_CELL_LINES, _N_GENES_PROFILED = _load_dataset_stats()
_DEPMAP_VERSION = open(os.path.join(_DATA_DIR, "depmap_version.txt")).read().strip()


# =============================================================================
# APP INITIALISATION
# =============================================================================
cyto.load_extra_layouts()

app = dash.Dash(__name__)
app.title = "DepMap Co-Essentiality Explorer"


# =============================================================================
# COLOUR HELPER
# Maps a correlation value in [-1, 1] to a hex colour on a
# blue (#0072B2) → neutral grey (#dcdcdc) → orange (#E69F00) gradient.
# Palette: Wong (2011) colorblind-safe — distinguishable under deuteranopia,
# protanopia, and tritanopia.
# The three-stop gradient is computed in Python because Cytoscape's mapData
# only interpolates between two colours, producing a muddy midpoint.
# =============================================================================
_NEUTRAL  = np.array([220, 220, 220])
_POSITIVE = np.array([  0, 114, 178])   # #0072B2 — blue  (positive co-essentiality)
_NEGATIVE = np.array([230, 159,   0])   # #E69F00 — orange (negative co-essentiality)

def corr_to_color(c):
    c = float(np.clip(c, -1.0, 1.0))
    # +1 → blue, 0 → neutral grey, -1 → orange
    rgb = (_NEUTRAL + c * (_POSITIVE - _NEUTRAL)).astype(int) if c >= 0 \
          else (_NEUTRAL + (-c) * (_NEGATIVE - _NEUTRAL)).astype(int)
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


# Flat (non-gradient) colours for the multi-gene network:
#  - input genes with at least one displayed pair are a bright, near-black
#    grey — the most prominent nodes
#  - input genes with no displayed pair are a light, dull grey
#  - genes added via degree-of-interaction expansion are a light, dull
#    purple — distinct in hue from the light grey above
#  - edges are coloured by the sign of the underlying correlation, using the
#    same Positive/Negative colours as the single-gene view
# Grey/purple shades are distinguished by lightness/hue independent of the
# blue/orange edge colours, so the palette stays colour-blind friendly.
_MULTI_NODE_COLOR         = "#222222"   # near-black — input gene with a pair (bright/prominent)
_MULTI_NODE_NEUTRAL_COLOR = "#D9D9D9"   # light grey — input gene with no pair (dull)
_MULTI_POSITIVE_COLOR     = "#0072B2"   # blue   — positive correlation edges
_MULTI_NEGATIVE_COLOR     = "#E69F00"   # orange — negative correlation edges
_MULTI_EXTENDED_COLOR     = "#D8BFD8"   # light purple — added by degree of interaction (dull)


# =============================================================================
# CYTOSCAPE STYLESHEETS
# Nodes use data(bg_color) — a hex string computed per-node in the callback —
# so colour reflects the correlation value without needing mapData colours.
# White text + dark outline stays readable on any background shade.
# Edge width is mapped from "weight" (-log10 adj-p, clamped to [0, 10]).
# =============================================================================
_NODE_BASE = {
    "label":             "data(id)",
    "font-size":         "17px",
    "font-weight":       "bold",
    "width":             "62px",
    "height":            "62px",
    "background-color":  "data(bg_color)",
    "color":             "#fff",
    "text-valign":       "center",
    "text-halign":       "center",
    "text-outline-color": "rgba(0,0,0,0.9)",
    "text-outline-width": 1,
    "text-wrap":         "wrap",
    "text-max-width":    "58px",
}

SINGLE_STYLESHEET = [
    {"selector": "node",       "style": _NODE_BASE},
    {"selector": "node.query", "style": {
        "background-color": "#D55E00",
        "width":  "74px",
        "height": "74px",
        "font-size": "18px",
        "border-width": 3,
        "border-color": "#333",
    }},
    {"selector": "edge", "style": {
        "width":      "mapData(weight, 0, 10, 1, 6)",
        "line-color": "#aaa",
        "opacity":    0.6,
    }},
]

# Node fill is set per-node via data(bg_color); edge colour reflects the
# sign of the underlying correlation via data(edge_color) — both computed
# in the callback.
MULTI_STYLESHEET = [
    {"selector": "node", "style": _NODE_BASE},
    {"selector": "edge", "style": {
        "width":      "mapData(weight, 0, 10, 1, 8)",
        "line-color": "data(edge_color)",
        "opacity":    0.65,
    }},
]

# Extra rules appended to MULTI_STYLESHEET when a module row is selected in
# the table — nodes/edges tagged data(module) == <cluster id> get an orange
# highlight; everything else keeps its normal correlation colouring.
def _module_highlight_rules(cluster_id):
    return [
        {"selector": f"node[module = {cluster_id}]", "style": {
            "border-width": 4, "border-color": "#F58518", "border-style": "solid",
        }},
        {"selector": f"edge[module = {cluster_id}]", "style": {
            "line-color": "#F58518", "opacity": 1,
            "width": "mapData(weight, 0, 10, 2, 10)",
        }},
    ]


# =============================================================================
# CO-ESSENTIAL MODULE DETECTION + GO:BP ANNOTATION  (Tab 2 "Gene-list network")
# Modules are connected components of the user's induced sub-network — same
# convention as create_coessential_clusters() in the pilot/sanger scripts —
# numbered by size (1 = largest). Detection is cheap (no API calls) and runs
# live as the user edits their gene list; GO:BP annotation calls the Enrichr
# API per module, so it is gated behind an explicit button to avoid firing a
# query on every keystroke.
# =============================================================================
MULTI_MAX_DEGREE              = 2     # max degree-of-interaction slider value — beyond this the
                                       # network grows too large to interpret
MULTI_MAX_EXPANDED_GENES      = 300   # safety cap on total network size after degree-of-interaction expansion
MULTI_MIN_MODULE_GENES        = 2      # smallest connected component counted as a "module"
MULTI_MIN_GENES_FOR_GO        = 4      # only modules with MORE than 3 genes are GO:BP-annotated
                                       # (a 2-3 gene module is too small for a meaningful enrichment test)
MULTI_MAX_MODULES_TO_ANNOTATE = 10    # cap on Enrichr queries per click (politeness + latency)
MULTI_MAX_GO_GENES            = 200   # modules larger than this are reported but not queried
GO_GENE_SETS                  = "GO_Biological_Process_2025"
GO_ORGANISM                   = "human"
GO_ADJ_P_THRESHOLD            = 0.05  # only GO:BP terms at or below this FDR are shown
GO_COVER_THRESHOLD            = 0.50  # WSC redundancy threshold: a term is redundant if
                                      # ≥50% of its genes are already covered by a more
                                      # significant selected term (mirrors WebGestalt default)

GO_TERM_PATTERN = re.compile(r"^(.*)\s\((GO:\d+)\)$")

def _split_go_term(term):
    """Split an Enrichr term string 'Name (GO:0000000)' into (name, go_id)."""
    match = GO_TERM_PATTERN.match(term)
    return match.groups() if match else (term, "")


def _weighted_set_cover(sig_df):
    """Remove redundant GO:BP terms using a greedy weighted set cover.

    Terms are processed from most to least significant (by Adjusted P-value).
    A term is kept if fewer than GO_COVER_THRESHOLD of its genes are already
    explained by previously selected terms. This mirrors the redundancy-reduction
    strategy used in WebGestalt's Weighted Set Cover output.

    Args:
        sig_df: DataFrame of significant GO:BP terms (Enrichr output, pre-filtered
                to adj. p ≤ GO_ADJ_P_THRESHOLD). Must have 'Genes' and
                'Adjusted P-value' columns.

    Returns:
        Non-redundant subset of sig_df, sorted by adjusted p-value.
    """
    ordered = sig_df.sort_values("Adjusted P-value").reset_index(drop=True)
    covered = set()
    kept = []

    for _, row in ordered.iterrows():
        genes = {g.strip() for g in str(row["Genes"]).split(";") if g.strip()}
        if not genes:
            continue
        already_covered = len(genes & covered) / len(genes)
        if already_covered < GO_COVER_THRESHOLD:
            kept.append(row)
            covered |= genes

    return pd.DataFrame(kept) if kept else pd.DataFrame(columns=sig_df.columns)


def _expand_by_degree(seed_genes, df, max_degree):
    """BFS-expand a gene set along co-essential edges, up to max_degree hops.

    seed_genes: the user's input genes (degree 0).
    df:         FDR-filtered network to search for neighbours (full df_all,
                not just the induced sub-network).

    Returns (all_genes, gene_degree) where gene_degree maps gene -> hop
    distance from the seed set (0 for seed genes). Expansion stops early if
    MULTI_MAX_EXPANDED_GENES is reached.
    """
    gene_degree = {g: 0 for g in seed_genes}
    all_genes = set(seed_genes)
    frontier = set(seed_genes)

    for d in range(1, max_degree + 1):
        if not frontier or len(all_genes) >= MULTI_MAX_EXPANDED_GENES:
            break
        mask = df["source"].isin(frontier) | df["target"].isin(frontier)
        neighbours = set(df.loc[mask, "source"]) | set(df.loc[mask, "target"])
        new_genes = neighbours - all_genes
        if not new_genes:
            break
        room = MULTI_MAX_EXPANDED_GENES - len(all_genes)
        if len(new_genes) > room:
            new_genes = set(list(new_genes)[:room])
        for g in new_genes:
            gene_degree[g] = d
        all_genes |= new_genes
        frontier = new_genes

    return all_genes, gene_degree


def _detect_modules(genes, edges_df):
    """Connected components of the induced sub-network, numbered by size.

    Returns (gene_module, modules) where gene_module maps gene -> cluster id
    (1 = largest; genes outside any >= MULTI_MIN_MODULE_GENES component map to
    0) and modules is an ordered list of {"cluster", "cluster_size", "genes"}.
    """
    graph = nx.Graph()
    graph.add_nodes_from(genes)
    graph.add_edges_from(edges_df[["source", "target"]].itertuples(index=False, name=None))

    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    components = [c for c in components if len(c) >= MULTI_MIN_MODULE_GENES]

    gene_module = {g: cid for cid, comp in enumerate(components, start=1) for g in comp}
    modules = [
        {"cluster": cid, "cluster_size": len(comp), "genes": sorted(comp)}
        for cid, comp in enumerate(components, start=1)
    ]
    return gene_module, modules


# =============================================================================
# DESIGN TOKENS  (EBI Perturbation Catalogue palette)
# =============================================================================
_G      = "#2E7D52"   # primary green
_G_DARK = "#1D5C3A"   # hover / pressed
_G_LITE = "#E8F5EE"   # tinted backgrounds
_BORDER = "#DEE2E6"
_BG     = "#F8F9FA"
_TEXT   = "#212121"
_MUTED  = "#6C757D"


def _panel(header_text, children, header_icon=""):
    """Green-header card matching EBI filter-panel style."""
    return html.Div([
        html.Div(
            f"{header_icon}  {header_text}" if header_icon else header_text,
            style={"background": _G, "color": "#fff", "fontWeight": "600",
                   "fontSize": "20px", "padding": "10px 16px",
                   "borderRadius": "6px 6px 0 0"},
        ),
        html.Div(children, style={
            "background": "#fff", "border": f"1px solid {_BORDER}",
            "borderTop": "none", "padding": "16px",
            "borderRadius": "0 0 6px 6px",
        }),
    ], style={"marginBottom": "24px"})


# =============================================================================
# LAYOUT
# =============================================================================
app.layout = html.Div([


    # ── Main content ─────────────────────────────────────────────────────────
    html.Main([

        # Page title
        html.Div([
            html.H1("Cancer Dependency Map (DepMap) Co-essential Genetic Interaction Network", style={
                "fontSize": "32px", "fontWeight": "700",
                "margin": "0 0 12px 0", "color": _TEXT}),
            html.P(
                "Explore gene co-essentiality relationships derived from DepMap CRISPR "
                "screens. Search for a single gene to see its partners, or paste a gene "
                "list to visualise the network among them.",
                style={"color": _MUTED, "fontSize": "20px",
                       "margin": "0 0 16px 0", "lineHeight": "1.7"}),
            html.P([
                "Genes that are co-essential (perturbation of either gene impairs fitness across many cancer cell lines) "
                "are likely to be functionally related and operate in the same pathway or complex.",
                html.Br(),
                html.Br(),
                "Co-essentiality is computed by applying Generalised Least Squares (GLS) "
                "regression to DepMap CRISPR gene effect scores to correct for cell-line "
                "covariance structure, following ",
                html.Em(
                    "A genome-wide atlas of co-essential modules assigns function to "
                    "uncharacterized genes",
                    style={"color": _TEXT},
                ),
                " — Wainberg et al., ",
                html.A(
                    "Nature Genetics 53, 638–649 (2021)",
                    href="https://doi.org/10.1038/s41588-021-00840-z",
                    target="_blank",
                    style={"color": _G, "textDecoration": "none",
                           "borderBottom": f"1px solid {_G}"},
                ),
                ".",
            ], style={
                "color": _MUTED, "fontSize": "18px",
                "margin": "0 0 20px 0",
                "padding": "14px 18px",
                "lineHeight": "1.7",
                "background": _G_LITE,
                "borderLeft": f"3px solid {_G}",
                "borderRadius": "0 4px 4px 0",
            }),

            # Dataset stats chips
            html.Div([
                html.Span(_DEPMAP_VERSION, style={
                    "fontSize": "17px", "color": _MUTED,
                    "border": f"1px solid {_BORDER}", "borderRadius": "20px",
                    "padding": "5px 14px",
                }),
                html.Span([
                    html.Strong(f"{_N_CELL_LINES:,}" if _N_CELL_LINES else "—",
                                style={"color": _G}),
                    "  cancer cell lines",
                ], style={
                    "fontSize": "17px", "color": _TEXT,
                    "background": _G_LITE, "border": f"1px solid {_G}",
                    "borderRadius": "20px", "padding": "5px 14px",
                }),
                *([ html.Span([
                    html.Strong(f"{_N_GENES_PROFILED:,}", style={"color": _G}),
                    "  genes profiled",
                ], style={
                    "fontSize": "17px", "color": _TEXT,
                    "background": _G_LITE, "border": f"1px solid {_G}",
                    "borderRadius": "20px", "padding": "5px 14px",
                }) ] if _N_GENES_PROFILED else []),
            ], style={"display": "flex", "flexWrap": "wrap", "gap": "8px",
                      "marginBottom": "28px"}),
        ]),

        # FDR filter row
        html.Div([
            html.Div([
                html.Span("FDR (Benjamini-Hochberg) threshold", style={
                    "fontWeight": "600", "fontSize": "18px",
                    "color": _TEXT, "marginRight": "10px"}),
                dcc.Dropdown(
                    id="fdr-filter",
                    options=[
                        {"label": "FDR ≤ 5%",  "value": 0.05},
                        {"label": "FDR ≤ 10%", "value": 0.10},
                    ],
                    value=0.05,
                    clearable=False,
                    style={"width": "150px", "fontSize": "18px"},
                ),
            ], style={"display": "flex", "alignItems": "center"}),
            html.Span(id="network-size-text", style={
                "fontSize": "18px", "color": _MUTED,
                "border": f"1px solid {_BORDER}", "borderRadius": "4px",
                "padding": "4px 12px", "background": "#fff"}),
        ], style={
            "display": "flex", "alignItems": "center",
            "justifyContent": "space-between",
            "padding": "12px 16px", "background": _G_LITE,
            "border": f"1px solid {_BORDER}", "borderRadius": "6px",
            "marginBottom": "28px",
        }),

        # Tabs
        dcc.Tabs(
            id="main-tabs",
            value="tab-single",
            colors={"border": _BORDER, "primary": _G, "background": _BG},
            style={"marginBottom": "0"},
            children=[

                # ── Tab 1: Single-gene explorer ──────────────────────────────
                dcc.Tab(
                    label="Single-gene explorer",
                    value="tab-single",
                    style={"fontWeight": "500", "fontSize": "18px",
                           "padding": "10px 20px", "color": _MUTED},
                    selected_style={"fontWeight": "600", "fontSize": "18px",
                                    "padding": "10px 20px", "color": _G,
                                    "borderTop": f"3px solid {_G}"},
                    children=[html.Div([

                        # Search row
                        html.Div([
                            html.Label("Search gene", style={
                                "fontWeight": "600", "fontSize": "18px",
                                "color": _TEXT, "marginBottom": "6px",
                                "display": "block"}),
                            html.Div([
                                dcc.Dropdown(
                                    id="gene-dropdown",
                                    options=[{"label": g, "value": g} for g in all_genes],
                                    placeholder="Type a gene symbol…",
                                    searchable=True,
                                    clearable=True,
                                    style={"flex": "1", "maxWidth": "500px",
                                           "fontSize": "18px"},
                                ),
                                html.Span("Try:", style={
                                    "fontSize": "16px", "color": _MUTED,
                                    "marginLeft": "12px", "whiteSpace": "nowrap",
                                }),
                                html.Button("TP53", id="example-btn-tp53", n_clicks=0,
                                    style={"background": "transparent", "color": _G,
                                           "border": f"1px solid {_G}", "borderRadius": "4px",
                                           "padding": "4px 10px", "fontSize": "16px",
                                           "cursor": "pointer"}),
                                html.Button("BRCA1", id="example-btn-brca1", n_clicks=0,
                                    style={"background": "transparent", "color": _G,
                                           "border": f"1px solid {_G}", "borderRadius": "4px",
                                           "padding": "4px 10px", "fontSize": "16px",
                                           "cursor": "pointer"}),
                                html.Button("KRAS", id="example-btn-kras", n_clicks=0,
                                    style={"background": "transparent", "color": _G,
                                           "border": f"1px solid {_G}", "borderRadius": "4px",
                                           "padding": "4px 10px", "fontSize": "16px",
                                           "cursor": "pointer"}),
                            ], style={"display": "flex", "alignItems": "center",
                                      "gap": "8px", "flexWrap": "wrap"}),
                        ], style={"marginTop": "24px", "marginBottom": "20px"}),

                        # Summary badge + download button
                        html.Div([
                            html.Div(id="summary-text", style={
                                "fontSize": "18px", "color": _MUTED,
                                "fontStyle": "italic", "flex": "1"}),
                            html.Button(
                                "Download all partners (CSV)",
                                id="single-download-btn",
                                n_clicks=0,
                                disabled=True,
                                style={
                                    "background": _G, "color": "#fff",
                                    "border": "none", "borderRadius": "4px",
                                    "padding": "8px 16px", "fontSize": "16px",
                                    "fontWeight": "600", "cursor": "pointer",
                                    "whiteSpace": "nowrap",
                                },
                            ),
                            dcc.Download(id="single-download"),
                        ], style={"display": "flex", "alignItems": "center",
                                  "gap": "16px", "marginBottom": "20px"}),

                        # Chart + table row
                        html.Div([
                            dcc.Loading(
                                dcc.Graph(id="partner-bar",
                                          config={"displayModeBar": False},
                                          style={"flex": "1", "minWidth": "0"}),
                                type="circle", color=_G,
                                style={"flex": "1"},
                            ),
                            html.Div([
                                html.P(
                                    "Co-essential partners — significance & correlation",
                                    style={"fontWeight": "600", "fontSize": "20px",
                                           "color": _TEXT, "margin": "0 0 8px 0"}),
                                html.P(
                                    "Each row is a gene whose CRISPR essentiality profile "
                                    "co-varies significantly with the query gene across cancer "
                                    "cell lines. Adj. p-value is BH-corrected; correlation "
                                    "reflects the direction of co-essentiality "
                                    "(positive = both essential together).",
                                    style={"fontSize": "18px", "color": _MUTED,
                                           "margin": "0 0 12px 0", "lineHeight": "1.5"}),
                                dash_table.DataTable(
                                    id="partner-table",
                                    columns=[
                                        {"name": "PARTNER GENE", "id": "partner"},
                                        {"name": "GLS P-VALUE",      "id": "pvalue",
                                         "type": "numeric",
                                         "format": {"specifier": ".2e"}},
                                        {"name": "GLS ADJ. P-VALUE", "id": "pvalue_adj",
                                         "type": "numeric",
                                         "format": {"specifier": ".2e"}},
                                        {"name": "CORRELATION",  "id": "corr_genes",
                                         "type": "numeric",
                                         "format": {"specifier": ".3f"}},
                                    ],
                                    sort_action="native",
                                    page_size=15,
                                    style_table={"overflowX": "auto", "border": "none"},
                                    style_cell={
                                        "textAlign": "left", "padding": "10px 12px",
                                        "fontSize": "18px", "border": "none",
                                        "borderBottom": f"1px solid {_BORDER}",
                                        "fontFamily": "inherit",
                                    },
                                    style_header={
                                        "fontWeight": "700", "fontSize": "18px",
                                        "color": _MUTED, "letterSpacing": "0.05em",
                                        "border": "none",
                                        "borderBottom": f"2px solid {_BORDER}",
                                        "background": "#fff",
                                    },
                                    style_data_conditional=[{
                                        "if": {"state": "selected"},
                                        "backgroundColor": _G_LITE,
                                        "border": f"1px solid {_G}",
                                    }],
                                ),
                            ], style={"flex": "1", "paddingLeft": "28px", "minWidth": "0"}),
                        ], style={"display": "flex", "alignItems": "flex-start",
                                  "marginBottom": "32px"}),

                        # Network panel
                        _panel("Co-essential genetic interaction network", [
                            html.P(
                                "Partner nodes coloured by "
                                "correlation (edge thickness = −log₁₀(adj. p-value)). ",
                                style={"color": _MUTED, "fontSize": "20px",
                                       "margin": "0 0 10px 0"}),
                            html.Div([
                                html.Span("■ Query gene", style={
                                    "fontSize": "18px", "color": "#D55E00",
                                    "fontWeight": "700", "marginRight": "20px"}),
                                html.Span("■ Positive", style={
                                    "fontSize": "18px", "color": "#0072B2",
                                    "fontWeight": "700", "marginRight": "4px"}),
                                html.Span("(co-essential pair)", style={
                                    "fontSize": "18px", "color": _MUTED,
                                    "marginRight": "20px"}),
                                html.Span("■ Neutral", style={
                                    "fontSize": "18px", "color": "#b0b0b0",
                                    "fontWeight": "700", "marginRight": "24px"}),
                                html.Span("■ Negative", style={
                                    "fontSize": "18px", "color": "#E69F00",
                                    "fontWeight": "700", "marginRight": "4px"}),
                                html.Span("(anti-correlated pair)", style={
                                    "fontSize": "18px", "color": _MUTED}),
                            ], style={"display": "flex", "alignItems": "center",
                                      "flexWrap": "wrap", "marginBottom": "12px",
                                      "padding": "6px 10px", "background": _G_LITE,
                                      "borderRadius": "4px"}),
                            cyto.Cytoscape(
                                id="cyto-graph",
                                layout={"name": "cose"},
                                style={"width": "100%", "height": "500px",
                                       "border": f"1px solid {_BORDER}",
                                       "borderRadius": "4px"},
                                stylesheet=SINGLE_STYLESHEET,
                                elements=[],
                            ),
                        ]),

                    ], style={"padding": "0 4px"})],
                ),

                # ── Tab 2: Gene-list network ──────────────────────────────────
                dcc.Tab(
                    label="Gene-list network",
                    value="tab-multi",
                    style={"fontWeight": "500", "fontSize": "18px",
                           "padding": "10px 20px", "color": _MUTED},
                    selected_style={"fontWeight": "600", "fontSize": "18px",
                                    "padding": "10px 20px", "color": _G,
                                    "borderTop": f"3px solid {_G}"},
                    children=[html.Div([

                        # Gene input panel
                        _panel("Gene list input", [
                            html.P(
                                "Paste or type gene symbols — one per line, or "
                                "comma-, tab-, or space-separated.",
                                style={"color": _MUTED, "fontSize": "18px",
                                       "margin": "0 0 10px 0"}),
                            dcc.Textarea(
                                id="multi-gene-input",
                                placeholder="BRCA1\nTP53, KRAS\nPTEN SMAD4",
                                style={
                                    "width": "100%", "height": "130px",
                                    "fontFamily": "monospace", "fontSize": "18px",
                                    "padding": "10px", "borderRadius": "4px",
                                    "border": f"1px solid {_BORDER}",
                                    "resize": "vertical", "boxSizing": "border-box",
                                    "outline": "none",
                                },
                            ),
                            html.Button(
                                "Load example gene list",
                                id="multi-example-btn",
                                n_clicks=0,
                                style={
                                    "marginTop": "8px",
                                    "background": "transparent", "color": _G,
                                    "border": f"1px solid {_G}", "borderRadius": "4px",
                                    "padding": "5px 12px", "fontSize": "16px",
                                    "cursor": "pointer",
                                },
                            ),
                        ]),

                        # Summary badge
                        html.Div([
                            html.Div(id="multi-summary-text", style={
                                "fontSize": "18px", "color": _MUTED,
                                "fontStyle": "italic", "flex": "1"}),
                            html.Button(
                                "Download pairs (CSV)",
                                id="multi-download-btn",
                                n_clicks=0,
                                disabled=True,
                                style={
                                    "background": _G, "color": "#fff",
                                    "border": "none", "borderRadius": "4px",
                                    "padding": "8px 16px", "fontSize": "16px",
                                    "fontWeight": "600", "cursor": "pointer",
                                    "whiteSpace": "nowrap",
                                },
                            ),
                            dcc.Download(id="multi-download"),
                        ], style={"display": "flex", "alignItems": "center",
                                  "gap": "16px", "marginBottom": "20px"}),

                        # Network (sticky) + co-essential modules (below)
                        html.Div([

                            # Network panel — sticks below the site header while
                            # the user scrolls down to interact with the modules table
                            html.Div(_panel("Co-essential genetic interaction network", [
                                html.P(
                                    "Edges coloured by the sign of the correlation "
                                    "(thickness = −log₁₀(adj. p-value)). "
                                    "Nodes coloured by gene group (see legend).",
                                    style={"color": _MUTED, "fontSize": "20px",
                                           "margin": "0 0 10px 0"}),
                                html.Div([
                                    html.Label("Degree of interaction (N) from input genes", style={
                                        "fontWeight": "600", "fontSize": "18px",
                                        "color": _TEXT, "marginBottom": "6px",
                                        "display": "block"}),
                                    html.P(
                                        "Expand the network to include genes that are "
                                        "co-essential with your input genes within N step(s). "
                                        "Added genes are shown in purple (see legend).",
                                        style={"color": _MUTED, "fontSize": "16px",
                                               "margin": "0 0 10px 0", "lineHeight": "1.5"}),
                                    dcc.Slider(
                                        id="multi-degree-slider",
                                        min=0, max=MULTI_MAX_DEGREE, step=1, value=0,
                                        marks={i: str(i) for i in range(MULTI_MAX_DEGREE + 1)},
                                    ),
                                ], style={"marginBottom": "16px", "maxWidth": "420px"}),
                                html.Div([
                                    html.Div([
                                        html.Span("Nodes:", style={
                                            "fontSize": "18px", "fontWeight": "700",
                                            "color": _TEXT, "marginRight": "10px"}),
                                        html.Span("■ Input gene", style={
                                            "fontSize": "18px", "color": _MULTI_NODE_COLOR,
                                            "fontWeight": "700", "marginRight": "4px"}),
                                        html.Span("(has a pair)", style={
                                            "fontSize": "18px", "color": _MUTED,
                                            "marginRight": "20px"}),
                                        html.Span("■ Input gene", style={
                                            "fontSize": "18px", "color": _MULTI_NODE_NEUTRAL_COLOR,
                                            "fontWeight": "700", "marginRight": "4px"}),
                                        html.Span("(no pair)", style={
                                            "fontSize": "18px", "color": _MUTED,
                                            "marginRight": "20px"}),
                                        html.Span("■ Extended", style={
                                            "fontSize": "18px", "color": _MULTI_EXTENDED_COLOR,
                                            "fontWeight": "700", "marginRight": "4px"}),
                                        html.Span("(added by degree of interaction)", style={
                                            "fontSize": "18px", "color": _MUTED}),
                                    ], style={"display": "flex", "alignItems": "center",
                                              "flexWrap": "wrap", "marginBottom": "6px"}),
                                    html.Div([
                                        html.Span("Edges:", style={
                                            "fontSize": "18px", "fontWeight": "700",
                                            "color": _TEXT, "marginRight": "10px"}),
                                        html.Span("━ Positive", style={
                                            "fontSize": "18px", "color": _MULTI_POSITIVE_COLOR,
                                            "fontWeight": "700", "marginRight": "4px"}),
                                        html.Span("(co-essential pair)", style={
                                            "fontSize": "18px", "color": _MUTED,
                                            "marginRight": "20px"}),
                                        html.Span("━ Negative", style={
                                            "fontSize": "18px", "color": _MULTI_NEGATIVE_COLOR,
                                            "fontWeight": "700", "marginRight": "4px"}),
                                        html.Span("(anti-correlated pair)", style={
                                            "fontSize": "18px", "color": _MUTED}),
                                    ], style={"display": "flex", "alignItems": "center",
                                              "flexWrap": "wrap"}),
                                ], style={"marginBottom": "12px",
                                          "padding": "6px 10px", "background": _G_LITE,
                                          "borderRadius": "4px"}),
                                cyto.Cytoscape(
                                    id="multi-cyto-graph",
                                    layout={"name": "cose"},
                                    style={"width": "100%", "height": "750px",
                                           "border": f"1px solid {_BORDER}",
                                           "borderRadius": "4px"},
                                    stylesheet=MULTI_STYLESHEET,
                                    elements=[],
                                ),
                            ]), style={
                                "position": "sticky",
                                "top": "56px",
                                "zIndex": "10",
                                "background": _BG,
                                "paddingBottom": "16px",
                                "marginBottom": "8px",
                            }),

                            # Co-essential modules + GO:BP annotation panel
                            html.Div(_panel("Co-essential modules", [
                                html.P(
                                    "Genes connected within the network are grouped into "
                                    "co-essential modules (connected components, numbered by "
                                    "size — 1 = largest). Modules with more than 3 genes are "
                                    "annotated with all significant GO Biological Process terms "
                                    "(adj. p ≤ 5%). Redundant terms are filtered using a "
                                    "Weighted Set Cover algorithm (coverage threshold 50%) — only non-redundant representative "
                                    "terms are shown. Smaller modules are listed without "
                                    "annotation. Click a row to highlight that module in the "
                                    "network above.",
                                    style={"color": _MUTED, "fontSize": "18px",
                                           "margin": "0 0 14px 0", "lineHeight": "1.5"}),
                                html.Div([
                                    html.Button(
                                        "Find modules & annotate with GO:BP",
                                        id="annotate-modules-btn",
                                        n_clicks=0,
                                        style={
                                            "background": _G, "color": "#fff", "border": "none",
                                            "borderRadius": "4px", "padding": "10px 18px",
                                            "fontSize": "18px", "fontWeight": "600",
                                            "cursor": "pointer",
                                        },
                                    ),
                                    dcc.Dropdown(
                                        id="module-cluster-filter",
                                        options=[],
                                        value=None,
                                        placeholder="Filter by cluster…",
                                        clearable=True,
                                        disabled=True,
                                        style={"width": "240px", "fontSize": "18px"},
                                    ),
                                    html.Button(
                                        "Download GO:BP terms (CSV)",
                                        id="multi-go-download-btn",
                                        n_clicks=0,
                                        disabled=True,
                                        style={
                                            "background": _G, "color": "#fff",
                                            "border": "none", "borderRadius": "4px",
                                            "padding": "10px 18px", "fontSize": "18px",
                                            "fontWeight": "600", "cursor": "pointer",
                                            "whiteSpace": "nowrap",
                                        },
                                    ),
                                    dcc.Download(id="multi-go-download"),
                                ], style={"display": "flex", "alignItems": "center",
                                          "gap": "16px", "marginBottom": "16px"}),
                                html.Div(id="multi-module-status", style={
                                    "marginBottom": "14px", "fontSize": "18px",
                                    "color": _MUTED, "fontStyle": "italic"}),
                                dcc.Loading(
                                    dash_table.DataTable(
                                        id="multi-module-table",
                                        columns=[
                                            {"name": "CLUSTER",        "id": "cluster"},
                                            {"name": "MODULE SIZE",    "id": "cluster_size"},
                                            {"name": "GO:ID",          "id": "go_id"},
                                            {"name": "GO:BP TERM",     "id": "go_term"},
                                            {"name": "P-VALUE",        "id": "p_value",
                                             "type": "numeric",
                                             "format": {"specifier": ".2e"}},
                                            {"name": "ADJ. P-VALUE",   "id": "p_value_adj",
                                             "type": "numeric",
                                             "format": {"specifier": ".2e"}},
                                        ],
                                        data=[],
                                        sort_action="native",
                                        page_size=10,
                                        style_table={"overflowX": "auto", "border": "none"},
                                        style_cell={
                                            "textAlign": "left", "padding": "10px 12px",
                                            "fontSize": "18px", "border": "none",
                                            "borderBottom": f"1px solid {_BORDER}",
                                            "fontFamily": "inherit",
                                        },
                                        style_header={
                                            "fontWeight": "700", "fontSize": "18px",
                                            "color": _MUTED, "letterSpacing": "0.05em",
                                            "border": "none",
                                            "borderBottom": f"2px solid {_BORDER}",
                                            "background": "#fff",
                                        },
                                        style_data_conditional=[{
                                            "if": {"state": "active"},
                                            "backgroundColor": _G_LITE,
                                            "border": f"1px solid {_G}",
                                        }],
                                    ),
                                    type="circle", color=_G,
                                ),
                            ])),

                        ]),

                        dcc.Store(id="multi-modules-store"),
                        dcc.Store(id="multi-go-rows-store"),
                        dcc.Store(id="multi-gene-list-version", data=0),

                    ], style={"padding": "0 4px"})],
                ),
            ],
        ),

    ], style={"maxWidth": "1200px", "margin": "0 auto", "padding": "32px 24px"}),

], style={
    "fontFamily": "'Segoe UI', Arial, sans-serif",
    "background": _BG, "minHeight": "100vh", "color": _TEXT,
})


# =============================================================================
# CALLBACKS
# =============================================================================

# --- Network size badge (updates whenever FDR changes) ---
@app.callback(
    Output("network-size-text", "children"),
    Input("fdr-filter", "value"),
)
def update_network_size(fdr):
    n = (df_all["pvalue_adj"] <= fdr).sum()
    pct = int(fdr * 100)
    return f"{n:,} significant gene pairs at FDR ≤ {pct}%"


# --- Single-gene explorer ---
@app.callback(
    Output("summary-text",          "children"),
    Output("partner-bar",           "figure"),
    Output("partner-table",         "data"),
    Output("cyto-graph",            "elements"),
    Output("single-download-btn",   "disabled"),
    Input("gene-dropdown",  "value"),
    Input("fdr-filter",     "value"),
)
def update_single(gene, fdr):
    df = df_all[df_all["pvalue_adj"] <= fdr]
    pct = int(fdr * 100)

    if not gene:
        placeholder = go.Figure()
        placeholder.update_layout(
            xaxis={"visible": False}, yaxis={"visible": False},
            paper_bgcolor="white", plot_bgcolor="white", height=500,
            margin={"l": 0, "r": 0, "t": 0, "b": 0},
            annotations=[{"text": "Select a gene above to see its co-essential partners",
                           "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5,
                           "showarrow": False,
                           "font": {"size": 15, "color": _MUTED}}],
        )
        return "Select a gene to explore its co-essential partners.", placeholder, [], [], True

    mask = (df["source"] == gene) | (df["target"] == gene)
    sub = df[mask].copy()
    sub["partner"] = sub.apply(
        lambda r: r["target"] if r["source"] == gene else r["source"], axis=1
    )
    sub = sub[["partner", "pvalue", "pvalue_adj", "corr_genes"]].sort_values("pvalue_adj")

    summary = f"{gene} has {len(sub):,} co-essential partner(s) at FDR ≤ {pct}%."

    top = sub.head(20).copy()
    top["-log10(FDR)"] = -np.log10(top["pvalue_adj"].clip(lower=1e-300))
    fig = px.bar(
        top, x="-log10(FDR)", y="partner", orientation="h",
        title=f"Top co-essential partners of {gene}  (FDR ≤ {pct}%)",
        labels={"partner": ""},
        height=500,
    )
    fig.update_traces(marker_color=_G)
    fig.update_layout(
        yaxis={"categoryorder": "total ascending",
               "tickfont": {"size": 16}, "ticksuffix": "  "},
        plot_bgcolor="white", paper_bgcolor="white",
        font={"family": "'Segoe UI', Arial, sans-serif", "color": _TEXT, "size": 16},
        title={"font": {"size": 18, "color": _TEXT}},
        xaxis={"gridcolor": _BORDER, "linecolor": _BORDER, "title": "−log₁₀(FDR)",
               "tickfont": {"size": 17}, "title_font": {"size": 18}},
        margin={"l": 10, "t": 52, "b": 52, "r": 10},
    )

    top_partners = set(top["partner"].tolist())

    # corr_genes per partner (used for node colour)
    partner_corr = sub.set_index("partner")["corr_genes"].to_dict()

    # Query node: orange via .query selector; bg_color unused but set for consistency
    nodes = [{"data": {"id": gene, "bg_color": "#D55E00"}, "classes": "query"}]
    for p in top_partners:
        nodes.append({"data": {"id": p, "bg_color": corr_to_color(partner_corr.get(p, 0.0))}})

    # Spoke edges: query → each partner (fixed mid-weight so they're visible)
    edges = [{"data": {"source": gene, "target": p, "weight": 3.0}}
             for p in top_partners]

    # Cross-edges among partners, weighted by significance
    partner_mask = df["source"].isin(top_partners) & df["target"].isin(top_partners)
    for _, row in df[partner_mask].iterrows():
        w = min(-np.log10(max(row["pvalue_adj"], 1e-300)), 10)
        edges.append({"data": {"source": row["source"], "target": row["target"],
                                "weight": w}})

    return summary, fig, sub.to_dict("records"), nodes + edges, False


# --- Tab 1: example-gene buttons load a gene into the dropdown ---
@app.callback(
    Output("gene-dropdown", "value"),
    Input("example-btn-tp53",  "n_clicks"),
    Input("example-btn-brca1", "n_clicks"),
    Input("example-btn-kras",  "n_clicks"),
    prevent_initial_call=True,
)
def set_example_gene(_tp53, _brca1, _kras):
    gene_map = {"example-btn-tp53": "TP53", "example-btn-brca1": "BRCA1",
                "example-btn-kras": "KRAS"}
    return gene_map.get(ctx.triggered_id, dash.no_update)


# --- Multi-gene network ---
_ALL_GENES_SET = set(all_genes)

_EXAMPLE_GENES = (
    "BRCA1\nBRCA2\nPALB2\nRAD51\nATM\nCHEK2\nTP53\nPTEN\nRB1\n"
    "FANCL\nFANCA\nFANCD2\nFANCG\nFANCI\nMDM2\nMDM4"
)

@app.callback(
    Output("multi-summary-text",       "children"),
    Output("multi-cyto-graph",         "elements"),
    Output("multi-modules-store",      "data"),
    Output("multi-gene-input",         "value"),
    Output("multi-gene-list-version",  "data"),
    Output("multi-download-btn",       "disabled"),
    Input("multi-gene-input",          "n_blur"),
    Input("multi-example-btn",         "n_clicks"),
    Input("fdr-filter",                "value"),
    Input("multi-degree-slider",       "value"),
    State("multi-gene-input",          "value"),
    State("multi-gene-list-version",   "data"),
)
def update_multi(_n_blur, _example_btn, fdr, degree, text, version):
    new_textarea = dash.no_update
    if ctx.triggered_id == "multi-example-btn":
        text = _EXAMPLE_GENES
        new_textarea = _EXAMPLE_GENES
    pct = int(fdr * 100)
    new_version = (version or 0) + 1

    # Parse textarea: split on any whitespace, comma, semicolon, or tab so
    # users can paste from Excel (tab/newline), a paper (comma/space), or
    # type one symbol per line — all formats work without reformatting first.
    if not text or not text.strip():
        return ("Paste or type gene symbols to see their network.", [], None,
                new_textarea, new_version, True)

    seen, genes, unknown = set(), [], []
    for token in re.split(r"[\s,;]+", text):
        g = token.strip()
        if not g:
            continue
        if g in seen:
            continue
        seen.add(g)
        if g in _ALL_GENES_SET:
            genes.append(g)
        else:
            unknown.append(g)

    if len(genes) < 2:
        msg = "Enter at least 2 valid gene symbols."
        if unknown:
            msg += f"  Not found: {', '.join(unknown[:10])}{'…' if len(unknown) > 10 else ''}."
        return msg, [], None, new_textarea, new_version, True

    df = df_all[df_all["pvalue_adj"] <= fdr]
    gene_set = set(genes)
    degree = int(degree or 0)

    # Degree-of-interaction expansion — BFS out from the input genes along
    # the FDR-filtered network. Expanded genes are tagged with their hop
    # distance so they can be coloured separately from the input genes.
    if degree > 0:
        all_genes, gene_degree = _expand_by_degree(gene_set, df, degree)
    else:
        all_genes, gene_degree = gene_set, {g: 0 for g in gene_set}
    all_genes_sorted = sorted(all_genes)

    mask = df["source"].isin(all_genes) & df["target"].isin(all_genes)
    edges_df = df[mask]

    summary = f"{len(genes)} input gene(s)"
    n_added = len(all_genes) - len(gene_set)
    if n_added:
        summary += f" + {n_added} gene(s) within {degree} degree(s) of interaction"
        if len(all_genes) >= MULTI_MAX_EXPANDED_GENES:
            summary += f" (capped at {MULTI_MAX_EXPANDED_GENES})"
    summary += f" — {len(edges_df):,} co-essential pair(s) at FDR ≤ {pct}%."
    if unknown:
        summary += (f"  Not found in network: "
                    f"{', '.join(unknown[:10])}{'…' if len(unknown) > 10 else ''}.")

    # Co-essential modules — connected components of the induced sub-network,
    # numbered by size (1 = largest). Cheap (no API calls), so it's safe to
    # recompute live on every keystroke; nodes/edges are tagged data(module)
    # so the highlight callback can select a module via Cytoscape selectors.
    gene_module, modules = _detect_modules(all_genes_sorted, edges_df)
    if modules:
        summary += (f"  {len(modules)} co-essential module(s) detected "
                    f"(≥ {MULTI_MIN_MODULE_GENES} genes)")

    # Edges are coloured by the sign of the underlying correlation.
    edges = []
    for _, row in edges_df.iterrows():
        w = min(-np.log10(max(row["pvalue_adj"], 1e-300)), 10)
        src_mod = gene_module.get(row["source"], 0)
        tgt_mod = gene_module.get(row["target"], 0)
        edge_mod = src_mod if src_mod == tgt_mod else 0
        edge_color = _MULTI_POSITIVE_COLOR if row["corr_genes"] >= 0 else _MULTI_NEGATIVE_COLOR
        edges.append({"data": {"source": row["source"], "target": row["target"],
                                "weight": w, "module": edge_mod,
                                "edge_color": edge_color}})

    # Input genes with at least one pair to ANOTHER INPUT GENE are dark
    # grey; input genes with no such pair are light grey — this stays fixed
    # regardless of the degree-of-interaction slider, even if expansion adds
    # edges from that gene to extended genes. Genes added via
    # degree-of-interaction expansion are purple.
    seed_mask = edges_df["source"].isin(gene_set) & edges_df["target"].isin(gene_set)
    seed_genes_with_pairs = (set(edges_df.loc[seed_mask, "source"])
                             | set(edges_df.loc[seed_mask, "target"]))
    nodes = []
    for g in all_genes_sorted:
        if gene_degree.get(g, 0) > 0:
            bg_color = _MULTI_EXTENDED_COLOR
        elif g in seed_genes_with_pairs:
            bg_color = _MULTI_NODE_COLOR
        else:
            bg_color = _MULTI_NODE_NEUTRAL_COLOR
        nodes.append({"data": {"id": g, "bg_color": bg_color,
                               "module": gene_module.get(g, 0),
                               "degree": gene_degree.get(g, 0)}})

    has_pairs = len(edges_df) > 0
    return summary, nodes + edges, modules, new_textarea, new_version, not has_pairs


# --- Multi-gene network: GO:BP annotation of detected modules (on demand) ---
# Gated behind a button (rather than firing on every keystroke like the live
# network draw above) because each module costs one Enrichr API round-trip —
# fine for a handful of clicks, not for a per-keystroke callback.
# All significant GO:BP terms (adj. p ≤ GO_ADJ_P_THRESHOLD) are returned for
# each module — one row per term. Results are written to multi-go-rows-store;
# the filter_go_table callback below handles the cluster filter dropdown.
@app.callback(
    Output("multi-go-rows-store",       "data"),
    Output("multi-module-status",       "children"),
    Output("module-cluster-filter",     "options"),
    Output("module-cluster-filter",     "value"),
    Output("module-cluster-filter",     "disabled"),
    Input("annotate-modules-btn",       "n_clicks"),
    Input("multi-gene-list-version",    "data"),
    State("multi-modules-store",        "data"),
    prevent_initial_call=True,
)
def annotate_multi_modules(n_clicks, _version, modules):
    if ctx.triggered_id == "multi-gene-list-version":
        return [], "", [], None, True

    if not modules:
        return ([], ("No co-essential modules detected in the current network "
                     f"(need a connected component of ≥ {MULTI_MIN_MODULE_GENES} genes — "
                     "try a larger or more interconnected gene list)."),
                [], None, True)

    to_annotate = modules[:MULTI_MAX_MODULES_TO_ANNOTATE]
    rows = []
    n_queried = n_too_small = 0

    for mod in to_annotate:
        cluster_id, gene_list = mod["cluster"], mod["genes"]
        base = {"cluster": cluster_id, "cluster_size": mod["cluster_size"],
                "go_id": "", "p_value": None, "p_value_adj": None}

        if mod["cluster_size"] < MULTI_MIN_GENES_FOR_GO:
            n_too_small += 1
            rows.append({**base, "go_term": (
                f"not annotated — module has {mod['cluster_size']} genes "
                f"(GO:BP annotation needs > {MULTI_MIN_GENES_FOR_GO - 1})")})
            continue

        if len(gene_list) > MULTI_MAX_GO_GENES:
            rows.append({**base, "go_term": (
                f"skipped — {len(gene_list)} genes > {MULTI_MAX_GO_GENES} "
                "(likely a giant-component artefact)")})
            continue

        n_queried += 1
        try:
            enrichment = gp.enrichr(gene_list=gene_list, gene_sets=GO_GENE_SETS,
                                    organism=GO_ORGANISM, outdir=None)
            terms = enrichment.results.sort_values("Adjusted P-value")
        except Exception as exc:
            rows.append({**base, "go_term": f"GO:BP query failed: {exc}"})
            continue

        sig = terms[terms["Adjusted P-value"] <= GO_ADJ_P_THRESHOLD]
        if sig.empty:
            rows.append({**base, "go_term": (
                f"no significant GO:BP terms (adj. p > {GO_ADJ_P_THRESHOLD})")})
        else:
            nr = _weighted_set_cover(sig)
            n_removed = len(sig) - len(nr)
            for _, t in nr.iterrows():
                name, go_id = _split_go_term(t["Term"])
                rows.append({**base,
                              "go_term":     name,
                              "go_id":       go_id,
                              "p_value":     float(t["P-value"]),
                              "p_value_adj": float(t["Adjusted P-value"]),
                              "n_removed":   n_removed})

    total_removed = sum(r.get("n_removed", 0) for r in rows if isinstance(r.get("n_removed"), int))
    status_bits = [f"Queried GO:BP for {n_queried} of {len(modules)} co-essential module(s)"]
    if total_removed:
        status_bits.append(f"{total_removed} redundant term(s) removed by Weighted Set Cover "
                           f"(threshold {int(GO_COVER_THRESHOLD * 100)}%)")
    if n_too_small:
        status_bits.append(f"{n_too_small} too small to annotate "
                           f"(≤ {MULTI_MIN_GENES_FOR_GO - 1} genes — see table)")
    if len(modules) > len(to_annotate):
        status_bits.append(f"capped at {MULTI_MAX_MODULES_TO_ANNOTATE} modules per click")
    status = " — ".join(status_bits) + "."

    cluster_options = (
        [{"label": "All modules", "value": "all"}] +
        [{"label": f"Module {m['cluster']}  ({m['cluster_size']} genes)", "value": m["cluster"]}
         for m in to_annotate]
    )

    return rows, status, cluster_options, "all", False


# --- Multi-gene network: filter GO:BP table by selected cluster ---
@app.callback(
    Output("multi-module-table",     "data"),
    Output("multi-go-download-btn",  "disabled"),
    Input("multi-go-rows-store",     "data"),
    Input("module-cluster-filter",   "value"),
)
def filter_go_table(rows, cluster_value):
    if not rows:
        return [], True
    if not cluster_value or cluster_value == "all":
        return rows, False
    filtered = [r for r in rows if r["cluster"] == cluster_value]
    return filtered, not bool(filtered)


# --- Multi-gene network: highlight a module on row click ---
@app.callback(
    Output("multi-cyto-graph",         "stylesheet"),
    Input("multi-module-table",        "active_cell"),
    Input("multi-gene-list-version",   "data"),
    State("multi-module-table",        "data"),
    prevent_initial_call=True,
)
def highlight_module(active_cell, _version, table_data):
    # Gene list, FDR, or degree-of-interaction changed — module numbering is
    # stale, so drop any highlight; the user must re-run GO:BP annotation
    # before highlighting a (new) module again.
    if ctx.triggered_id == "multi-gene-list-version":
        return MULTI_STYLESHEET
    if not active_cell or not table_data:
        return MULTI_STYLESHEET
    row = table_data[active_cell["row"]]
    return MULTI_STYLESHEET + _module_highlight_rules(row["cluster"])


# --- Single-gene: download all co-essential partners as CSV ---
@app.callback(
    Output("single-download", "data"),
    Input("single-download-btn", "n_clicks"),
    State("gene-dropdown", "value"),
    State("fdr-filter",    "value"),
    prevent_initial_call=True,
)
def download_single_partners(_n_clicks, gene, fdr):
    if not gene:
        return dash.no_update
    df = df_all[df_all["pvalue_adj"] <= fdr]
    mask = (df["source"] == gene) | (df["target"] == gene)
    sub = df[mask].copy()
    sub["partner"] = sub.apply(
        lambda r: r["target"] if r["source"] == gene else r["source"], axis=1
    )
    sub = (sub[["partner", "pvalue", "pvalue_adj", "corr_genes"]]
           .sort_values("pvalue_adj")
           .rename(columns={
               "partner":    "Partner Gene",
               "pvalue":     "GLS P-value",
               "pvalue_adj": "GLS Adj. P-value (FDR)",
               "corr_genes": "Correlation",
           }))
    sub["Correlation Notes"] = sub["Correlation"].apply(
        lambda c: "co-essential" if c > 0 else "anti-correlated"
    )
    pct = int(fdr * 100)
    return dcc.send_data_frame(sub.to_csv, f"{gene}_coessential_partners_FDR{pct}pct.csv",
                               index=False)


# --- Multi-gene: download co-essential pairs as CSV ---
# Mirrors the network shown on screen, including any genes added via the
# degree-of-interaction slider. "Pair Type" + the two degree columns let
# users tell input-gene pairs apart from pairs involving expanded genes.
@app.callback(
    Output("multi-download", "data"),
    Input("multi-download-btn",  "n_clicks"),
    State("multi-gene-input",    "value"),
    State("fdr-filter",          "value"),
    State("multi-degree-slider", "value"),
    prevent_initial_call=True,
)
def download_multi_pairs(_n_clicks, text, fdr, degree):
    if not text or not text.strip():
        return dash.no_update
    seen, genes = set(), []
    for token in re.split(r"[\s,;]+", text):
        g = token.strip()
        if g and g not in seen and g in _ALL_GENES_SET:
            seen.add(g)
            genes.append(g)
    if len(genes) < 2:
        return dash.no_update

    df = df_all[df_all["pvalue_adj"] <= fdr]
    gene_set = set(genes)
    degree = int(degree or 0)
    if degree > 0:
        all_genes, gene_degree = _expand_by_degree(gene_set, df, degree)
    else:
        all_genes, gene_degree = gene_set, {g: 0 for g in gene_set}

    mask = df["source"].isin(all_genes) & df["target"].isin(all_genes)
    out = df[mask][["source", "target", "pvalue", "pvalue_adj", "corr_genes"]].copy()
    out["Gene Degree of Interaction"] = out["source"].map(gene_degree)
    out["Partner Gene Degree of Interaction"] = out["target"].map(gene_degree)
    out["Pair Type"] = np.where(
        (out["Gene Degree of Interaction"] == 0) & (out["Partner Gene Degree of Interaction"] == 0),
        "Original input pair", "Degree-of-interaction pair",
    )
    out = (out.sort_values("pvalue_adj")
           .rename(columns={
               "source":     "Gene",
               "target":     "Partner Gene",
               "pvalue":     "GLS P-value",
               "pvalue_adj": "GLS Adj. P-value (FDR)",
               "corr_genes": "Correlation",
           }))
    out["Correlation Notes"] = out["Correlation"].apply(
        lambda c: "co-essential" if c > 0 else "anti-correlated"
    )
    out = out[["Gene", "Partner Gene", "GLS P-value", "GLS Adj. P-value (FDR)",
               "Correlation", "Correlation Notes",
               "Gene Degree of Interaction", "Partner Gene Degree of Interaction", "Pair Type"]]
    pct = int(fdr * 100)
    return dcc.send_data_frame(out.to_csv, f"coessential_pairs_FDR{pct}pct.csv", index=False)


# --- Multi-gene: download GO:BP enrichment results as CSV ---
# Mirrors the table currently shown (respects the cluster filter dropdown
# beside this button) and adds each module's gene list for context.
@app.callback(
    Output("multi-go-download", "data"),
    Input("multi-go-download-btn", "n_clicks"),
    State("multi-module-table",    "data"),
    State("multi-modules-store",   "data"),
    prevent_initial_call=True,
)
def download_multi_go_terms(_n_clicks, rows, modules):
    if not rows:
        return dash.no_update
    gene_lists = {m["cluster"]: "; ".join(m["genes"]) for m in (modules or [])}
    out = pd.DataFrame(rows)
    out["genes"] = out["cluster"].map(gene_lists)
    out = out.rename(columns={
        "cluster":      "Cluster",
        "cluster_size": "Module Size",
        "genes":        "Genes",
        "go_id":        "GO:ID",
        "go_term":      "GO:BP Term",
        "p_value":      "P-value",
        "p_value_adj":  "Adj. P-value (FDR)",
    })
    out = out[["Cluster", "Module Size", "Genes", "GO:ID", "GO:BP Term", "P-value", "Adj. P-value (FDR)"]]
    return dcc.send_data_frame(out.to_csv, "coessential_modules_GO_BP.csv", index=False)


# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    app.run(debug=True)
