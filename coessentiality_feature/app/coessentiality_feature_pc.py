# =============================================================================
# IMPORTS
# =============================================================================
import json
import logging
import os
import re
import numpy as np
import pandas as pd
import dash
import dash_cytoscape as cyto
from dash import dcc, html, dash_table, Input, Output, State, ctx
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


def _load_all_profiled_genes():
    """All genes that passed NA filtering and entered the GLS analysis — the
    full searchable gene set, not just genes with a significant partner in the
    FDR 10% network (which would silently exclude profiled-but-isolated genes)."""
    version_raw = open(os.path.join(_DATA_DIR, "depmap_version.txt")).read().strip()
    m = re.search(r"(\d+Q\d+)", version_raw, re.IGNORECASE)
    version = m.group(1) if m else None

    if version:
        genes_path = os.path.join(_DATA_DIR, f"depmap_{version}_genes.txt")
        if os.path.isfile(genes_path):
            with open(genes_path) as fh:
                return sorted(line.strip() for line in fh if line.strip())

    # Fallback: derive from the network CSV alone (excludes isolated genes).
    return sorted(set(df_all["source"]).union(set(df_all["target"])))


# Full profiled gene set — populates both search dropdowns. Used separately
# from df_all, which only carries gene *pairs/relationships*.
all_genes = _load_all_profiled_genes()


def _load_dataset_stats():
    """Return (n_cell_lines, n_genes_profiled) from the metadata JSON written
    by step2_gls_coessentiality.py — avoids scanning the full CRISPRGeneEffect
    CSV (hundreds of MB) just to count rows at app startup."""
    version_raw = open(os.path.join(_DATA_DIR, "depmap_version.txt")).read().strip()
    m = re.search(r"(\d+Q\d+)", version_raw, re.IGNORECASE)
    version = m.group(1) if m else None

    if version:
        metadata_path = os.path.join(_DATA_DIR, f"depmap_{version}_metadata.json")
        if os.path.isfile(metadata_path):
            with open(metadata_path) as fh:
                metadata = json.load(fh)
            return metadata.get("n_cell_lines"), metadata.get("n_genes_profiled")

    return None, None


_N_CELL_LINES, _N_GENES_PROFILED = _load_dataset_stats()
_DEPMAP_VERSION = open(os.path.join(_DATA_DIR, "depmap_version.txt")).read().strip()


# =============================================================================
# LOGGING
# Server-side only — error details (incl. tracebacks) go to stdout/stderr,
# which Cloud Run captures automatically. Never shown to users; UI gets a
# generic message instead (see annotate_multi_modules).
# =============================================================================
logger = logging.getLogger(__name__)


# =============================================================================
# APP INITIALISATION
# =============================================================================
cyto.load_extra_layouts()

app = dash.Dash(__name__)
app.title = "DepMap Co-Essentiality Explorer"

# DataTable wraps markdown cells (the Direction badge / Statistical
# Confidence bar) in a <p> with the browser's default margin, which adds
# unwanted vertical spacing inside the cell. Injected here via index_string
# rather than a separate assets/ CSS file, since it's the only global rule
# needed and this keeps everything in one file.
app.index_string = app.index_string.replace(
    "</head>",
    "<style>#partner-table .dash-cell-value p { margin: 0; }</style></head>",
)


# =============================================================================
# COLOUR HELPER
# Maps a co-essentiality direction (+1/-1) to a hex colour on a
# blue (#0072B2) → neutral grey (#dcdcdc) → orange (#E69F00) gradient.
# Palette: Wong (2011) colorblind-safe — distinguishable under deuteranopia,
# protanopia, and tritanopia. Direction is intentionally kept on this palette
# rather than the Perturbation Catalogue brand green/red, since green+red is
# the classic problematic combination for red-green colourblindness.
# The three-stop gradient is computed in Python because Cytoscape's mapData
# only interpolates between two colours, producing a muddy midpoint.
# =============================================================================
_NEUTRAL  = np.array([220, 220, 220])
_POSITIVE = np.array([  0, 114, 178])   # #0072B2 — blue  (positive co-essentiality)
_NEGATIVE = np.array([230, 159,   0])   # #E69F00 — orange (negative co-essentiality)

def direction_to_color(c):
    c = float(np.clip(c, -1.0, 1.0))
    # +1 → blue, 0 → neutral grey, -1 → orange
    rgb = (_NEUTRAL + c * (_POSITIVE - _NEUTRAL)).astype(int) if c >= 0 \
          else (_NEUTRAL + (-c) * (_NEGATIVE - _NEUTRAL)).astype(int)
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


# Flat (non-gradient) colours for the multi-gene network:
#  - input genes are a bright, near-black grey — the most prominent nodes
#  - genes added via degree-of-interaction expansion are a light, dull
#    purple — distinct in hue from the grey above
#  - edges are coloured by the underlying co-essentiality direction, using the
#    same Positive/Negative colours as the single-gene view and the table
# Grey/purple shades are distinguished by lightness/hue independent of the
# blue/orange edge colours, so the palette stays colour-blind friendly.
_MULTI_NODE_COLOR         = "#222222"   # near-black — input gene
_MULTI_POSITIVE_COLOR     = "#0072B2"   # blue   — positive-direction edges
_MULTI_NEGATIVE_COLOR     = "#E69F00"   # orange — negative-direction edges
_MULTI_EXTENDED_COLOR     = "#D8BFD8"   # light purple — added by degree of interaction (dull)


# =============================================================================
# CYTOSCAPE STYLESHEETS
# Nodes use data(bg_color) — a hex string computed per-node in the callback —
# so colour reflects the co-essentiality direction without needing mapData colours.
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
# numbered by size (1 = largest). Detection is cheap and runs live as the user
# edits their gene list; GO:BP annotation runs a local hypergeometric test per
# module (no network calls), gated behind an explicit button to avoid
# recomputing on every keystroke.
# =============================================================================
MULTI_MAX_DEGREE              = 2     # max degree-of-interaction slider value — beyond this the
                                       # network grows too large to interpret
MULTI_MAX_EXPANDED_GENES      = 500   # safety ceiling on total rendered network size (raw pasted input
                                       # AND degree-of-interaction expansion) — only kicks in when
                                       # exceeded; truncation keeps the most statistically significant
                                       # genes (PR #365 comment 9), not an arbitrary/random subset
MULTI_MIN_MODULE_GENES        = 2      # smallest connected component counted as a "module"
MULTI_MIN_GENES_FOR_GO        = 4      # only modules with MORE than 3 genes are GO:BP-annotated
                                       # (a 2-3 gene module is too small for a meaningful enrichment test)
MULTI_MAX_MODULES_TO_ANNOTATE = 10    # cap on enrichment runs per click (latency, not API politeness anymore)
MULTI_MAX_GO_GENES            = 200   # modules larger than this are reported but not queried
GO_GENE_SET_PATH              = os.path.join(_DATA_DIR, "GO_Biological_Process_2025.gmt")
GO_ADJ_P_THRESHOLD            = 0.05  # only GO:BP terms at or below this FDR are shown
GO_COVER_THRESHOLD            = 0.50  # WSC redundancy threshold: a term is redundant if
                                      # ≥50% of its genes are already covered by a more
                                      # significant selected term (mirrors WebGestalt default)

# Local GO:BP gene-set library — downloaded once via:
#   gp.get_library(name="GO_Biological_Process_2025", organism="Human",
#                   save="required_data/GO_Biological_Process_2025.gmt")
# Pinned to a specific Enrichr-packaged version (these are refreshed roughly
# every 1-2 years) for reproducibility; bump deliberately, not automatically.
# Enrichment runs fully offline against this — no per-query network calls,
# no dependency on Enrichr's uptime/bandwidth.
_GO_LIBRARY = gp.get_library(name=GO_GENE_SET_PATH)

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


def _truncate_genes_by_significance(gene_set, df, cap):
    """Truncate gene_set down to at most `cap` genes if it exceeds the cap —
    deterministically, keeping genes that participate in the most
    statistically significant pairs (smallest adjusted p-value) within the
    set first, rather than an arbitrary/random subset (PR #365 comment 9).
    No-op if gene_set is already within the cap. Returns (kept, truncated).
    """
    if len(gene_set) <= cap:
        return gene_set, False

    mask = df["source"].isin(gene_set) & df["target"].isin(gene_set)
    induced = df.loc[mask, ["source", "target", "pvalue_adj"]].sort_values("pvalue_adj")

    kept = set()
    for src, tgt, _ in induced.itertuples(index=False):
        new = {g for g in (src, tgt) if g not in kept}
        if not new:
            continue
        if len(kept) + len(new) > cap:
            break
        kept |= new

    if len(kept) < cap:
        # Fill any remaining room with leftover genes that have no
        # significant edge within this set at all (e.g. isolated inputs),
        # in deterministic alphabetical order.
        leftover = sorted(gene_set - kept)
        kept |= set(leftover[: cap - len(kept)])

    return kept, True


def _expand_by_degree(seed_genes, df, max_degree):
    """BFS-expand a gene set along co-essential edges, up to max_degree hops.

    seed_genes: the user's input genes (degree 0).
    df:         FDR-filtered network to search for neighbours (full df_all,
                not just the induced sub-network).

    Returns (all_genes, gene_degree, truncated) where gene_degree maps
    gene -> hop distance from the seed set (0 for seed genes), and
    truncated is True if MULTI_MAX_EXPANDED_GENES was reached and some
    candidate genes had to be dropped.
    """
    gene_degree = {g: 0 for g in seed_genes}
    all_genes = set(seed_genes)
    frontier = set(seed_genes)
    truncated = False

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
            # Rank candidates by their best (smallest) adjusted p-value edge
            # to the current frontier — deterministic and principled,
            # instead of list(set)[:room]'s hash-order-dependent arbitrary cut.
            candidate_mask = mask & (df["source"].isin(new_genes) | df["target"].isin(new_genes))
            best_p = {}
            for src, tgt, p in df.loc[candidate_mask, ["source", "target", "pvalue_adj"]].itertuples(index=False):
                for g in (src, tgt):
                    if g in new_genes and (g not in best_p or p < best_p[g]):
                        best_p[g] = p
            ranked = sorted(new_genes, key=lambda g: (best_p.get(g, 1.0), g))
            new_genes = set(ranked[:room])
            truncated = True
        for g in new_genes:
            gene_degree[g] = d
        all_genes |= new_genes
        frontier = new_genes

    return all_genes, gene_degree, truncated


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


def _partner_table_columns():
    """Columns for the co-essential partners table."""
    return [
        {"name": "PARTNER GENE", "id": "partner"},
        {"name": "DIRECTION", "id": "direction_badge", "presentation": "markdown"},
        {"name": "GLS P-VALUE", "id": "pvalue", "type": "numeric",
         "format": {"specifier": ".2e"}},
        {"name": "STATISTICAL CONFIDENCE (BH-ADJUSTED P-VALUE)",
         "id": "strength_bar", "presentation": "markdown"},
    ]


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

                        # Co-essential partners table — single table, one row per
                        # gene (replaces the former side-by-side bar chart +
                        # table, which showed the same gene list twice with no
                        # guaranteed row alignment between them, and squeezed
                        # the table into only half the page width). Bar length
                        # encodes statistical confidence (BH-adjusted), kept
                        # distinct from biological effect size, which GLS does
                        # not report (see PR #365 review discussion). Direction
                        # badge replaces a raw +/-1.0 number; a legend above
                        # the table spells out what it means without needing
                        # to read the paragraph below.
                        html.Div([
                            dcc.Loading(html.Div([
                                html.P(
                                    "Co-essential partners — significance & direction",
                                    style={"fontWeight": "600", "fontSize": "20px",
                                           "color": _TEXT, "margin": "0 0 8px 0"}),
                                html.P(
                                    "Each row is a gene whose CRISPR essentiality profile "
                                    "co-varies significantly with the query gene across cancer "
                                    "cell lines. Bar length shows statistical confidence after "
                                    "Benjamini-Hochberg multiple-testing correction (the grey "
                                    "number is −log₁₀ of the adjusted p-value — higher means more "
                                    "confident) — a longer bar means we are more confident the "
                                    "relationship is real, not necessarily that it is a bigger "
                                    "biological effect.",
                                    style={"fontSize": "18px", "color": _MUTED,
                                           "margin": "0 0 12px 0", "lineHeight": "1.5"}),
                                html.Div([
                                    html.Div("Direction:", style={
                                        "fontSize": "16px", "fontWeight": "700",
                                        "color": _TEXT, "marginBottom": "8px"}),
                                    html.Div([
                                        html.Span("+ Positive", style={
                                            "display": "inline-block", "padding": "3px 11px",
                                            "borderRadius": "13px", "fontSize": "15px",
                                            "fontWeight": "600", "color": "#fff",
                                            "whiteSpace": "nowrap", "marginRight": "8px",
                                            "background": "#0072B2"}),
                                        html.Span(
                                            "co-essential gene pair",
                                            style={"fontSize": "16px", "color": _MUTED}),
                                    ], style={"display": "flex", "alignItems": "center",
                                              "marginBottom": "6px"}),
                                    html.Div([
                                        html.Span("− Negative", style={
                                            "display": "inline-block", "padding": "3px 11px",
                                            "borderRadius": "13px", "fontSize": "15px",
                                            "fontWeight": "600", "color": "#fff",
                                            "whiteSpace": "nowrap", "marginRight": "8px",
                                            "background": "#E69F00"}),
                                        html.Span(
                                            "anti-correlated gene pair",
                                            style={"fontSize": "16px", "color": _MUTED}),
                                    ], style={"display": "flex", "alignItems": "center"}),
                                ], style={"marginBottom": "16px"}),
                                html.Div(
                                    "Please select a gene from the Search gene box above "
                                    "to see its co-essential partners.",
                                    id="partner-table-placeholder",
                                    style={"display": "none", "textAlign": "center",
                                           "padding": "48px 20px", "fontSize": "18px",
                                           "color": _MUTED, "background": _BG,
                                           "borderRadius": "6px"},
                                ),
                                html.Div(dash_table.DataTable(
                                    id="partner-table",
                                    columns=_partner_table_columns(),
                                    markdown_options={"html": True},
                                    sort_action="native",
                                    page_size=15,
                                    style_table={"overflowX": "auto", "border": "none"},
                                    style_cell={
                                        "textAlign": "left", "padding": "10px 12px",
                                        "fontSize": "18px", "border": "none",
                                        "borderBottom": f"1px solid {_BORDER}",
                                        "fontFamily": "inherit",
                                        "verticalAlign": "middle",
                                    },
                                    style_cell_conditional=[
                                        {"if": {"column_id": "partner"}, "width": "16%"},
                                        {"if": {"column_id": "direction_badge"}, "width": "14%"},
                                        {"if": {"column_id": "pvalue"}, "width": "14%",
                                         "textAlign": "right"},
                                        {"if": {"column_id": "strength_bar"}, "width": "56%"},
                                    ],
                                    style_header={
                                        "fontWeight": "700", "fontSize": "16px",
                                        "color": _MUTED, "letterSpacing": "0.05em",
                                        "border": "none",
                                        "borderBottom": f"2px solid {_BORDER}",
                                        "background": "#fff",
                                    },
                                    style_header_conditional=[
                                        {"if": {"column_id": "pvalue"}, "textAlign": "right"},
                                    ],
                                    style_data_conditional=[{
                                        "if": {"state": "selected"},
                                        "backgroundColor": _G_LITE,
                                        "border": f"1px solid {_G}",
                                    }],
                                ), id="partner-table-wrapper"),
                            ]), type="circle", color=_G),
                        ], style={"marginBottom": "32px"}),

                        # Network panel
                        _panel("Co-essential genetic interaction network", [
                            html.P(
                                "Partner nodes coloured by "
                                "direction (edge thickness = −log₁₀(adj. p-value)). ",
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
                                    "Edges coloured by co-essentiality direction "
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
                                            "fontWeight": "700", "marginRight": "20px"}),
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


# Inline styles (not a separate CSS file) for the Direction badge and
# Statistical Confidence bar rendered inside DataTable markdown cells.
# Colours use the Wong (2011) colourblind-safe blue/orange pair — kept
# deliberately separate from the Perturbation Catalogue brand green/red,
# since green+red is the classic problematic combination for red-green
# colourblindness.
_BADGE_STYLE = "display:inline-block;padding:3px 11px;border-radius:13px;font-size:15px;font-weight:600;color:#fff;white-space:nowrap;"
_BADGE_POS_STYLE = _BADGE_STYLE + "background:#0072B2;"
_BADGE_NEG_STYLE = _BADGE_STYLE + "background:#E69F00;"


def _direction_badge_html(direction_value):
    if direction_value >= 0:
        return f'<span style="{_BADGE_POS_STYLE}">+ Positive</span>'
    return f'<span style="{_BADGE_NEG_STYLE}">&minus; Negative</span>'


def _strength_bar_html(pvalue_adj, max_strength):
    # Fill is a light/semi-transparent tint (not solid colour) so the dark
    # label text reads fine whether it lands on the filled or unfilled part
    # of the bar — same trick as tskir's mockup screenshot. The -log10(FDR)
    # scale itself is shown once, as axis tick numbers in the column header
    # (see _partner_table_columns), not repeated per row.
    strength = -np.log10(max(pvalue_adj, 1e-300))
    pct = min(100, (strength / max_strength) * 100) if max_strength > 0 else 0
    track_style = "position:relative;height:30px;background:#eef2ef;border-radius:3px;overflow:hidden;min-width:160px;"
    fill_style = (f"position:absolute;top:0;left:0;bottom:0;width:{pct:.1f}%;"
                  "background:rgba(0,123,83,0.35);border-radius:3px 0 0 3px;"
                  "display:flex;align-items:center;")
    label_style = "color:#212121;font-size:18px;font-weight:600;padding-left:8px;white-space:nowrap;"
    return (
        f'<div style="{track_style}"><div style="{fill_style}">'
        f'<span style="{label_style}">{pvalue_adj:.2e}</span></div></div>'
    )


_PLACEHOLDER_STYLE_BASE = {
    "textAlign": "center", "padding": "48px 20px", "fontSize": "18px",
    "color": _MUTED, "background": _BG, "borderRadius": "6px",
}
_PLACEHOLDER_SHOWN  = {**_PLACEHOLDER_STYLE_BASE, "display": "block"}
_PLACEHOLDER_HIDDEN = {**_PLACEHOLDER_STYLE_BASE, "display": "none"}
_TABLE_SHOWN  = {"display": "block"}
_TABLE_HIDDEN = {"display": "none"}


# --- Single-gene explorer ---
@app.callback(
    Output("summary-text",              "children"),
    Output("partner-table",             "data"),
    Output("partner-table",             "columns"),
    Output("cyto-graph",                "elements"),
    Output("single-download-btn",       "disabled"),
    Output("partner-table-placeholder", "children"),
    Output("partner-table-placeholder", "style"),
    Output("partner-table-wrapper",     "style"),
    Input("gene-dropdown",  "value"),
    Input("fdr-filter",     "value"),
)
def update_single(gene, fdr):
    df = df_all[df_all["pvalue_adj"] <= fdr]
    pct = int(fdr * 100)

    if not gene:
        placeholder_text = ("Please select a gene from the Search gene box above "
                             "to see its co-essential partners.")
        return ("Select a gene to explore its co-essential partners.", [], _partner_table_columns(), [], True,
                placeholder_text, _PLACEHOLDER_SHOWN, _TABLE_HIDDEN)

    mask = (df["source"] == gene) | (df["target"] == gene)
    sub = df[mask].copy()
    sub["partner"] = sub.apply(
        lambda r: r["target"] if r["source"] == gene else r["source"], axis=1
    )
    sub = sub[["partner", "pvalue", "pvalue_adj", "direction"]].sort_values("pvalue_adj")

    if sub.empty:
        summary = f"No co-essential partners found for {gene} at FDR ≤ {pct}%."
        query_node = [{"data": {"id": gene, "bg_color": "#D55E00"}, "classes": "query"}]
        placeholder_text = f"No co-essential partners found for {gene} at FDR ≤ {pct}%."
        return (summary, [], _partner_table_columns(), query_node, True,
                placeholder_text, _PLACEHOLDER_SHOWN, _TABLE_HIDDEN)

    summary = f"{gene} has {len(sub):,} co-essential partner(s) at FDR ≤ {pct}%."

    # Bar length is scaled relative to this gene's own strongest partner, so
    # the scale stays meaningful across every page of this gene's table.
    max_strength = -np.log10(max(sub["pvalue_adj"].min(), 1e-300))
    sub["direction_badge"] = sub["direction"].apply(_direction_badge_html)
    sub["strength_bar"] = sub["pvalue_adj"].apply(_strength_bar_html, max_strength=max_strength)

    top = sub.head(20)
    top_partners = set(top["partner"].tolist())

    # direction per partner (used for node colour)
    partner_direction = sub.set_index("partner")["direction"].to_dict()

    # Query node: orange via .query selector; bg_color unused but set for consistency
    nodes = [{"data": {"id": gene, "bg_color": "#D55E00"}, "classes": "query"}]
    for p in top_partners:
        nodes.append({"data": {"id": p, "bg_color": direction_to_color(partner_direction.get(p, 0.0))}})

    # Spoke edges: query → each partner (fixed mid-weight so they're visible)
    edges = [{"data": {"source": gene, "target": p, "weight": 3.0}}
             for p in top_partners]

    # Cross-edges among partners, weighted by significance
    partner_mask = df["source"].isin(top_partners) & df["target"].isin(top_partners)
    for _, row in df[partner_mask].iterrows():
        w = min(-np.log10(max(row["pvalue_adj"], 1e-300)), 10)
        edges.append({"data": {"source": row["source"], "target": row["target"],
                                "weight": w}})

    return (summary, sub.to_dict("records"), _partner_table_columns(), nodes + edges, False,
            "", _PLACEHOLDER_HIDDEN, _TABLE_SHOWN)


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

    # Truncate the raw pasted input itself if it alone exceeds the ceiling —
    # previously only expansion was capped, so pasting a huge list directly
    # (degree=0) went straight to Cytoscape uncapped and could crash the
    # browser (PR #365 comment 9).
    gene_set, input_truncated = _truncate_genes_by_significance(gene_set, df, MULTI_MAX_EXPANDED_GENES)

    # Degree-of-interaction expansion — BFS out from the input genes along
    # the FDR-filtered network. Expanded genes are tagged with their hop
    # distance so they can be coloured separately from the input genes.
    if degree > 0:
        all_genes, gene_degree, expansion_truncated = _expand_by_degree(gene_set, df, degree)
    else:
        all_genes, gene_degree, expansion_truncated = gene_set, {g: 0 for g in gene_set}, False
    all_genes_sorted = sorted(all_genes)

    mask = df["source"].isin(all_genes) & df["target"].isin(all_genes)
    edges_df = df[mask]

    summary = f"{len(genes)} input gene(s)"
    n_added = len(all_genes) - len(gene_set)
    if n_added:
        summary += f" + {n_added} gene(s) within {degree} degree(s) of interaction"
    summary += f" — {len(edges_df):,} co-essential pair(s) at FDR ≤ {pct}%."
    if input_truncated or expansion_truncated:
        summary += (f"  Network truncated to the {MULTI_MAX_EXPANDED_GENES:,} most "
                    f"statistically significant genes (sorted by GLS p-value) "
                    f"to preserve performance.")

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
        edge_color = _MULTI_POSITIVE_COLOR if row["direction"] >= 0 else _MULTI_NEGATIVE_COLOR
        edges.append({"data": {"source": row["source"], "target": row["target"],
                                "weight": w, "module": edge_mod,
                                "edge_color": edge_color}})

    # All input genes share one colour, regardless of whether they have a
    # pair to another input gene specifically — whether a gene is connected
    # is already visible from the rendered graph itself, and distinguishing
    # "has a pair" vs "no pair" became actively misleading once
    # degree-of-interaction expansion could give a previously-isolated input
    # gene real partners while it kept showing the "no pair" colour. Genes
    # added via degree-of-interaction expansion are purple.
    nodes = []
    for g in all_genes_sorted:
        if gene_degree.get(g, 0) > 0:
            bg_color = _MULTI_EXTENDED_COLOR
        else:
            bg_color = _MULTI_NODE_COLOR
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
            enrichment = gp.enrich(gene_list=gene_list, gene_sets=_GO_LIBRARY,
                                    background=all_genes, outdir=None, no_plot=True)
            terms = enrichment.results.sort_values("Adjusted P-value")
        except Exception:
            logger.exception(
                "GO:BP enrichment failed for module %s (%d genes)",
                cluster_id, len(gene_list),
            )
            rows.append({**base, "go_term": "GO:BP enrichment failed — please try again."})
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
    sub = (sub[["partner", "pvalue", "pvalue_adj", "direction"]]
           .sort_values("pvalue_adj")
           .rename(columns={
               "partner":    "Partner Gene",
               "pvalue":     "GLS P-value",
               "pvalue_adj": "GLS Adj. P-value (FDR)",
               "direction":  "Direction",
           }))
    sub["Direction Notes"] = sub["Direction"].apply(
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
    # Same truncation as the on-screen network (update_multi), so the
    # download matches what's actually displayed.
    gene_set, _ = _truncate_genes_by_significance(gene_set, df, MULTI_MAX_EXPANDED_GENES)
    if degree > 0:
        all_genes, gene_degree, _ = _expand_by_degree(gene_set, df, degree)
    else:
        all_genes, gene_degree = gene_set, {g: 0 for g in gene_set}

    mask = df["source"].isin(all_genes) & df["target"].isin(all_genes)
    out = df[mask][["source", "target", "pvalue", "pvalue_adj", "direction"]].copy()
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
               "direction":  "Direction",
           }))
    out["Direction Notes"] = out["Direction"].apply(
        lambda c: "co-essential" if c > 0 else "anti-correlated"
    )
    out = out[["Gene", "Partner Gene", "GLS P-value", "GLS Adj. P-value (FDR)",
               "Direction", "Direction Notes",
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
