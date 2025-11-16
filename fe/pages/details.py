"""New Target Details page implementing the reusable data table component."""

from __future__ import annotations

import copy
import json
from typing import Any, Dict, List, Optional

import dash
from dash import ALL, MATCH, Input, Output, State, callback, dcc, html
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from components.target_data_table import TargetDataTable
from utils import (
    COLORS,
    fetch_dataset_rows,
    fetch_modality_datasets,
)


dash.register_page(
    __name__,
    path_template="/target/<target_name>",
    name="Target Details (Grid)",
    title="Target Details",
)

TARGET_NAME_STORE = "target-details-v2-name"
SECTION_CONFIGS = [
    {
        "id": "crispr",
        "title": "CRISPR screen data",
        "modality": "crispr-screen",
        "filter_field": "perturbation_gene_name",
    },
    {
        "id": "mave",
        "title": "MAVE data",
        "modality": "mave",
        "filter_field": "perturbation_gene_name",
    },
    {
        "id": "perturb_seq_perturbed",
        "title": "Perturb-Seq (Perturbed)",
        "modality": "perturb-seq",
        "filter_field": "perturbation_gene_name",
    },
    {
        "id": "perturb_seq_affected",
        "title": "Perturb-Seq (Affected)",
        "modality": "perturb-seq",
        "filter_field": "effect_gene_name",
    },
]
SECTION_LOOKUP = {config["id"]: config for config in SECTION_CONFIGS}
MODALITY_BUTTON_GROUP = "target-modality-buttons"
SCROLL_VISIBILITY_STORE = "target-scroll-visible"
SCROLL_TOP_STYLE = {
    "position": "fixed",
    "bottom": "1.5rem",
    "right": "1.5rem",
    "zIndex": 1050,
    "display": "none",
}

DATASET_LIMIT = 4
ROWS_PER_DATASET_LIMIT = 5
DATASET_LOAD_MORE_SIZE = 5


def layout(target_name: Optional[str] = None, **kwargs):
    """Page layout with data stores and section shells."""
    stores = [
        dcc.Store(id=TARGET_NAME_STORE, data=target_name),
        dcc.Store(id=SCROLL_VISIBILITY_STORE, data=False),
    ]
    stores.extend(
        dcc.Store(
            id={"type": "target-section-store", "section": config["id"]},
            data={"section": config["id"], "dataset_search": ""},
        )
        for config in SECTION_CONFIGS
    )

    heading = target_name or "Target"

    sections = []
    for config in SECTION_CONFIGS:
        sections.append(
            html.Section(
                [
                    html.H3(config["title"], className="fw-semibold mb-3"),
                    html.Div(
                        [
                            dcc.Input(
                                id={"type": "dataset-search", "section": config["id"]},
                                type="text",
                                placeholder="Search datasets…",
                                debounce=True,
                                className="form-control form-control-sm",
                                style={
                                    "maxWidth": "240px",
                                    "position": "absolute",
                                    "top": "2.5rem",
                                    "left": "0.75rem",
                                    "zIndex": 5,
                                    "marginLeft": "15rem",
                                },
                            ),
                            dcc.Loading(
                                html.Div(
                                    id={
                                        "type": "target-table-container",
                                        "section": config["id"],
                                    }
                                ),
                                type="default",
                                color=COLORS["primary"],
                            ),
                        ],
                        className="position-relative",
                        style={"paddingTop": "2.5rem"},
                    ),
                ],
                id=f"section-{config['id']}",
                className="mb-5",
            )
        )

    return dbc.Container(
        stores
        + [
            html.Div(id="page-top"),
            html.Div(
                [
                    html.H1(
                        f"Target: {heading}",
                        className="display-5 fw-bold mb-2 text-center",
                        style={"color": COLORS["primary"]},
                    ),
                    html.P(
                        "Explore perturbation datasets across modalities.",
                        className="text-muted text-center mb-4",
                    ),
                    html.Div(
                        id=MODALITY_BUTTON_GROUP,
                        className="d-flex flex-wrap gap-3 justify-content-center mb-4",
                    ),
                ]
            ),
            *sections,
            dbc.Button(
                "↑ Back to top",
                id="target-scroll-top",
                href="#page-top",
                color="primary",
                className="position-fixed shadow",
                style=SCROLL_TOP_STYLE,
            ),
        ],
        fluid=True,
        className="py-4 target-details-v2",
    )


def _get_triggered_id():
    """Return the parsed triggered ID for Dash versions with/without triggered_id."""
    ctx = dash.callback_context
    triggered = getattr(ctx, "triggered_id", None)
    if triggered is not None:
        return triggered
    if not ctx.triggered:
        return None
    prop_id = ctx.triggered[0]["prop_id"].split(".")[0]
    try:
        return json.loads(prop_id)
    except json.JSONDecodeError:
        return prop_id


def _empty_store(config: Dict[str, Any], error: Optional[str] = None) -> Dict[str, Any]:
    return {
        "section": config["id"],
        "modality": config["modality"],
        "title": config["title"],
        "datasets": [],
        "total_datasets_count": 0,
        "filters": {},
        "rows_per_dataset_limit": ROWS_PER_DATASET_LIMIT,
        "error": error,
        "target_name": None,
        "dataset_search": "",
    }


@callback(
    Output({"type": "target-table-container", "section": MATCH}, "children"),
    Input({"type": "target-section-store", "section": MATCH}, "data"),
)
def render_section(store_data: Optional[Dict[str, Any]]):
    """Render each section's summary plus data table."""
    triggered = _get_triggered_id()
    section_id = triggered.get("section") if isinstance(triggered, dict) else None

    if not store_data:
        return html.Div("No data available.", className="text-muted fst-italic")

    datasets = store_data.get("datasets") or []
    total_datasets = store_data.get("total_datasets_count", len(datasets))
    rows_limit = store_data.get("rows_per_dataset_limit", ROWS_PER_DATASET_LIMIT)

    summary = html.Div(
        [
            html.Span(
                f"{len(datasets)} dataset(s) shown (of {total_datasets})",
                className="fw-semibold text-muted",
            ),
            html.Span(
                f"Top {rows_limit} rows per dataset shown by default",
                className="text-muted small",
            ),
        ],
        className="d-flex flex-column flex-md-row gap-2 justify-content-between mb-3",
    )

    def control_factory(dataset_id: str, dataset_entry: Dict[str, Any]):
        if (
            not dataset_id
            or not dataset_entry.get("has_more")
            or not section_id
        ):
            return None
        return dbc.Button(
            "Load more rows",
            id={
                "type": "target-table-load-more",
                "section": section_id,
                "dataset": dataset_id,
            },
            color="secondary",
            outline=True,
            size="sm",
        )

    section_id = store_data.get("section")

    table = TargetDataTable(
        data=datasets,
        modality=store_data.get("modality", ""),
        rows_per_dataset_limit=rows_limit,
        error_message=store_data.get("error"),
        dataset_control_factory=control_factory,
        effect_gene_source="perturbation"
        if section_id == "perturb_seq_affected"
        else "effect",
    )

    return html.Div([summary, table])


@callback(
    Output(MODALITY_BUTTON_GROUP, "children"),
    [
        Input({"type": "target-section-store", "section": config["id"]}, "data")
        for config in SECTION_CONFIGS
    ],
)
def render_modality_buttons(*stores):
    """Display per-modality dataset counts as anchor buttons."""
    buttons = []
    color_map = {
        "crispr-screen": "secondary",
        "mave": "danger",
        "perturb-seq": "success",
    }

    for store, config in zip(stores, SECTION_CONFIGS):
        total = store.get("total_datasets_count") if store else None
        count_text = "…" if total is None else total
        buttons.append(
            html.A(
                dbc.Button(
                    f"{config['title']} ({count_text})",
                    id={"type": "modality-nav-button", "section": config["id"]},
                    color=color_map.get(config["modality"], "primary"),
                    outline=True,
                    className="px-3 py-2 fw-semibold",
                    style={"borderWidth": "2px"},
                ),
                href=f"#section-{config['id']}",
                className="text-decoration-none",
            )
        )

    if not buttons:
        return html.Div("Loading modality counts...", className="text-muted fst-italic")

    return buttons


@callback(
    Output(SCROLL_VISIBILITY_STORE, "data"),
    [
        Input({"type": "modality-nav-button", "section": config["id"]}, "n_clicks")
        for config in SECTION_CONFIGS
    ],
    State(SCROLL_VISIBILITY_STORE, "data"),
    prevent_initial_call=True,
)
def mark_scroll_visible(*args):
    *clicks, current = args
    if any(click and click > 0 for click in clicks):
        return True
    return current


@callback(
    Output("target-scroll-top", "style"),
    Input(SCROLL_VISIBILITY_STORE, "data"),
)
def toggle_scroll_button(visible: bool):
    style = SCROLL_TOP_STYLE.copy()
    if visible:
        style["display"] = "block"
    return style


@callback(
    Output({"type": "target-section-store", "section": MATCH}, "data"),
    Input(TARGET_NAME_STORE, "data"),
    Input({"type": "target-table-load-more", "section": MATCH, "dataset": ALL}, "n_clicks"),
    Input({"type": "dataset-search", "section": MATCH}, "value"),
    State({"type": "target-section-store", "section": MATCH}, "data"),
    State({"type": "target-section-store", "section": MATCH}, "id"),
)
def update_section_store(
    target_name: Optional[str],
    _load_more_clicks,
    search_value: Optional[str],
    store_data: Optional[Dict[str, Any]],
    store_id: Optional[Dict[str, Any]],
):
    """Hydrate or append rows for a section store."""
    if not store_id:
        raise PreventUpdate

    section_id = store_id.get("section")
    config = SECTION_LOOKUP.get(section_id)
    if not config:
        raise PreventUpdate

    if not target_name:
        return _empty_store(config, "No target specified.")

    triggered = _get_triggered_id()

    search_changed = False
    previous_query = ""
    if store_data is not None:
        previous_query = store_data.get("dataset_search", "")
        search_changed = (search_value or "") != (previous_query or "")

    should_refresh = triggered == TARGET_NAME_STORE or not store_data or store_data.get(
        "target_name"
    ) != target_name

    if isinstance(triggered, dict) and triggered.get("type") == "dataset-search":
        should_refresh = True

    if should_refresh or search_changed:
        current_search = (search_value or "") if search_changed else previous_query
        return _fetch_section_payload(config, target_name, current_search)

    if isinstance(triggered, dict) and triggered.get("type") == "target-table-load-more":
        dataset_id = triggered.get("dataset")
        if not dataset_id:
            raise PreventUpdate
        return _append_dataset_rows(store_data, dataset_id)

    return store_data


def _fetch_section_payload(
    config: Dict[str, Any],
    target_name: str,
    dataset_search: Optional[str] = None,
) -> Dict[str, Any]:
    filters = {config["filter_field"]: target_name}
    if dataset_search:
        cleaned = dataset_search.strip()
        if cleaned:
            filters["dataset_metadata"] = cleaned
    response = fetch_modality_datasets(
        config["modality"],
        filters=filters,
        dataset_limit=DATASET_LIMIT,
        rows_per_dataset_limit=ROWS_PER_DATASET_LIMIT,
    )
    datasets = response.get("datasets") or []
    for dataset in datasets:
        results = dataset.get("results") or []
        dataset["rows_per_dataset_limit"] = ROWS_PER_DATASET_LIMIT
        dataset["truncated"] = len(results) == ROWS_PER_DATASET_LIMIT
        dataset["has_more"] = dataset["truncated"]
        dataset["next_offset"] = len(results)

    return {
        "section": config["id"],
        "modality": config["modality"],
        "title": config["title"],
        "datasets": datasets,
        "total_datasets_count": response.get("total_datasets_count", 0),
        "filters": filters,
        "rows_per_dataset_limit": ROWS_PER_DATASET_LIMIT,
        "error": response.get("error"),
        "target_name": target_name,
        "dataset_search": dataset_search or "",
    }


def _append_dataset_rows(store_data: Dict[str, Any], dataset_id: str) -> Dict[str, Any]:
    updated_store = copy.deepcopy(store_data)
    datasets = updated_store.get("datasets") or []
    target_dataset = next(
        (
            entry
            for entry in datasets
            if (entry.get("dataset") or {}).get("dataset_id") == dataset_id
        ),
        None,
    )

    if not target_dataset:
        return updated_store

    current_results = target_dataset.get("results") or []
    next_offset = target_dataset.get("next_offset", len(current_results))

    target_symbol = updated_store.get("target_name")
    filters = {"perturbation_gene_name": target_symbol} if target_symbol else {}

    response = fetch_dataset_rows(
        updated_store.get("modality", ""),
        dataset_id,
        filters=filters,
        offset=next_offset,
        limit=DATASET_LOAD_MORE_SIZE,
    )

    new_results = response.get("results") or []
    if not new_results:
        target_dataset["has_more"] = False
        if response.get("error"):
            updated_store["error"] = response["error"]
        return updated_store

    target_dataset["results"] = current_results + new_results
    total_count = response.get("total_rows_count") or len(target_dataset["results"])
    target_dataset["has_more"] = len(target_dataset["results"]) < total_count
    target_dataset["truncated"] = False
    target_dataset["next_offset"] = next_offset + len(new_results)

    return updated_store

