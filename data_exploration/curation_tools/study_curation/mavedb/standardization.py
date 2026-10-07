"""Declarative term mappings and metadata validation for MaveDB curation."""

from __future__ import annotations

import copy
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import pandera as pa
from curation_tools.perturbseq_anndata_schema import _TREATMENT_TYPE_LABELS, ObsSchema
from curation_tools.study_curation.paths import (
    DEFAULT_MAPPING_PATH,
    DEFAULT_OVERRIDE_PATH,
)
from IPython.display import display
from pandera.errors import SchemaError, SchemaErrors

MAPPING_COLUMNS = {
    "source_field",
    "target_field",
    "source_value",
    "target_value",
    "target_is_missing",
    "dataset_id_prefix",
    "mapping_order",
}
OVERRIDE_COLUMNS = {
    "rule_id",
    "phase",
    "action",
    "scope_type",
    "scope_column",
    "scope_value",
    "target_column",
    "value",
    "value_type",
    "reason",
}
LOWERCASE_LABEL_COLUMNS = {
    # These fields have no finite allowed-label set in perturbseq_anndata_schema.py.
    "cell_line_label",
    "cell_type_label",
    "disease_label",
    "tissue_label",
    "treatment_label",
}
QC_COLUMNS = [
    "source_field",
    "target_field",
    "unmapped_value",
    "row_count",
    "issue_type",
]


def _read_resource(path: str | Path, required_columns: set[str]) -> pd.DataFrame:
    resource = pd.read_csv(path, dtype="string", keep_default_na=False)
    missing_columns = required_columns.difference(resource.columns)
    if missing_columns:
        raise ValueError(
            f"{path} is missing required columns: {sorted(missing_columns)}"
        )
    return resource


def _as_missing_boolean(values: pd.Series) -> pd.Series:
    return values.str.strip().str.lower().eq("true")


def _mapping_mask(
    source_values: pd.Series,
    dataset_ids: pd.Series,
    mapping_group: pd.DataFrame,
) -> pd.Series:
    """Match source values and any dataset-prefix restriction for a batch."""
    source_case_sensitive = (
        "source_case_sensitive" in mapping_group
        and _as_missing_boolean(mapping_group["source_case_sensitive"]).iloc[0]
    )
    normalized_values = source_values.astype("string")
    normalized_terms = mapping_group["source_value"].astype("string")
    if not source_case_sensitive:
        normalized_values = normalized_values.str.casefold()
        normalized_terms = normalized_terms.str.casefold()
    mask = source_values.notna() & normalized_values.isin(normalized_terms)
    prefix = str(mapping_group["dataset_id_prefix"].iloc[0])
    if prefix:
        mask &= dataset_ids.astype("string").str.startswith(prefix, na=False)
    return mask


def _lowercase_label_fields(metadata: pd.DataFrame) -> pd.DataFrame:
    """Lowercase labels not enumerated in the schema, preserving controlled terms."""
    for column in LOWERCASE_LABEL_COLUMNS:
        if column in metadata.columns:
            metadata[column] = metadata[column].astype("string").str.lower()
    return metadata


@lru_cache(maxsize=1)
def _controlled_label_spellings() -> dict[str, dict[str, str]]:
    """Read canonical label spelling from schema enums and token validation."""
    spellings = {
        "treatment_type_label": {
            label.casefold(): label for label in _TREATMENT_TYPE_LABELS
        }
    }
    for column, field in ObsSchema.to_schema().columns.items():
        if not column.endswith("_label"):
            continue
        for check in field.checks:
            if check.name == "isin":
                spellings[column] = {
                    label.casefold(): label
                    for label in check.statistics["allowed_values"]
                }
    return spellings


def _canonicalize_controlled_label_fields(metadata: pd.DataFrame) -> pd.DataFrame:
    """Normalize recognized labels, preserving unknown values and missing cells."""
    for column, spellings in _controlled_label_spellings().items():
        if column not in metadata.columns:
            continue
        values = metadata[column].astype("string")
        replacements = {}
        for value in values.dropna().unique():
            if column == "treatment_type_label":
                replacements[value] = "|".join(
                    spellings.get(token.casefold(), token) for token in value.split("|")
                )
            else:
                replacements[value] = spellings.get(value.casefold(), value)
        metadata[column] = values.map(replacements).astype("string")
    return metadata


def _dataset_prefixes_overlap(first: str, second: str) -> bool:
    """Return whether two dataset ID prefix scopes can match the same ID."""
    if not first or not second:
        return True
    return first.startswith(second) or second.startswith(first)


def _validate_mapping_order(mappings: pd.DataFrame) -> None:
    """Require source rewrites to follow matching mapping consumers."""
    normalized_source_values = mappings["source_value"].astype("string").str.casefold()
    rewrites_value = _as_missing_boolean(mappings["target_is_missing"]) | mappings[
        "source_value"
    ].ne(mappings["target_value"])
    self_mappings = mappings.loc[
        mappings["source_field"].eq(mappings["target_field"]) & rewrites_value
    ]
    other_mappings = mappings.loc[mappings["source_field"].ne(mappings["target_field"])]

    for _, rewrite in self_mappings.iterrows():
        consumers = other_mappings.loc[
            other_mappings["source_field"].eq(rewrite["source_field"])
            & normalized_source_values.loc[other_mappings.index].eq(
                str(rewrite["source_value"]).casefold()
            )
        ]
        for _, consumer in consumers.iterrows():
            rewrite_prefix = str(rewrite["dataset_id_prefix"])
            consumer_prefix = str(consumer["dataset_id_prefix"])
            if not _dataset_prefixes_overlap(rewrite_prefix, consumer_prefix):
                continue
            if rewrite["mapping_order"] <= consumer["mapping_order"]:
                rewrite_scope = rewrite_prefix or "<all datasets>"
                consumer_scope = consumer_prefix or "<all datasets>"
                raise ValueError(
                    f"Mapping for {rewrite['source_field']} value "
                    f"{rewrite['source_value']!r} rewrites the source at "
                    f"mapping_order={rewrite['mapping_order']} before its "
                    f"mapping to {consumer['target_field']} at "
                    f"mapping_order={consumer['mapping_order']}; dataset ID "
                    f"prefix scopes overlap ({rewrite_scope!r}, "
                    f"{consumer_scope!r})"
                )


def apply_term_mappings(
    metadata: pd.DataFrame,
    mapping_path: str | Path = DEFAULT_MAPPING_PATH,
    inplace: bool = False,
    include_rule_audit: bool = False,
    collect_qc: bool = True,
) -> (
    tuple[pd.DataFrame, pd.DataFrame] | tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
):
    """Apply ordered mappings and optionally collect term and missing-ID QC.

    ``standardize_metadata`` defers QC until its final overrides and label
    normalization are complete. Direct callers receive QC for the mapped data.
    """
    mappings = _read_resource(mapping_path, MAPPING_COLUMNS)
    mappings["mapping_order"] = pd.to_numeric(mappings["mapping_order"])
    mappings = mappings.sort_values("mapping_order", kind="stable")
    _validate_mapping_order(mappings)
    result = metadata if inplace else metadata.copy()
    dataset_ids = result.get(
        "dataset_id", pd.Series(pd.NA, index=result.index, dtype="string")
    )
    rule_audit_rows: list[dict[str, object]] = []

    group_columns = [
        "mapping_order",
        "source_field",
        "target_field",
        "dataset_id_prefix",
    ]
    if "source_case_sensitive" in mappings:
        group_columns.append("source_case_sensitive")
    for group_key, group in mappings.groupby(group_columns, sort=False, dropna=False):
        _, source_field, target_field, dataset_id_prefix = group_key[:4]
        has_rule_audit = "rule_id" in group and group["rule_id"].ne("").any()
        if source_field not in result:
            if has_rule_audit and target_field not in result:
                result[target_field] = pd.NA
            if include_rule_audit and has_rule_audit:
                rule_audit_rows.extend(_mapping_group_audit(group, None, None))
            continue

        if target_field not in result:
            result[target_field] = pd.NA
        source_values = result[source_field].astype("string")
        mask = _mapping_mask(source_values, dataset_ids, group)
        if include_rule_audit and has_rule_audit:
            rule_audit_rows.extend(_mapping_group_audit(group, source_values, mask))

        if mask.any():
            source_case_sensitive = (
                "source_case_sensitive" in group
                and _as_missing_boolean(group["source_case_sensitive"]).iloc[0]
            )
            normalized_source_values = (
                source_values if source_case_sensitive else source_values.str.casefold()
            )
            mapped_source_values = group["source_value"].astype("string")
            if not source_case_sensitive:
                mapped_source_values = mapped_source_values.str.casefold()
            output_values = dict(
                zip(
                    mapped_source_values,
                    zip(
                        group["target_value"],
                        _as_missing_boolean(group["target_is_missing"]),
                    ),
                )
            )
            mapped = normalized_source_values.loc[mask].map(
                {
                    source_value: (pd.NA if is_missing else target_value)
                    for source_value, (
                        target_value,
                        is_missing,
                    ) in output_values.items()
                }
            )
            if not pd.api.types.is_object_dtype(result[target_field].dtype):
                result[target_field] = result[target_field].astype("object")
            result.loc[mask, target_field] = mapped

    standardized = _lowercase_label_fields(result)
    unmapped = (
        _collect_metadata_qc(standardized, mappings)
        if collect_qc
        else pd.DataFrame(columns=QC_COLUMNS)
    )
    if include_rule_audit:
        rule_audit = pd.DataFrame(
            rule_audit_rows,
            columns=[
                "rule_id",
                "phase",
                "action",
                "target_column",
                "rows_affected",
                "reason",
            ],
        )
        return standardized, unmapped, rule_audit
    return standardized, unmapped


def _mapping_group_audit(
    group: pd.DataFrame,
    source_values: pd.Series | None,
    group_mask: pd.Series | None,
) -> list[dict[str, object]]:
    """Create audit rows for direct mapping-table rules."""
    audit_rows: list[dict[str, object]] = []
    if "rule_id" not in group:
        return audit_rows
    for _, rule in group.loc[group["rule_id"].ne("")].iterrows():
        missing_target = str(rule["target_is_missing"]).strip().lower() == "true"
        rows_affected = 0
        if source_values is not None and group_mask is not None:
            case_sensitive = (
                "source_case_sensitive" in rule.index
                and str(rule["source_case_sensitive"]).strip().lower() == "true"
            )
            values = source_values if case_sensitive else source_values.str.casefold()
            rule_value = str(rule["source_value"])
            if not case_sensitive:
                rule_value = rule_value.casefold()
            rows_affected = int((group_mask & values.eq(rule_value)).sum())
        audit_rows.append(
            {
                "rule_id": str(rule.get("rule_id", "")),
                "phase": "mapping",
                "action": "clear" if missing_target else "set",
                "target_column": str(rule["target_field"]),
                "rows_affected": rows_affected,
                "reason": str(rule.get("reason", "")),
            }
        )
    return audit_rows


def _override_mask(metadata: pd.DataFrame, rule: pd.Series) -> pd.Series:
    scope_type = str(rule["scope_type"])
    if scope_type == "all":
        return pd.Series(True, index=metadata.index)
    scope_column = str(rule["scope_column"])
    if scope_column not in metadata:
        return pd.Series(False, index=metadata.index)
    values = metadata[scope_column].astype("string")
    scope_value = str(rule["scope_value"])
    case_insensitive = scope_column.endswith("_label") or scope_column == "species"
    comparison_values = values.str.casefold() if case_insensitive else values
    comparison_scope = scope_value.casefold() if case_insensitive else scope_value
    if scope_type == "equals":
        return comparison_values.eq(comparison_scope).fillna(False)
    if scope_type == "not_equals":
        return ~comparison_values.eq(comparison_scope).fillna(False)
    if scope_type == "prefix":
        return comparison_values.str.startswith(comparison_scope, na=False)
    if scope_type == "missing":
        return metadata[scope_column].isna() | values.eq("").fillna(False)
    raise ValueError(f"Unsupported override scope type: {scope_type}")


def _override_value(rule: pd.Series) -> Any:
    value_type = str(rule["value_type"])
    value = str(rule["value"])
    if value_type == "missing":
        return pd.NA
    if value_type == "integer":
        return int(value)
    if value_type == "string":
        return value
    raise ValueError(f"Unsupported override value type: {value_type}")


def apply_metadata_overrides(
    metadata: pd.DataFrame,
    phase: str,
    override_path: str | Path = DEFAULT_OVERRIDE_PATH,
    inplace: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply ordered, auditable rules from the MaveDB override table."""
    overrides = _read_resource(override_path, OVERRIDE_COLUMNS)
    result = metadata if inplace else metadata.copy()
    audit_rows: list[dict[str, object]] = []

    for _, rule in overrides.loc[overrides["phase"].eq(phase)].iterrows():
        action = str(rule["action"])
        mask = _override_mask(result, rule)
        affected = int(mask.sum())
        target_column = str(rule["target_column"])
        if action == "drop_rows":
            if inplace:
                result.drop(index=result.index[mask], inplace=True)
            else:
                result = result.loc[~mask].copy()
        elif action in {"set", "clear", "fill_missing"}:
            if not target_column:
                raise ValueError(f"Rule {rule['rule_id']} requires a target column")
            if target_column not in result:
                result[target_column] = pd.NA
            if action == "fill_missing":
                current = result[target_column].astype("string")
                mask &= result[target_column].isna() | current.eq("").fillna(False)
                affected = int(mask.sum())
            value = pd.NA if action == "clear" else _override_value(rule)
            if not pd.api.types.is_object_dtype(result[target_column].dtype):
                result[target_column] = result[target_column].astype("object")
            result.loc[mask, target_column] = value
        else:
            raise ValueError(f"Unsupported override action: {action}")

        audit_rows.append(
            {
                "rule_id": str(rule["rule_id"]),
                "phase": phase,
                "action": action,
                "target_column": target_column,
                "rows_affected": affected,
                "reason": str(rule["reason"]),
            }
        )

    return result, pd.DataFrame(
        audit_rows,
        columns=[
            "rule_id",
            "phase",
            "action",
            "target_column",
            "rows_affected",
            "reason",
        ],
    )


def _qc_tokens(value: object, field: str) -> list[str]:
    if pd.isna(value) or str(value) == "":
        return []
    text = str(value)
    return text.split("|") if field.startswith("treatment_") else [text]


def _qc_normalized_value(value: str, field: str) -> str:
    spellings = _controlled_label_spellings().get(field, {})
    tokens = [
        spellings.get(token.casefold(), token) for token in _qc_tokens(value, field)
    ]
    normalized = "|".join(tokens)
    return normalized.lower() if field in LOWERCASE_LABEL_COLUMNS else normalized


def _qc_case_sensitive(rule: dict[str, Any]) -> bool:
    return str(rule.get("source_case_sensitive", "false")).strip().lower() == "true"


def _qc_mapping_missing(rule: dict[str, Any]) -> bool:
    return (
        str(rule["target_is_missing"]).strip().lower() == "true"
        or rule["target_value"] == ""
    )


def _qc_scope_mask(metadata: pd.DataFrame, rule: dict[str, Any]) -> pd.Series:
    dataset_ids = metadata.get(
        "dataset_id", pd.Series(pd.NA, index=metadata.index, dtype="string")
    ).astype("string")
    mask = pd.Series(True, index=metadata.index)
    prefix = str(rule["dataset_id_prefix"])
    if prefix:
        mask &= dataset_ids.str.startswith(prefix, na=False)
    if rule["source_field"] == "dataset_id":
        values = dataset_ids if _qc_case_sensitive(rule) else dataset_ids.str.casefold()
        term = str(rule["source_value"])
        mask &= values.eq(term if _qc_case_sensitive(rule) else term.casefold()).fillna(
            False
        )
    return mask


def _qc_term_mask(values: pd.Series, term: str, case_sensitive: bool) -> pd.Series:
    compared = values if case_sensitive else values.str.casefold()
    return compared.eq(term if case_sensitive else term.casefold()).fillna(False)


def _qc_id_aliases(
    rule: dict[str, Any], rewrites: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Carry an ID lookup's expectation through scoped label renaming chains."""
    aliases = [rule]
    seen = {(rule["source_value"], rule["dataset_id_prefix"], _qc_case_sensitive(rule))}
    for alias in aliases:
        for rewrite in rewrites:
            term = str(alias["source_value"])
            source = str(rewrite["source_value"])
            if not _qc_case_sensitive(rewrite):
                term, source = term.casefold(), source.casefold()
            if term != source:
                continue
            first = str(alias["dataset_id_prefix"])
            second = str(rewrite["dataset_id_prefix"])
            if not _dataset_prefixes_overlap(first, second):
                continue
            prefix = first if len(first) >= len(second) else second
            value = _qc_normalized_value(
                str(rewrite["target_value"]), str(rule["source_field"])
            )
            key = (value, prefix, True)
            if not value or key in seen:
                continue
            seen.add(key)
            aliases.append(
                {
                    **rule,
                    "source_value": value,
                    "dataset_id_prefix": prefix,
                    "source_case_sensitive": "true",
                }
            )
    return aliases


def _qc_field_values(
    metadata: pd.DataFrame, field: str, id_column: str
) -> pd.DataFrame:
    """Expand compact metadata into corresponding label/ID token positions."""
    indices: list[int] = []
    values: list[str] = []
    missing_ids: list[bool] = []
    for index, row in metadata.iterrows():
        ids = _qc_tokens(row.get(id_column, pd.NA), id_column)
        for position, value in enumerate(_qc_tokens(row[field], field)):
            indices.append(index)
            values.append(value)
            missing_ids.append(position >= len(ids) or ids[position] == "")
    result = metadata.loc[indices].reset_index(drop=True)
    result["_qc_row"] = indices
    result["_qc_value"] = pd.Series(values, dtype="string")
    result["_qc_id_missing"] = missing_ids
    return result


def _qc_issue_rows(
    values: pd.DataFrame,
    mask: pd.Series,
    field: str,
    target: str,
    issue_type: str,
) -> list[dict[str, Any]]:
    # Repeated treatment tokens count once per affected metadata row and term.
    affected = values.loc[mask].drop_duplicates(["_qc_row", "_qc_value"])
    counts = affected.groupby("_qc_value", sort=False)["_qc_row_count"].sum()
    return [
        {
            "source_field": field,
            "target_field": target,
            "unmapped_value": value,
            "row_count": int(count),
            "issue_type": issue_type,
        }
        for value, count in counts.items()
    ]


def _collect_metadata_qc(
    metadata: pd.DataFrame,
    mappings: pd.DataFrame,
    overrides: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Check final terms and expected IDs independently of mapping execution."""
    if metadata.empty:
        return pd.DataFrame(columns=QC_COLUMNS)
    rules = mappings.sort_values("mapping_order", kind="stable").to_dict("records")
    override_rules = [] if overrides is None else overrides.to_dict("records")
    allowed = {
        field: set(spellings.values())
        for field, spellings in _controlled_label_spellings().items()
    }
    for field, column in ObsSchema.to_schema().columns.items():
        for check in column.checks:
            if check.name == "isin":
                allowed[field] = set(check.statistics["allowed_values"])
    source_fields = {str(rule["source_field"]) for rule in rules} - {"dataset_id"}
    fields = (
        source_fields
        | {
            str(rule["target_field"])
            for rule in rules
            if str(rule["target_field"]).endswith("_label")
        }
        | {field for field in allowed if field.endswith("_label")}
    )
    fields |= {
        str(rule["target_column"])
        for rule in override_rules
        if str(rule["target_column"]).endswith("_label")
    }
    issues: list[dict[str, Any]] = []
    for field in sorted(fields):
        if field not in metadata:
            if field in source_fields:
                prefixes = {
                    str(rule["dataset_id_prefix"])
                    for rule in rules
                    if rule["source_field"] == field
                }
                row_count = len(metadata) if "" in prefixes else 0
                if "" not in prefixes and "dataset_id" in metadata:
                    row_count = int(
                        metadata["dataset_id"]
                        .astype("string")
                        .str.startswith(tuple(prefixes), na=False)
                        .sum()
                    )
                if row_count:
                    issues.append(
                        {
                            "source_field": field,
                            "target_field": field,
                            "unmapped_value": "<source column missing>",
                            "row_count": row_count,
                            "issue_type": "missing_source_column",
                        }
                    )
            continue
        id_column = (
            field.removesuffix("_label") + "_id" if field.endswith("_label") else ""
        )
        scope_columns = {
            str(rule["scope_column"])
            for rule in override_rules
            if rule["target_column"] in {field, id_column}
        }
        columns = sorted(
            ({field, id_column, "dataset_id"} | scope_columns) & set(metadata.columns)
        )
        # Group one field and its scope/ID columns at a time. This bounds memory
        # for score tables with millions of repeated metadata rows, and retains
        # their exact counts without copying the full metadata table.
        compact = (
            metadata.groupby(columns, dropna=False, sort=False, observed=True)
            .size()
            .reset_index(name="_qc_row_count")
        )
        compact[columns] = compact[columns].astype("string")
        values = _qc_field_values(compact, field, id_column)
        if values.empty:
            continue
        terms = values["_qc_value"]
        recognized = pd.Series(False, index=values.index)
        expected_id = pd.Series(field == "treatment_type_label", index=values.index)
        rewrites = [
            rule
            for rule in rules
            if rule["source_field"] == field
            and rule["target_field"] == field
            and not _qc_mapping_missing(rule)
        ]

        def apply_id_overrides(phase: str) -> None:
            for override in override_rules:
                if override["phase"] != phase or override["target_column"] != id_column:
                    continue
                if override["action"] not in {"set", "clear", "fill_missing"}:
                    continue
                mask = _override_mask(values, pd.Series(override))
                intentional_missing = (
                    override["action"] == "clear"
                    or override["value_type"] == "missing"
                    or override["value"] == ""
                )
                expected_id.loc[mask] = not intentional_missing

        apply_id_overrides("pre_mapping")
        for rule in rules:
            source = str(rule["source_field"])
            target = str(rule["target_field"])
            missing = _qc_mapping_missing(rule)
            if source == field:
                scope = _qc_scope_mask(values, rule)
                for term in _qc_tokens(rule["source_value"], field):
                    recognized |= scope & _qc_term_mask(
                        terms, term, _qc_case_sensitive(rule)
                    )
            if target == field and not missing:
                scope = _qc_scope_mask(values, rule)
                normalized = _qc_normalized_value(str(rule["target_value"]), field)
                for term in _qc_tokens(normalized, field):
                    recognized |= scope & _qc_term_mask(terms, term, True)
            if target != id_column or source not in {field, "dataset_id"}:
                continue
            candidates = (
                [rule] if source == "dataset_id" else _qc_id_aliases(rule, rewrites)
            )
            for candidate in candidates:
                mask = _qc_scope_mask(values, candidate)
                if source != "dataset_id":
                    matches = pd.Series(False, index=values.index)
                    for term in _qc_tokens(candidate["source_value"], field):
                        matches |= _qc_term_mask(
                            terms, term, _qc_case_sensitive(candidate)
                        )
                    mask &= matches
                expected_id.loc[mask] = not missing
        apply_id_overrides("post_mapping")
        for override in override_rules:
            if override["target_column"] != field or override["action"] not in {
                "set",
                "fill_missing",
            }:
                continue
            if override["value_type"] == "missing":
                continue
            scope = _override_mask(values, pd.Series(override))
            normalized = _qc_normalized_value(str(override["value"]), field)
            for term in _qc_tokens(normalized, field):
                recognized |= scope & _qc_term_mask(terms, term, True)
        if field in allowed:
            # A source alias is not a valid final controlled label. Check the
            # actual canonical spelling against the schema after normalization.
            recognized = terms.isin(allowed[field])
        issues.extend(
            _qc_issue_rows(values, ~recognized, field, field, "unrecognized_label")
        )
        if id_column:
            issues.extend(
                _qc_issue_rows(
                    values,
                    recognized & expected_id & values["_qc_id_missing"],
                    field,
                    id_column,
                    "missing_expected_id",
                )
            )
    return pd.DataFrame(issues, columns=QC_COLUMNS)


def standardize_metadata(
    metadata: pd.DataFrame,
    mapping_path: str | Path = DEFAULT_MAPPING_PATH,
    override_path: str | Path = DEFAULT_OVERRIDE_PATH,
    inplace: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run pre-mapping corrections, term mappings, then final corrections.

    Set ``inplace=True`` to avoid duplicating large combined MaveDB metadata tables.
    """
    result = metadata if inplace else metadata.copy()
    before_mapping, before_audit = apply_metadata_overrides(
        result, "pre_mapping", override_path, inplace=True
    )
    mapped, _, mapping_audit = apply_term_mappings(
        before_mapping,
        mapping_path,
        inplace=True,
        include_rule_audit=True,
        collect_qc=False,
    )
    standardized, after_audit = apply_metadata_overrides(
        mapped, "post_mapping", override_path, inplace=True
    )
    _canonicalize_controlled_label_fields(standardized)
    _lowercase_label_fields(standardized)
    mappings = _read_resource(mapping_path, MAPPING_COLUMNS)
    mappings["mapping_order"] = pd.to_numeric(mappings["mapping_order"])
    overrides = _read_resource(override_path, OVERRIDE_COLUMNS)
    unmapped = _collect_metadata_qc(standardized, mappings, overrides)
    audit = pd.concat([before_audit, mapping_audit, after_audit], ignore_index=True)
    return standardized, unmapped, audit


def filter_curated_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Retain rows marked with a curation agent, matching the prior notebook."""
    if "curation_agent_type" not in metadata.columns:
        raise ValueError("Metadata is missing the curation_agent_type column")
    return metadata.loc[metadata["curation_agent_type"].notna()].copy()


def validate_with_unique_values(
    schema: pa.DataFrameSchema,
    dataframe: pd.DataFrame,
    n_failure_cases: int | None = 1,
    **validate_kwargs: Any,
) -> pd.DataFrame:
    """Validate a dataframe and display a compact summary of failed values."""
    if n_failure_cases is not None and n_failure_cases < 1:
        raise ValueError("n_failure_cases must be positive or None")

    schema_to_validate = (
        schema.to_schema() if hasattr(schema, "to_schema") else copy.deepcopy(schema)
    )
    checks = list(schema_to_validate.checks)
    for column in schema_to_validate.columns.values():
        checks.extend(column.checks)
    if schema_to_validate.index is not None:
        checks.extend(getattr(schema_to_validate.index, "checks", []))
        for index in getattr(schema_to_validate.index, "indexes", []):
            checks.extend(index.checks)
    for check in checks:
        check.n_failure_cases = n_failure_cases

    try:
        return schema_to_validate.validate(dataframe, **validate_kwargs)
    except (SchemaError, SchemaErrors) as exc:
        errors = exc.schema_errors if isinstance(exc, SchemaErrors) else [exc]
        data_for_rendering = (
            exc.data
            if isinstance(getattr(exc, "data", None), pd.DataFrame)
            else dataframe
        )
        row_indices: list[Any] = []
        errors_by_index: dict[Any, list[str]] = {}
        for error in errors:
            failure_cases = error.failure_cases
            if (
                not hasattr(failure_cases, "drop_duplicates")
                or "index" not in failure_cases.columns
            ):
                continue
            check = error.check
            check_label = (
                getattr(check, "error", None)
                or getattr(check, "name", None)
                or str(check)
            )
            for row_index in failure_cases["index"].dropna().drop_duplicates():
                if row_index not in errors_by_index:
                    row_indices.append(row_index)
                    errors_by_index[row_index] = []
                errors_by_index[row_index].append(check_label)

        if row_indices:
            error_rows = data_for_rendering.loc[row_indices].copy()
            error_rows.insert(
                0,
                "_validation_errors",
                [
                    "; ".join(dict.fromkeys(errors_by_index.get(index, [])))
                    for index in error_rows.index
                ],
            )
            display(error_rows)

        reduced_errors = []
        for error in errors:
            failure_cases = error.failure_cases
            if (
                hasattr(failure_cases, "drop_duplicates")
                and "failure_case" in failure_cases.columns
            ):
                key_columns = [
                    column
                    for column in (
                        "schema_context",
                        "column",
                        "check",
                        "check_number",
                        "failure_case",
                    )
                    if column in failure_cases.columns
                ]
                failure_cases = failure_cases.drop_duplicates(
                    subset=key_columns, keep="first"
                )
                unique_values = failure_cases["failure_case"].tolist()
                prefix, separator, _ = str(error).partition("failure cases:")
                message = (
                    f"{prefix}{separator} {', '.join(map(repr, unique_values))}"
                    if separator
                    else str(error)
                )
            else:
                message = str(error)

            reduced_errors.append(
                SchemaError(
                    schema=error.schema,
                    data=getattr(error, "data", None),
                    message=message,
                    failure_cases=failure_cases,
                    check=error.check,
                    check_index=error.check_index,
                    check_output=error.check_output,
                    parser=error.parser,
                    parser_index=getattr(error, "parser_index", None),
                    parser_output=error.parser_output,
                    reason_code=error.reason_code,
                    column_name=error.column_name,
                )
            )

        if isinstance(exc, SchemaErrors):
            raise SchemaErrors(
                schema=exc.schema,
                schema_errors=reduced_errors,
                data=exc.data,
            ) from exc
        raise reduced_errors[0] from exc
