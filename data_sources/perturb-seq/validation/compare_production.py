#!/usr/bin/env python3
"""Compare a bounded unified-pipeline sample with live production results."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import urllib.parse
import urllib.request
from urllib.error import HTTPError
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def ordered(values: set[str]) -> list[str]:
    return sorted(values, key=lambda value: hashlib.sha256(value.encode()).hexdigest())


def scan_unique(path: str, columns: list[str]) -> dict[str, set[str]]:
    found = {column: set() for column in columns}
    parquet = pq.ParquetFile(path)
    for batch in parquet.iter_batches(columns=columns, batch_size=1_000_000):
        for column in columns:
            for value in pc.unique(batch.column(column)).to_pylist():
                if value is not None and str(value):
                    found[column].add(str(value))
    return found


def filtered_rows(
    path: str,
    columns: list[str],
    target_values: set[str],
    effect_values: set[str] | None = None,
) -> list[dict]:
    parquet = pq.ParquetFile(path)
    target_array = pa.array(sorted(target_values))
    effect_array = (
        pa.array(sorted(effect_values)) if effect_values is not None else None
    )
    rows: list[dict] = []
    for batch in parquet.iter_batches(columns=columns, batch_size=1_000_000):
        mask = pc.is_in(batch.column("perturbed_target_symbol"), value_set=target_array)
        if effect_array is not None:
            mask = pc.and_(
                mask,
                pc.is_in(batch.column("effect_gene_symbol"), value_set=effect_array),
            )
        selected = pa.Table.from_batches([batch]).filter(mask)
        rows.extend(selected.to_pylist())
    return rows


def json_request(url: str) -> object:
    request = urllib.request.Request(
        url, headers={"User-Agent": "PerturbationCatalogue-validation/1.0"}
    )
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.load(response)
        except HTTPError as exc:
            detail = exc.read().decode("utf-8", "replace")[:1000]
            if attempt == 2:
                raise RuntimeError(f"HTTP {exc.code} for {url}: {detail}") from exc
        except Exception:
            if attempt == 2:
                raise
    raise AssertionError("unreachable")


def query_dea(base_url: str, target: str) -> list[dict]:
    rows: list[dict] = []
    offset = 0
    # The production service returns HTTP 500 for some targets at 20k rows;
    # smaller pages are equivalent and keep each response bounded.
    limit = 5_000
    while True:
        print(f"Fetching DEA target {target} offset {offset}", flush=True)
        query = urllib.parse.urlencode(
            {"limit": limit, "offset": offset, "perturbation_gene_name": target}
        )
        payload = json_request(
            f"{base_url}/v1/perturb-seq/replogle_2022_k562_gw_normalized/search?{query}"
        )
        page = payload["results"]
        rows.extend(page)
        if len(rows) >= payload["total_rows_count"] or not page:
            return rows
        offset += len(page)


def query_gsea(base_url: str, target: str) -> list[dict]:
    query = urllib.parse.urlencode(
        {
            "dataset_id": "replogle_2022_k562_gw_normalized",
            "perturbation_gene_name": target,
        }
    )
    payload = json_request(f"{base_url}/v1/perturb-seq-gsea?{query}")
    return [
        effect
        | {"perturbed_target_ensg": result["perturbation"]["perturbed_target_ensg"]}
        for result in payload
        for effect in result.get("effects", [])
    ]


def ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    result = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        result[order[start:end]] = (start + end - 1) / 2 + 1
        start = end
    return result


def correlation(
    left: np.ndarray, right: np.ndarray
) -> tuple[float | None, float | None]:
    if len(left) < 2 or np.all(left == left[0]) or np.all(right == right[0]):
        return None, None
    pearson = float(np.corrcoef(left, right)[0, 1])
    spearman = float(np.corrcoef(ranks(left), ranks(right))[0, 1])
    return pearson, spearman


def paired_stats(
    left: list[float], right: list[float], transform=lambda value: value
) -> dict:
    pairs = [
        (transform(a), transform(b))
        for a, b in zip(left, right)
        if a is not None
        and b is not None
        and math.isfinite(transform(a))
        and math.isfinite(transform(b))
    ]
    if not pairs:
        return {"n": 0}
    a = np.asarray([pair[0] for pair in pairs], dtype=float)
    b = np.asarray([pair[1] for pair in pairs], dtype=float)
    pearson, spearman = correlation(a, b)
    delta = np.abs(a - b)
    return {
        "n": int(len(a)),
        "pearson": pearson,
        "spearman": spearman,
        "mae": float(np.mean(delta)),
        "median_abs_error": float(np.median(delta)),
        "max_abs_error": float(np.max(delta)),
    }


def sign_agreement(left: list[float], right: list[float]) -> float | None:
    pairs = [(a, b) for a, b in zip(left, right) if a is not None and b is not None]
    if not pairs:
        return None
    return float(np.mean([np.sign(a) == np.sign(b) for a, b in pairs]))


def dea_metrics(cluster: list[dict], production: list[dict]) -> dict:
    def cluster_key(row):
        return (row["perturbed_target_ensg"], row["effect_gene_ensg"])

    def production_key(row):
        return (
            row["perturbation"]["perturbed_target_ensg"],
            row["effect"]["effect_gene_ensg"],
        )

    cluster_map = {cluster_key(row): row for row in cluster}
    production_map = {production_key(row): row for row in production}
    common = sorted(cluster_map.keys() & production_map.keys())
    crows = [cluster_map[key] for key in common]
    prows = [production_map[key] for key in common]
    metric_map = {
        "log2foldchange": ("log2foldchange", "log2fc", lambda value: value),
        "score_value": ("score_value", "score_value", lambda value: value),
        "minus_log10_padj": (
            "padj",
            "padj",
            lambda value: -math.log10(max(value, 1e-300)),
        ),
    }
    metrics = {}
    for name, (cluster_field, production_field, transform) in metric_map.items():
        left = [row[cluster_field] for row in crows]
        right = [row["effect"][production_field] for row in prows]
        metrics[name] = paired_stats(left, right, transform)
        if name in {"log2foldchange", "score_value"}:
            metrics[name]["sign_agreement"] = sign_agreement(left, right)
    return {
        "cluster_rows": len(cluster),
        "production_rows": len(production),
        "common_rows": len(common),
        "cluster_to_production_coverage": (
            len(common) / len(cluster) if cluster else None
        ),
        "production_to_cluster_coverage": (
            len(common) / len(production) if production else None
        ),
        "duplicate_cluster_keys": len(cluster) - len(cluster_map),
        "duplicate_production_keys": len(production) - len(production_map),
        "numeric_metrics": metrics,
    }


def gsea_metrics(cluster: list[dict], production: list[dict]) -> dict:
    def key(row):
        return (row["perturbed_target_ensg"], row["term"])

    cluster_map = {key(row): row for row in cluster}
    production_map = {key(row): row for row in production}
    common = sorted(cluster_map.keys() & production_map.keys())
    crows = [cluster_map[item] for item in common]
    prows = [production_map[item] for item in common]
    metrics = {}
    for name, field, transform in (
        ("es", "es", lambda value: value),
        ("nes", "nes", lambda value: value),
        ("minus_log10_pval", "pval", lambda value: -math.log10(max(value, 1e-300))),
        ("minus_log10_fdr", "fdr", lambda value: -math.log10(max(value, 1e-300))),
    ):
        left = [row[field] for row in crows]
        right = [row[field] for row in prows]
        metrics[name] = paired_stats(left, right, transform)
        if name in {"es", "nes"}:
            metrics[name]["sign_agreement"] = sign_agreement(left, right)
    jaccards = []
    for left, right in zip(crows, prows):
        a = set(left.get("leading_edge") or [])
        b = set(right.get("leading_edge") or [])
        if a or b:
            jaccards.append(len(a & b) / len(a | b))
    metrics["leading_edge_jaccard"] = {
        "n": len(jaccards),
        "median": float(np.median(jaccards)) if jaccards else None,
        "mean": float(np.mean(jaccards)) if jaccards else None,
    }
    return {
        "cluster_rows": len(cluster),
        "production_rows": len(production),
        "common_rows": len(common),
        "cluster_to_production_coverage": (
            len(common) / len(cluster) if cluster else None
        ),
        "production_to_cluster_coverage": (
            len(common) / len(production) if production else None
        ),
        "duplicate_cluster_keys": len(cluster) - len(cluster_map),
        "duplicate_production_keys": len(production) - len(production_map),
        "numeric_metrics": metrics,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-dea", required=True)
    parser.add_argument("--cluster-gsea", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--production-url",
        default="https://perturbation-catalogue-be-asmuzum42q-nw.a.run.app",
    )
    parser.add_argument("--effect-fraction", type=float, default=0.02)
    parser.add_argument("--target-count", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dea_columns = [
        "perturbed_target_symbol",
        "perturbed_target_ensg",
        "effect_gene_symbol",
        "effect_gene_ensg",
        "padj",
        "log2foldchange",
        "score_value",
    ]
    gsea_columns = [
        "perturbed_target_symbol",
        "perturbed_target_ensg",
        "term",
        "es",
        "nes",
        "pval",
        "fdr",
        "leading_edge",
    ]
    unique = scan_unique(
        args.cluster_dea, ["perturbed_target_symbol", "effect_gene_symbol"]
    )
    target_symbols = ordered(unique["perturbed_target_symbol"])[: args.target_count]
    effect_count = max(
        1, math.ceil(len(unique["effect_gene_symbol"]) * args.effect_fraction)
    )
    effect_symbols = ordered(unique["effect_gene_symbol"])[:effect_count]
    cluster_dea = filtered_rows(
        args.cluster_dea, dea_columns, set(target_symbols), set(effect_symbols)
    )
    cluster_gsea = filtered_rows(args.cluster_gsea, gsea_columns, set(target_symbols))
    target_ensgs = {
        row["perturbed_target_ensg"]
        for row in cluster_dea
        if row["perturbed_target_ensg"]
    }
    effect_ensgs = {
        row["effect_gene_ensg"] for row in cluster_dea if row["effect_gene_ensg"]
    }

    production_dea = []
    production_gsea = []
    for target in target_symbols:
        production_dea.extend(query_dea(args.production_url.rstrip("/"), target))
        production_gsea.extend(query_gsea(args.production_url.rstrip("/"), target))
    production_dea = [
        row
        for row in production_dea
        if row["effect"]["effect_gene_ensg"] in effect_ensgs
        and row["perturbation"]["perturbed_target_ensg"] in target_ensgs
    ]
    production_gsea = [
        row for row in production_gsea if row["perturbed_target_ensg"] in target_ensgs
    ]

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset_id": "replogle_2022_k562_gw_normalized",
        "production_api": args.production_url.rstrip("/"),
        "sampling": {
            "effect_gene_fraction": args.effect_fraction,
            "effect_genes_in_full_dea": len(unique["effect_gene_symbol"]),
            "effect_genes_sampled": len(effect_symbols),
            "perturbation_targets_sampled": len(target_symbols),
            "target_symbols": target_symbols,
            "effect_symbols_sha256_ordered_sample": effect_symbols,
        },
        "dea": dea_metrics(cluster_dea, production_dea),
        "gsea": gsea_metrics(cluster_gsea, production_gsea),
    }
    json_path = outdir / "replogle-production-comparison.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n")
    markdown = [
        "# Replogle unified pipeline vs production API",
        "",
        f"Generated `{report['created_at']}` from live `{report['production_api']}`.",
        "The cluster sample uses 2% of the 12,420 effect genes (deterministically selected by SHA-256 order) across 12 deterministic perturbation targets.",
        "",
        "## Results",
        "",
        "| Product | Cluster rows | Production rows | Common keys | Cluster→production | Production→cluster |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name in ("dea", "gsea"):
        result = report[name]
        markdown.append(
            f"| {name.upper()} | {result['cluster_rows']:,} | {result['production_rows']:,} | {result['common_rows']:,} | {result['cluster_to_production_coverage']:.3%} | {result['production_to_cluster_coverage']:.3%} |"
        )
    markdown.extend(
        [
            "",
            "Detailed machine-readable metrics are in `replogle-production-comparison.json`.",
            "",
        ]
    )
    (outdir / "replogle-production-comparison.md").write_text("\n".join(markdown))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
