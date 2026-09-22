#!/usr/bin/env python3
"""Add a pipeline-ready perturbation column to the Norman author H5AD."""

import argparse
from pathlib import Path
import shutil

import h5py
import numpy as np
import pandas as pd


CONTROL_PREFIX = "NegCtrl"
DATASET_DIR = Path(__file__).resolve().parent
GENE_MAP = DATASET_DIR / "gene_ensg.tsv"


def text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def load_gene_map():
    table = pd.read_csv(GENE_MAP, sep="\t")
    expected = {"symbol", "ensembl_gene_id"}
    if set(table.columns) != expected:
        raise ValueError(f"Expected columns {sorted(expected)} in {GENE_MAP}")
    mapping = dict(
        zip(table["symbol"].astype(str), table["ensembl_gene_id"].astype(str))
    )
    if len(mapping) != len(table):
        raise ValueError(f"Duplicate symbols in {GENE_MAP}")
    return mapping


def standard_label(value, gene_map):
    identity = text(value).strip()
    if identity.lower() in {"", "nan", "none"}:
        return ""
    pair = identity.split("__", 1)[0]
    targets = [
        target
        for target in pair.split("_")
        if target and not target.startswith(CONTROL_PREFIX)
    ]
    return (
        ";".join(
            f"{target}__{gene_map[target]}" if target in gene_map else target
            for target in targets
        )
        if targets
        else "non-targeting"
    )


def prepare(source, output):
    source = Path(source).resolve()
    output = Path(output).resolve()
    temporary = output.with_name(output.name + ".partial")
    if not source.is_file():
        raise FileNotFoundError(source)
    if output.exists() or temporary.exists():
        raise FileExistsError(output)

    gene_map = load_gene_map()
    with h5py.File(source, "r") as handle:
        obs = handle["obs"]
        categories = handle["uns"]["guide_identity_categories"][:]
        values = obs[:]
        codes = values["guide_identity"]
        labels = np.asarray(
            [
                (
                    standard_label(categories[int(code)], gene_map)
                    if 0 <= int(code) < len(categories)
                    else ""
                )
                for code in codes
            ],
            dtype=object,
        )
        fields = [(name, obs.dtype.fields[name][0]) for name in obs.dtype.names]
        fields.append(("perturbation", h5py.string_dtype(encoding="utf-8")))
        expanded = np.empty(values.shape, dtype=np.dtype(fields))
        for name in obs.dtype.names:
            expanded[name] = values[name]
        expanded["perturbation"] = labels
        chunks = obs.chunks
        compression = obs.compression
        compression_opts = obs.compression_opts

    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, temporary)
    try:
        with h5py.File(temporary, "r+") as handle:
            del handle["obs"]
            handle.create_dataset(
                "obs",
                data=expanded,
                chunks=chunks,
                compression=compression,
                compression_opts=compression_opts,
            )
        temporary.replace(output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    print(f"Wrote standardized curated H5AD to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("output")
    args = parser.parse_args()
    prepare(args.source, args.output)
