#!/usr/bin/env python3
"""Check every UPR feature target resolves through the pipeline's GTF rules."""

import argparse
import importlib.util
import re
from pathlib import Path
import sys


LEGACY_IDS = {
    "AARS": "ENSG00000090861",
    "ATP5B": "ENSG00000110955",
    "CARS": "ENSG00000110619",
    "DARS": "ENSG00000115866",
    "HARS": "ENSG00000170445",
    "MARS": "ENSG00000166986",
    "QARS": "ENSG00000172053",
    "SARS": "ENSG00000031698",
    "SLMO2": "ENSG00000101166",
    "SRPR": "ENSG00000182934",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True)
    parser.add_argument(
        "--features",
        default=str(
            Path(__file__).with_name("adamson_2016_upr_perturb_seq_features.tsv")
        ),
    )
    parser.add_argument(
        "--comparison",
        default=str(Path(__file__).parents[2] / "comparison" / "comparison.py"),
    )
    args = parser.parse_args()

    comparison_path = Path(args.comparison).resolve()
    sys.path.insert(0, str(comparison_path.parent))
    sys.argv = [
        str(comparison_path),
        "--dataset-id",
        "adamson_2016_upr_perturb_seq",
        "--curated-h5ad",
        "unused.h5ad",
        "--reprocessed-h5ad",
        "unused.h5ad",
        "--gtf",
        args.gtf,
    ]
    spec = importlib.util.spec_from_file_location("comparison", comparison_path)
    comparison = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(comparison)
    comparison.GUIDE_TARGET_ENSG_BY_SYMBOL = comparison.load_gtf_gene_ensgs_by_symbol(
        args.gtf
    )

    by_target = {}
    labels_checked = 0
    with open(args.features) as handle:
        for row_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 2, f"Malformed feature row {row_number}"
            label = fields[1]
            target = comparison.guide_target_name(label)
            if target == comparison.CONTROL_TARGET_SYMBOL:
                continue
            labels_checked += 1

            explicit_ids = re.findall(r"ENSG\d+(?:\.\d+)?", label)
            assert len(explicit_ids) <= 1, f"Multiple ENSGs in {label}"
            ensg = comparison.guide_target_ensg(label)
            assert re.fullmatch(r"ENSG\d+", ensg), f"Unresolved target: {label}"
            if explicit_ids:
                assert target in LEGACY_IDS, f"Unexpected explicit ENSG in {label}"
                assert (
                    ensg == LEGACY_IDS[target]
                ), f"Wrong stable ID for {label}: {ensg}"
            by_target.setdefault(target, set()).add(ensg)

    assert set(LEGACY_IDS) <= set(by_target), "A legacy target is missing from features"
    ambiguous = {target: ids for target, ids in by_target.items() if len(ids) != 1}
    assert not ambiguous, f"Targets resolve to multiple ENSGs: {ambiguous}"
    print(f"Validated {labels_checked} guide labels across {len(by_target)} targets")


if __name__ == "__main__":
    main()
