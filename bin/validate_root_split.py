#!/usr/bin/env python3
"""Assert that two rooted Newick trees encode the same root bipartition."""

import argparse
import csv
from pathlib import Path

from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", required=True, type=Path)
    parser.add_argument("--observed", required=True, type=Path)
    parser.add_argument("--qc-output", required=True, type=Path)
    return parser.parse_args()


def tip_names(tree, label):
    names = [tip.name for tip in tree.get_terminals()]
    if any(not name for name in names):
        raise ValueError(f"{label} contains an unnamed tip")
    if len(names) != len(set(names)):
        raise ValueError(f"{label} contains duplicate tip names")
    return set(names)


def root_partition(tree, label):
    children = list(tree.root.clades)
    if len(children) != 2:
        raise ValueError(
            f"{label} must have exactly two root children; found {len(children)}"
        )
    sides = [frozenset(tip.name for tip in child.get_terminals()) for child in children]
    if not all(sides):
        raise ValueError(f"{label} has an empty root side")
    return frozenset(sides)


def validate(expected_path, observed_path, qc_path):
    expected = Phylo.read(expected_path, "newick")
    observed = Phylo.read(observed_path, "newick")
    expected_tips = tip_names(expected, "Expected rooted tree")
    observed_tips = tip_names(observed, "Observed rooted tree")
    if expected_tips != observed_tips:
        raise ValueError(
            "Final SpeciesRax tip set changed: "
            f"missing={sorted(expected_tips - observed_tips)[:10]}, "
            f"extra={sorted(observed_tips - expected_tips)[:10]}"
        )
    expected_partition = root_partition(expected, "Expected rooted tree")
    observed_partition = root_partition(observed, "Observed rooted tree")
    if expected_partition != observed_partition:
        raise ValueError(
            "Final SpeciesRax root split differs from the rooted MiniNJ input"
        )

    existing = []
    if qc_path.exists():
        with qc_path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            existing = list(reader)
    if not existing:
        existing = [["metric", "value"]]
    elif existing[0] != ["metric", "value"]:
        raise ValueError(f"Invalid root QC header: {qc_path}")
    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(existing)
        writer.writerow(("final_speciesrax_tips", len(observed_tips)))
        writer.writerow(("final_root_split_validated", "true"))


def main():
    args = parse_args()
    validate(args.expected, args.observed, args.qc_output)


if __name__ == "__main__":
    main()
