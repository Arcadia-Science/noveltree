#!/usr/bin/env python3
"""Root an inferred species-tree topology on an explicit monophyletic outgroup."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import Phylo


def canonical(name: str) -> str:
    return name.strip().replace("_", "-").replace(" ", "-")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-tree", required=True, type=Path)
    parser.add_argument("--outgroups", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--qc-output", required=True, type=Path)
    args = parser.parse_args()

    tree = Phylo.read(args.species_tree, "newick")
    tips = {canonical(tip.name): tip for tip in tree.get_terminals()}
    outgroups = {canonical(name) for name in args.outgroups.split(",") if name.strip()}
    missing = sorted(outgroups - set(tips))
    if missing:
        raise ValueError("Outgroup taxa absent from MiniNJ tree: " + ", ".join(missing))
    if not outgroups or len(outgroups) == len(tips):
        raise ValueError("Outgroups must be a nonempty proper subset of species-tree tips")

    matching = []
    all_taxa = set(tips)
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        side = {canonical(tip.name) for tip in clade.get_terminals()}
        if side == outgroups:
            matching.append(clade)
        elif all_taxa - side == outgroups:
            matching.append(clade)
    if not matching:
        raise ValueError("Explicit outgroups do not define an edge in the MiniNJ topology")

    target = min(matching, key=lambda clade: len(clade.get_terminals()))
    target_side = {canonical(tip.name) for tip in target.get_terminals()}
    if target_side != outgroups:
        # The represented side is the ingroup; root on it produces the same
        # edge while preserving the requested complementary outgroup split.
        targets = target.get_terminals()
    else:
        targets = [tips[name] for name in sorted(outgroups)]
    tree.root_with_outgroup(targets, outgroup_branch_length=0.0)
    tree.rooted = True

    root_children = list(tree.root.clades)
    if len(root_children) != 2:
        raise ValueError("Outgroup rooting did not produce a bifurcating root")
    observed = [
        {canonical(tip.name) for tip in child.get_terminals()}
        for child in root_children
    ]
    if outgroups not in observed:
        raise ValueError("Outgroup root split failed post-root validation")
    Phylo.write(tree, args.output, "newick")

    with args.qc_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerow(("species_tree_tips", len(tips)))
        writer.writerow(("outgroup_tips", len(outgroups)))
        writer.writerow(("outgroups", ",".join(sorted(outgroups))))
        writer.writerow(("root_split_validated", "true"))


if __name__ == "__main__":
    main()
