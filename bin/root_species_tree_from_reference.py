#!/usr/bin/env python3
"""Transfer the root split of a trusted chronogram to an unrooted species tree."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from Bio import Phylo


def canonical_taxon(name: str) -> str:
    return name.strip().replace("_", "-").replace(" ", "-")


def terminal_map(tree, label: str):
    mapped = {}
    for terminal in tree.get_terminals():
        if terminal.name is None:
            raise ValueError(f"{label} contains an unnamed tip")
        key = canonical_taxon(terminal.name)
        if key in mapped:
            raise ValueError(f"{label} contains duplicate normalized tip {key!r}")
        mapped[key] = terminal
    return mapped


def normalized_descendants(clade) -> set[str]:
    return {canonical_taxon(tip.name) for tip in clade.get_terminals()}


def root_sides(tree, shared: set[str], label: str) -> tuple[set[str], set[str]]:
    children = list(tree.root.clades)
    if len(children) != 2:
        raise ValueError(
            f"{label} must encode a bifurcating root; found {len(children)} root children"
        )
    sides = tuple(normalized_descendants(child) & shared for child in children)
    if not sides[0] or not sides[1] or sides[0] | sides[1] != shared:
        raise ValueError(f"{label} root does not partition the shared taxa")
    return sides


def assert_reference_chronogram(reference, tolerance: float) -> tuple[float, float]:
    if len(reference.root.clades) != 2:
        raise ValueError(
            "Reference chronogram must have exactly two children at its encoded root"
        )
    depths = reference.depths()
    tip_depths = [depths[tip] for tip in reference.get_terminals()]
    if any(not math.isfinite(depth) for depth in tip_depths):
        raise ValueError("Reference chronogram contains non-finite root-to-tip distances")
    maximum = max(tip_depths)
    minimum = min(tip_depths)
    allowed = max(tolerance, tolerance * max(1.0, abs(maximum)))
    if maximum - minimum > allowed:
        raise ValueError(
            "Reference tree is not ultrametric at its encoded root: "
            f"root-to-tip range={maximum - minimum:.8g}, tolerance={allowed:.8g}"
        )
    return minimum, maximum


def choose_matching_edge(tree, shared: set[str], reference_sides):
    all_taxa = set(terminal_map(tree, "Species tree"))
    candidates = []
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        full_side = normalized_descendants(clade)
        shared_side = full_side & shared
        for reference_index, reference_side in enumerate(reference_sides):
            if shared_side != reference_side:
                continue
            expected_fraction = len(reference_side) / len(shared)
            observed_fraction = len(full_side) / len(all_taxa)
            candidates.append(
                (
                    abs(observed_fraction - expected_fraction),
                    len(full_side),
                    reference_index,
                    clade,
                    full_side,
                )
            )

    if not candidates:
        sizes = ", ".join(str(len(side)) for side in reference_sides)
        raise ValueError(
            "The species-tree topology has no edge compatible with the reference "
            f"root split among {len(shared)} shared taxa (reference sides: {sizes})"
        )

    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    return candidates[0], len(candidates)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-tree", required=True, type=Path)
    parser.add_argument("--reference-chronogram", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--qc-output", required=True, type=Path)
    parser.add_argument("--ultrametric-tolerance", type=float, default=1e-5)
    args = parser.parse_args()

    species_tree = Phylo.read(args.species_tree, "newick")
    reference = Phylo.read(args.reference_chronogram, "newick")
    species_terminals = terminal_map(species_tree, "Species tree")
    reference_terminals = terminal_map(reference, "Reference chronogram")
    shared = set(species_terminals) & set(reference_terminals)
    if len(shared) < 3:
        raise ValueError(
            f"At least three shared taxa are required to transfer a root; found {len(shared)}"
        )

    min_depth, max_depth = assert_reference_chronogram(
        reference, args.ultrametric_tolerance
    )
    reference_sides = root_sides(reference, shared, "Reference chronogram")
    candidate, n_candidates = choose_matching_edge(
        species_tree, shared, reference_sides
    )
    _, _, reference_index, _, selected_taxa, = candidate

    # Root on the exact unrooted edge. The selected full-taxon side is a clade
    # in the input representation and may include taxa absent from TimeTree.
    outgroup_targets = [species_terminals[name] for name in sorted(selected_taxa)]
    species_tree.root_with_outgroup(outgroup_targets, outgroup_branch_length=0.0)
    species_tree.rooted = True

    observed_sides = root_sides(species_tree, shared, "Rooted species tree")
    if not (
        (observed_sides[0] == reference_sides[0] and observed_sides[1] == reference_sides[1])
        or (observed_sides[0] == reference_sides[1] and observed_sides[1] == reference_sides[0])
    ):
        raise ValueError("Root transfer failed post-root validation")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(species_tree, args.output, "newick")

    with args.qc_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerow(("species_tree_tips", len(species_terminals)))
        writer.writerow(("reference_chronogram_tips", len(reference_terminals)))
        writer.writerow(("shared_tips", len(shared)))
        writer.writerow(("reference_root_side_1_shared_tips", len(reference_sides[0])))
        writer.writerow(("reference_root_side_2_shared_tips", len(reference_sides[1])))
        writer.writerow(("matching_edges", n_candidates))
        writer.writerow(("selected_reference_side", reference_index + 1))
        writer.writerow(("selected_species_tree_side_tips", len(selected_taxa)))
        writer.writerow(("reference_min_root_to_tip", f"{min_depth:.12g}"))
        writer.writerow(("reference_max_root_to_tip", f"{max_depth:.12g}"))
        writer.writerow(("root_split_validated", "true"))

    selected_shared = sorted(selected_taxa & shared)
    print(
        "Transferred reference root split using "
        f"{len(shared)} shared taxa; selected side has {len(selected_taxa)} total tips "
        f"({len(selected_shared)} shared): {', '.join(selected_shared[:10])}"
        + (" ..." if len(selected_shared) > 10 else "")
    )


if __name__ == "__main__":
    main()
