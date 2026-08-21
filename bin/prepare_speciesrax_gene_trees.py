#!/usr/bin/env python3
"""Resolve and validate local SpeciesRax gene-tree copies.

The input trees are extracted copies inside the monolithic SpeciesRax task.
This script never edits the stored gene-tree outputs produced upstream.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

from resolve_polytomies import parse_newick, resolve_polytomies, to_newick


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest-glob",
        default="speciesrax_inputs_*.tsv",
        help="Glob selecting the extracted SpeciesRax shard manifests",
    )
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def walk(node):
    yield node
    for child in node["children"]:
        yield from walk(child)


def leaf_labels(node):
    return Counter(
        current["label"]
        for current in walk(node)
        if not current["children"]
    )


def nonbinary_nodes(node):
    return [
        len(current["children"])
        for current in walk(node)
        if current["children"] and len(current["children"]) != 2
    ]


def resolve_tree(tree_path):
    original = tree_path.read_text().strip()
    if not original:
        raise ValueError(f"Empty gene tree: {tree_path}")
    if not original.endswith(";"):
        raise ValueError(f"Gene tree lacks a terminal semicolon: {tree_path}")

    tree = parse_newick(original)
    leaves_before = leaf_labels(tree)
    invalid_before = nonbinary_nodes(tree)
    unary_before = [degree for degree in invalid_before if degree < 2]
    if unary_before:
        raise ValueError(
            f"Gene tree contains {len(unary_before)} unary internal nodes, "
            f"which cannot be resolved as polytomies: {tree_path}"
        )

    inserted_edges = sum(degree - 2 for degree in invalid_before)
    if invalid_before:
        resolve_polytomies(tree)

    invalid_after = nonbinary_nodes(tree)
    if invalid_after:
        raise ValueError(
            f"Gene tree is still not strictly binary after resolution "
            f"({invalid_after[:10]}): {tree_path}"
        )
    if leaf_labels(tree) != leaves_before:
        raise ValueError(f"Gene-tree leaf set changed during resolution: {tree_path}")

    if invalid_before:
        temporary = tree_path.with_name(f"{tree_path.name}.resolved.tmp")
        temporary.write_text(to_newick(tree) + ";\n")
        temporary.replace(tree_path)

    return sum(leaves_before.values()), len(invalid_before), inserted_edges


def manifest_tree_rows(manifest):
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"orthogroup", "gene_tree", "mapping"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"SpeciesRax manifest has an invalid header: {manifest}"
            )
        yield from reader


def main():
    args = parse_args()
    manifests = sorted(Path(".").glob(args.manifest_glob))
    if not manifests:
        raise ValueError(
            f"No SpeciesRax manifests matched: {args.manifest_glob}"
        )

    rows = []
    seen_trees = set()
    for manifest in manifests:
        for row in manifest_tree_rows(manifest):
            orthogroup = row["orthogroup"]
            tree_path = Path(row["gene_tree"])
            if tree_path.is_absolute() or ".." in tree_path.parts:
                raise ValueError(
                    f"Unsafe gene-tree path for {orthogroup}: {tree_path}"
                )
            if tree_path in seen_trees:
                raise ValueError(f"Duplicate SpeciesRax gene tree: {tree_path}")
            seen_trees.add(tree_path)
            if not tree_path.is_file():
                raise ValueError(
                    f"Missing SpeciesRax gene tree for {orthogroup}: {tree_path}"
                )

            leaves, multifurcations, inserted_edges = resolve_tree(tree_path)
            rows.append(
                {
                    "orthogroup": orthogroup,
                    "gene_tree": str(tree_path),
                    "leaves": leaves,
                    "multifurcations_resolved": multifurcations,
                    "zero_length_edges_inserted": inserted_edges,
                    "status": "resolved" if multifurcations else "binary",
                }
            )

    with args.report.open("w", newline="") as handle:
        fieldnames = [
            "orthogroup",
            "gene_tree",
            "leaves",
            "multifurcations_resolved",
            "zero_length_edges_inserted",
            "status",
        ]
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    changed = sum(row["status"] == "resolved" for row in rows)
    multifurcations = sum(row["multifurcations_resolved"] for row in rows)
    inserted_edges = sum(row["zero_length_edges_inserted"] for row in rows)
    print(
        "SpeciesRax gene-tree preparation: "
        f"{len(rows)} checked; {changed} trees resolved; "
        f"{multifurcations} multifurcations resolved; "
        f"{inserted_edges} zero-length edges inserted"
    )


if __name__ == "__main__":
    main()
