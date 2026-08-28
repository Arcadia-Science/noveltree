#!/usr/bin/env python3
"""Validate and prepare one inferred gene tree for reconciliation."""

import argparse
import csv
from collections import Counter
from pathlib import Path

from resolve_polytomies import (
    normalize_branch_lengths,
    parse_newick,
    resolve_polytomies,
    to_newick,
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroup", required=True)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    parser.add_argument("--min-branch-length", type=float, default=1e-6)
    parser.add_argument("--max-branch-length", type=float, default=100.0)
    return parser.parse_args()


def walk(node):
    yield node
    for child in node["children"]:
        yield from walk(child)


def leaf_labels(node):
    return Counter(
        current["label"] for current in walk(node) if not current["children"]
    )


def fasta_labels(path):
    labels = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                label = line[1:].strip().split(maxsplit=1)[0]
                if not label:
                    raise ValueError(f"Empty FASTA ID at {path}:{line_number}")
                labels.append(label)
    if not labels:
        raise ValueError(f"Empty alignment: {path}")
    return Counter(labels)


def mapping_labels(path):
    labels = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 2:
                raise ValueError(f"Invalid mapping row at {path}:{line_number}")
            protein, species = fields
            if not protein.startswith(f"{species}_"):
                raise ValueError(
                    f"Canonical mapping mismatch at {path}:{line_number}: "
                    f"{protein!r} is not prefixed by {species!r}"
                )
            labels.append(protein)
    if not labels:
        raise ValueError(f"Empty mapping: {path}")
    return Counter(labels)


def prepare(args):
    newick = args.tree.read_text().strip()
    if not newick or not newick.endswith(";"):
        raise ValueError(f"Invalid or unterminated Newick tree: {args.tree}")
    tree = parse_newick(newick)
    tree_leaves = leaf_labels(tree)
    alignment_leaves = fasta_labels(args.alignment)
    mapping_leaves = mapping_labels(args.mapping)
    for name, leaves in (
        ("alignment", alignment_leaves),
        ("mapping", mapping_leaves),
    ):
        if leaves != tree_leaves:
            missing = sorted((tree_leaves - leaves).elements())
            extra = sorted((leaves - tree_leaves).elements())
            raise ValueError(
                f"{args.orthogroup} tree/{name} leaf mismatch: "
                f"missing={missing[:10]}, extra={extra[:10]}"
            )

    nonbinary = [
        len(node["children"])
        for node in walk(tree)
        if node["children"] and len(node["children"]) != 2
    ]
    if any(degree < 2 for degree in nonbinary):
        raise ValueError(f"{args.orthogroup} contains a unary internal node")
    inserted_edges = sum(degree - 2 for degree in nonbinary)
    resolve_polytomies(tree)
    branch_qc = normalize_branch_lengths(
        tree, args.min_branch_length, args.max_branch_length
    )
    if leaf_labels(tree) != tree_leaves:
        raise ValueError(f"{args.orthogroup} leaf set changed during tree preparation")

    args.output.write_text(to_newick(tree) + ";\n")
    row = {
        "orthogroup": args.orthogroup,
        "leaves": sum(tree_leaves.values()),
        "multifurcations_resolved": len(nonbinary),
        "zero_length_edges_inserted": inserted_edges,
        **branch_qc,
        "minimum_allowed": args.min_branch_length,
        "maximum_allowed": args.max_branch_length,
    }
    with args.qc.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(row), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)


def main():
    prepare(parse_args())


if __name__ == "__main__":
    main()
