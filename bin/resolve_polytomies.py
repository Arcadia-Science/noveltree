#!/usr/bin/env python3
"""Resolve polytomies in a Newick tree to make it strictly binary.

Matches the behavior of ape::multi2di(tree, random=FALSE):
  - Children are taken in parse order (left-to-right in Newick).
  - First child peels off at each level, forming a right-ladder.
  - New internal edges get branch length 0.
  - Original branch lengths and node labels are preserved.
  - Root polytomies are treated the same as internal ones.

Usage:
    resolve_polytomies.py <input.newick> <output.newick>
        [--min-branch-length FLOAT --max-branch-length FLOAT]
        [--branch-length-qc FILE]
"""

import argparse
import math


def parse_newick(s):
    """Minimal recursive-descent Newick parser. Returns nested dicts."""
    pos = [0]
    s = s.strip().rstrip(";")

    def skip_ws():
        while pos[0] < len(s) and s[pos[0]] in " \t\n\r":
            pos[0] += 1

    def parse_node():
        skip_ws()
        node = {"children": [], "label": "", "brlen": None}

        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            node["children"].append(parse_node())
            while pos[0] < len(s) and s[pos[0]] == ",":
                pos[0] += 1
                node["children"].append(parse_node())
            if pos[0] < len(s) and s[pos[0]] == ")":
                pos[0] += 1

        # Label (may include NHX comments in brackets)
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ":,[);":
            pos[0] += 1
        node["label"] = s[start:pos[0]].strip()

        # Branch length
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            start = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",[);":
                pos[0] += 1
            node["brlen"] = s[start:pos[0]].strip()

        skip_ws()
        return node

    return parse_node()


def resolve_polytomies(node):
    """Recursively resolve polytomies using right-ladder (ape::multi2di).

    For a node with children [A, B, C, D], produces:
        (A, (B, (C, D):0):0)
    First child peels off at each level; last two share the deepest node.
    """
    # Resolve children first (bottom-up)
    for child in node["children"]:
        resolve_polytomies(child)

    # While this node has more than 2 children, peel off the first child
    # and push the last two into a new internal node (right-ladder, bottom-up)
    while len(node["children"]) > 2:
        # Group the last two children under a new internal node
        c_last = node["children"].pop()
        c_prev = node["children"].pop()
        new_node = {"children": [c_prev, c_last], "label": "", "brlen": "0"}
        node["children"].append(new_node)


def normalize_branch_lengths(node, minimum, maximum):
    """Clamp every non-root branch to likelihood-library bounds.

    GeneRax converts non-positive branch lengths to 0.1 internally but leaves
    positive values below libpll's 1e-6 optimization bound untouched. Applying
    one consistent bound after polytomy resolution avoids both start states.
    The root's Newick length, if present, is not an edge and is left unchanged.
    """
    stats = {
        "branches_total": 0,
        "branches_adjusted": 0,
        "branches_missing": 0,
        "branches_nonfinite_or_invalid": 0,
        "branches_negative": 0,
        "branches_zero": 0,
        "branches_positive_below_minimum": 0,
        "branches_above_maximum": 0,
    }

    def bounded_text(value):
        return format(value, ".12g")

    def visit(parent):
        for child in parent["children"]:
            stats["branches_total"] += 1
            raw = child["brlen"]
            replacement = None

            if raw is None or not raw.strip():
                stats["branches_missing"] += 1
                replacement = minimum
            else:
                try:
                    value = float(raw)
                except ValueError:
                    stats["branches_nonfinite_or_invalid"] += 1
                    replacement = minimum
                else:
                    if not math.isfinite(value):
                        stats["branches_nonfinite_or_invalid"] += 1
                        replacement = minimum
                    elif value < 0:
                        stats["branches_negative"] += 1
                        replacement = minimum
                    elif value == 0:
                        stats["branches_zero"] += 1
                        replacement = minimum
                    elif value < minimum:
                        stats["branches_positive_below_minimum"] += 1
                        replacement = minimum
                    elif value > maximum:
                        stats["branches_above_maximum"] += 1
                        replacement = maximum

            if replacement is not None:
                child["brlen"] = bounded_text(replacement)
                stats["branches_adjusted"] += 1

            visit(child)

    visit(node)
    return stats


def to_newick(node):
    """Convert tree dict back to Newick string."""
    parts = []
    if node["children"]:
        child_strs = ",".join(to_newick(c) for c in node["children"])
        parts.append(f"({child_strs})")
    parts.append(node["label"])
    if node["brlen"] is not None:
        parts.append(f":{node['brlen']}")
    return "".join(parts)


def write_branch_length_qc(path, stats, minimum, maximum):
    with open(path, "w") as handle:
        handle.write("metric\tvalue\n")
        for metric, value in stats.items():
            handle.write(f"{metric}\t{value}\n")
        handle.write(f"minimum_allowed\t{format(minimum, '.12g')}\n")
        handle.write(f"maximum_allowed\t{format(maximum, '.12g')}\n")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_newick")
    parser.add_argument("output_newick")
    parser.add_argument("--min-branch-length", type=float)
    parser.add_argument("--max-branch-length", type=float)
    parser.add_argument("--branch-length-qc")
    args = parser.parse_args()

    bounds = (args.min_branch_length, args.max_branch_length)
    if (bounds[0] is None) != (bounds[1] is None):
        parser.error(
            "--min-branch-length and --max-branch-length must be supplied together"
        )
    if bounds[0] is not None:
        if not all(math.isfinite(value) for value in bounds):
            parser.error("branch-length bounds must be finite")
        if bounds[0] <= 0:
            parser.error("--min-branch-length must be positive")
        if bounds[1] < bounds[0]:
            parser.error(
                "--max-branch-length must be greater than or equal to the minimum"
            )
    if args.branch_length_qc and bounds[0] is None:
        parser.error("--branch-length-qc requires branch-length bounds")
    return args


def main():
    args = parse_args()

    with open(args.input_newick) as f:
        newick = f.read().strip()

    tree = parse_newick(newick)
    resolve_polytomies(tree)

    if args.min_branch_length is not None:
        stats = normalize_branch_lengths(
            tree, args.min_branch_length, args.max_branch_length
        )
        print(
            "GeneRax branch-length normalization: "
            f"{stats['branches_total']} checked; "
            f"{stats['branches_adjusted']} adjusted "
            f"(missing={stats['branches_missing']}, "
            f"nonfinite_or_invalid={stats['branches_nonfinite_or_invalid']}, "
            f"negative={stats['branches_negative']}, "
            f"zero={stats['branches_zero']}, "
            "positive_below_minimum="
            f"{stats['branches_positive_below_minimum']}, "
            f"above_maximum={stats['branches_above_maximum']})"
        )
        if args.branch_length_qc:
            write_branch_length_qc(
                args.branch_length_qc,
                stats,
                args.min_branch_length,
                args.max_branch_length,
            )

    with open(args.output_newick, "w") as f:
        f.write(to_newick(tree) + ";\n")


if __name__ == "__main__":
    main()
