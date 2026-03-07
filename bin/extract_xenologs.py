#!/usr/bin/env python3
"""Extract xenolog (HGT) pairs from GeneRax reconciliated NHX files.

For each internal node annotated with H=Y (horizontal transfer), the two
child clades are xenologous to each other. This script emits all pairwise
combinations of leaves between the two child clades.

NHX format produced by GeneRax:
  [&&NHX:S=species:D=Y/N:H=Y/N:B=brlen]
"""

import re
import sys
from itertools import product


def parse_nhx_tree(nhx_string):
    """Minimal recursive-descent NHX parser.

    Returns a nested dict tree:
        {"children": [...], "name": str, "nhx": {key: val, ...}}
    Leaf nodes have no "children" key (or empty list).
    """
    pos = [0]  # mutable index

    def skip_ws():
        while pos[0] < len(nhx_string) and nhx_string[pos[0]] in " \t\n\r":
            pos[0] += 1

    def parse_node():
        skip_ws()
        node = {"children": [], "name": "", "nhx": {}}

        # Children
        if pos[0] < len(nhx_string) and nhx_string[pos[0]] == "(":
            pos[0] += 1  # skip '('
            node["children"].append(parse_node())
            while pos[0] < len(nhx_string) and nhx_string[pos[0]] == ",":
                pos[0] += 1  # skip ','
                node["children"].append(parse_node())
            if pos[0] < len(nhx_string) and nhx_string[pos[0]] == ")":
                pos[0] += 1  # skip ')'

        # Label (until ':', '[', ',', ')', ';')
        start = pos[0]
        while pos[0] < len(nhx_string) and nhx_string[pos[0]] not in ":,[);":
            pos[0] += 1
        node["name"] = nhx_string[start : pos[0]].strip()

        # Branch length
        if pos[0] < len(nhx_string) and nhx_string[pos[0]] == ":":
            pos[0] += 1
            bl_start = pos[0]
            while pos[0] < len(nhx_string) and nhx_string[pos[0]] not in "[,);":
                pos[0] += 1
            # branch length ignored for our purposes

        # NHX comment
        if pos[0] < len(nhx_string) and nhx_string[pos[0]] == "[":
            bracket_start = pos[0]
            depth = 1
            pos[0] += 1
            while pos[0] < len(nhx_string) and depth > 0:
                if nhx_string[pos[0]] == "[":
                    depth += 1
                elif nhx_string[pos[0]] == "]":
                    depth -= 1
                pos[0] += 1
            comment = nhx_string[bracket_start + 1 : pos[0] - 1]
            if comment.startswith("&&NHX:"):
                for pair in comment[6:].split(":"):
                    if "=" in pair:
                        k, v = pair.split("=", 1)
                        node["nhx"][k] = v

        skip_ws()
        return node

    tree = parse_node()
    return tree


def get_leaves(node):
    """Return list of leaf names under this node."""
    if not node["children"]:
        return [node["name"]]
    leaves = []
    for child in node["children"]:
        leaves.extend(get_leaves(child))
    return leaves


def species_from_tip(tip_label):
    """Extract species name by removing the protein ID after the last underscore."""
    return re.sub(r"_[^_]+$", "", tip_label)


def find_transfer_nodes(node):
    """Yield nodes where H=Y (horizontal transfer)."""
    h_val = node["nhx"].get("H", "N")
    if h_val.startswith("Y"):
        yield node
    for child in node["children"]:
        yield from find_transfer_nodes(child)


def main():
    if len(sys.argv) != 3:
        print("Usage: extract_xenologs.py <nhx_file> <og_name>", file=sys.stderr)
        sys.exit(1)

    nhx_file = sys.argv[1]
    og = sys.argv[2]

    with open(nhx_file) as f:
        nhx_string = f.read().strip()

    tree = parse_nhx_tree(nhx_string)

    # Header
    print("gene1\tgene2\tspecies1\tspecies2\tog")

    for tnode in find_transfer_nodes(tree):
        if len(tnode["children"]) < 2:
            continue
        # The two child clades are xenologous to each other
        left_leaves = get_leaves(tnode["children"][0])
        right_leaves = get_leaves(tnode["children"][1])
        for g1, g2 in product(left_leaves, right_leaves):
            sp1 = species_from_tip(g1)
            sp2 = species_from_tip(g2)
            print(f"{g1}\t{g2}\t{sp1}\t{sp2}\t{og}")


if __name__ == "__main__":
    main()
