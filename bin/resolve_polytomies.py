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
"""

import sys


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


def main():
    if len(sys.argv) != 3:
        print("Usage: resolve_polytomies.py <input.newick> <output.newick>",
              file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        newick = f.read().strip()

    tree = parse_newick(newick)
    resolve_polytomies(tree)

    with open(sys.argv[2], "w") as f:
        f.write(to_newick(tree) + ";\n")


if __name__ == "__main__":
    main()
