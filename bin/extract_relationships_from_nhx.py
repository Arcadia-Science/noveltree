#!/usr/bin/env python3
"""Extract ortholog, paralog, xenolog pairs and HOG membership from GeneRax NHX.

Replaces the orthoxml-tools pipeline (from-nhx + export-pairs + extract_hog_membership)
with a single lightweight script that parses the NHX tree directly. No XML DOM,
O(n) memory for the tree structure, pairs streamed to disk.

GeneRax NHX format:
  [&&NHX:S=species:D=Y/N:H=Y/N:B=brlen]

Event classification at each internal node:
  - D=N, H=N  → Speciation: cross-child pairs are orthologs
  - D=Y       → Duplication: cross-child pairs are paralogs
  - H=Y       → Transfer: cross-child pairs are xenologs

HOG membership (hierarchical):
  Each speciation node defines a HOG keyed by its species tree node (S= tag).
  A gene belongs to HOGs at ALL speciation ancestors, from root to its deepest
  speciation node. HOG IDs use a short numeric index (N0, N1, ...) derived from
  a deterministic mapping built from the species tree. A lookup table mapping
  these indices to the full species tree node names is also emitted.

Usage:
  extract_relationships_from_nhx.py <nhx_file> <species_tree> <og_name> <outprefix>

Outputs:
  <outprefix>_orthologs.tsv       (gene1, gene2, og)
  <outprefix>_paralogs.tsv        (gene1, gene2, og)
  <outprefix>_xenologs.tsv        (gene1, gene2, species1, species2, og)
  <outprefix>_hog_membership.tsv  (protein_id, species, hog_id, og)
  spp_tree_node_lookup.tsv        (hog_id, spp_tree_node — same for all OGs)
"""

import re
import sys


# ---------------------------------------------------------------------------
# Newick / NHX parser
# ---------------------------------------------------------------------------
def parse_nhx_tree(nhx_string):
    """Minimal recursive-descent NHX parser.

    Returns a nested dict tree:
        {"children": [...], "name": str, "nhx": {key: val, ...}}
    """
    pos = [0]

    def skip_ws():
        while pos[0] < len(nhx_string) and nhx_string[pos[0]] in " \t\n\r":
            pos[0] += 1

    def parse_node():
        skip_ws()
        node = {"children": [], "name": "", "nhx": {}}

        if pos[0] < len(nhx_string) and nhx_string[pos[0]] == "(":
            pos[0] += 1
            node["children"].append(parse_node())
            while pos[0] < len(nhx_string) and nhx_string[pos[0]] == ",":
                pos[0] += 1
                node["children"].append(parse_node())
            if pos[0] < len(nhx_string) and nhx_string[pos[0]] == ")":
                pos[0] += 1

        start = pos[0]
        while pos[0] < len(nhx_string) and nhx_string[pos[0]] not in ":,[);":
            pos[0] += 1
        node["name"] = nhx_string[start : pos[0]].strip()

        if pos[0] < len(nhx_string) and nhx_string[pos[0]] == ":":
            pos[0] += 1
            while pos[0] < len(nhx_string) and nhx_string[pos[0]] not in "[,);":
                pos[0] += 1

        if pos[0] < len(nhx_string) and nhx_string[pos[0]] == "[":
            depth = 1
            bracket_start = pos[0]
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

    return parse_node()


def parse_newick(newick_string):
    """Parse a plain Newick tree (no NHX) into nested dicts with node labels."""
    pos = [0]
    s = newick_string.strip().rstrip(";")

    def skip_ws():
        while pos[0] < len(s) and s[pos[0]] in " \t\n\r":
            pos[0] += 1

    def parse_node():
        skip_ws()
        node = {"children": [], "name": ""}

        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            node["children"].append(parse_node())
            while pos[0] < len(s) and s[pos[0]] == ",":
                pos[0] += 1
                node["children"].append(parse_node())
            if pos[0] < len(s) and s[pos[0]] == ")":
                pos[0] += 1

        # Label
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ":,[);":
            pos[0] += 1
        node["name"] = s[start : pos[0]].strip()

        # Branch length
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            while pos[0] < len(s) and s[pos[0]] not in ",[);":
                pos[0] += 1

        skip_ws()
        return node

    return parse_node()


def species_from_tip(tip_label):
    """Extract species name by removing the protein ID after the last underscore."""
    return re.sub(r"_[^_]+$", "", tip_label)


# ---------------------------------------------------------------------------
# Species tree node index — deterministic mapping
# ---------------------------------------------------------------------------
def build_spp_tree_node_index(spp_tree_file):
    """Read the species tree and build an internal node name → short ID mapping.

    Follows OrthoFinder convention: N0 = root, then depth-first traversal order.
    Only internal nodes get N-indices; tip species keep their names.

    Returns: dict {node_name: "N0", ...}
    """
    with open(spp_tree_file) as f:
        newick = f.read().strip()

    tree = parse_newick(newick)

    # Depth-first traversal: root = N0, then children left-to-right
    index = {}
    counter = [0]

    def dfs_label(node):
        if node["children"]:
            # Internal node — assign next N-index
            if node["name"]:
                index[node["name"]] = f"N{counter[0]}"
                counter[0] += 1
            for child in node["children"]:
                dfs_label(child)

    dfs_label(tree)
    return index


# ---------------------------------------------------------------------------
# Leaf collection (cached per subtree to avoid recomputation)
# ---------------------------------------------------------------------------
def cache_leaves(node):
    """Recursively cache the leaf list for every node (bottom-up)."""
    if not node["children"]:
        node["_leaves"] = [node["name"]]
    else:
        leaves = []
        for child in node["children"]:
            cache_leaves(child)
            leaves.extend(child["_leaves"])
        node["_leaves"] = leaves


# ---------------------------------------------------------------------------
# Pair extraction — single traversal, streaming output
# ---------------------------------------------------------------------------
def extract_pairs(node, og, ortho_fh, para_fh, xeno_fh):
    """Traverse tree; at each internal node emit cross-child pairs."""
    if not node["children"]:
        return

    # Determine event type
    is_dup = node["nhx"].get("D", "N").startswith("Y")
    is_transfer = node["nhx"].get("H", "N").startswith("Y")

    # Emit pairs between all child combinations (handles polytomies)
    children = node["children"]
    if len(children) >= 2:
        for i in range(len(children)):
            for j in range(i + 1, len(children)):
                left_leaves = children[i]["_leaves"]
                right_leaves = children[j]["_leaves"]

                if is_transfer:
                    # Xenolog pairs
                    for g1 in left_leaves:
                        sp1 = species_from_tip(g1)
                        for g2 in right_leaves:
                            sp2 = species_from_tip(g2)
                            xeno_fh.write(f"{g1}\t{g2}\t{sp1}\t{sp2}\t{og}\n")
                elif is_dup:
                    # Paralog pairs
                    for g1 in left_leaves:
                        for g2 in right_leaves:
                            para_fh.write(f"{g1}\t{g2}\t{og}\n")
                else:
                    # Speciation: cross-child pairs are orthologs (including
                    # same-species co-orthologs from post-speciation duplications)
                    for g1 in left_leaves:
                        for g2 in right_leaves:
                            ortho_fh.write(f"{g1}\t{g2}\t{og}\n")

    # Recurse into children
    for child in children:
        extract_pairs(child, og, ortho_fh, para_fh, xeno_fh)


# ---------------------------------------------------------------------------
# HOG membership — hierarchical, keyed by species tree node
# ---------------------------------------------------------------------------
def extract_hog_membership(node, og, node_index, rows=None, ancestor_hogs=None):
    """Assign each gene to ALL speciation-defined HOGs in its ancestry.

    Each speciation node defines a HOG named by its species tree node (S= tag),
    mapped to a short index via node_index. A gene gets one row per speciation
    ancestor, from the root down to its deepest speciation node.

    Duplication/transfer nodes do NOT define new HOGs — they inherit.

    Returns a list of (protein_id, species, hog_id, og) tuples.
    """
    if rows is None:
        rows = []
    if ancestor_hogs is None:
        ancestor_hogs = []

    is_dup = node["nhx"].get("D", "N").startswith("Y")
    is_transfer = node["nhx"].get("H", "N").startswith("Y")

    # Speciation node → new HOG level
    if node["children"] and not is_dup and not is_transfer:
        spp_tree_node = node["nhx"].get("S", "unknown")
        hog_id = node_index.get(spp_tree_node, spp_tree_node)
        ancestor_hogs = ancestor_hogs + [hog_id]

    if not node["children"]:
        # Leaf node — collect one row per ancestral HOG level
        species = species_from_tip(node["name"])
        for hog_id in ancestor_hogs:
            rows.append((node["name"], species, hog_id, og))
        return rows

    for child in node["children"]:
        extract_hog_membership(child, og, node_index, rows, ancestor_hogs)

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) != 5:
        print(
            "Usage: extract_relationships_from_nhx.py <nhx_file> <species_tree> <og_name> <outprefix>",
            file=sys.stderr,
        )
        sys.exit(1)

    nhx_file = sys.argv[1]
    spp_tree_file = sys.argv[2]
    og = sys.argv[3]
    prefix = sys.argv[4]

    # Increase recursion limit for deep trees
    sys.setrecursionlimit(100_000)

    # Build deterministic species tree node index
    node_index = build_spp_tree_node_index(spp_tree_file)

    # Write the lookup table (identical for every OG — will be deduplicated downstream)
    with open("spp_tree_node_lookup.tsv", "w") as lut_fh:
        lut_fh.write("hog_id\tspp_tree_node\n")
        for name, short_id in sorted(node_index.items(), key=lambda x: x[1]):
            lut_fh.write(f"{short_id}\t{name}\n")

    # Parse the gene family NHX tree
    with open(nhx_file) as f:
        nhx_string = f.read().strip()

    tree = parse_nhx_tree(nhx_string)

    # Cache leaves bottom-up (O(n) total memory for leaf lists)
    cache_leaves(tree)

    # Extract pairs and HOG membership — streaming to files
    ortho_fh = open(f"{prefix}_orthologs.tsv", "w")
    para_fh = open(f"{prefix}_paralogs.tsv", "w")
    xeno_fh = open(f"{prefix}_xenologs.tsv", "w")
    hog_fh = open(f"{prefix}_hog_membership.tsv", "w")
    try:
        # Headers
        ortho_fh.write("gene1\tgene2\tog\n")
        para_fh.write("gene1\tgene2\tog\n")
        xeno_fh.write("gene1\tgene2\tspecies1\tspecies2\tog\n")
        hog_fh.write("protein_id\tspecies\thog_id\tog\n")

        # Pairs
        extract_pairs(tree, og, ortho_fh, para_fh, xeno_fh)

        # HOG membership (hierarchical, sorted by HOG ID)
        hog_rows = extract_hog_membership(tree, og, node_index)
        hog_rows.sort(key=lambda r: (r[2], r[1], r[0]))  # sort by hog_id, species, protein
        for protein_id, species, hog_id, og_name in hog_rows:
            hog_fh.write(f"{protein_id}\t{species}\t{hog_id}\t{og_name}\n")
    finally:
        ortho_fh.close()
        para_fh.close()
        xeno_fh.close()
        hog_fh.close()


if __name__ == "__main__":
    main()
