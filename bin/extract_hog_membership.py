#!/usr/bin/env python3
"""Extract HOG membership table from an OrthoXML file.

Walks the orthologGroup/paralogGroup hierarchy and assigns each protein
to its deepest (most specific) HOG.

Output: TSV with columns protein_id, species, hog_id, og
"""

import sys
from lxml import etree


# OrthoXML namespace
NS = "http://orthoXML.org/2011/"


def main():
    if len(sys.argv) != 3:
        print("Usage: extract_hog_membership.py <orthoxml_file> <og_name>", file=sys.stderr)
        sys.exit(1)

    orthoxml_file = sys.argv[1]
    og = sys.argv[2]

    tree = etree.parse(orthoxml_file)
    root = tree.getroot()

    # Build gene ID -> (protId, species) lookup from <species><gene> elements
    gene_lookup = {}
    for species_el in root.findall(f"{{{NS}}}species"):
        species_name = species_el.get("name")
        for gene_el in species_el.iter(f"{{{NS}}}gene"):
            gene_id = gene_el.get("id")
            prot_id = gene_el.get("protId", gene_el.get("geneId", gene_id))
            gene_lookup[gene_id] = (prot_id, species_name)

    # Walk orthologGroup/paralogGroup tree, tracking deepest HOG per gene
    # Each gene gets assigned the most specific (deepest) HOG it belongs to
    gene_hog = {}  # gene_id -> hog_id

    def walk_group(group_el, parent_hog=None):
        # Get HOG id from this group element
        hog_id = group_el.get("id", parent_hog)

        # Process direct gene references at this level
        for gene_ref in group_el.findall(f"{{{NS}}}geneRef"):
            gid = gene_ref.get("id")
            # Always overwrite — deeper groups are visited after shallower ones
            gene_hog[gid] = hog_id

        # Recurse into child groups (orthologs or paralogs)
        for child in group_el:
            tag = etree.QName(child).localname
            if tag in ("orthologGroup", "paralogGroup"):
                walk_group(child, parent_hog=hog_id)

    # Find all top-level orthologGroups (HOGs)
    for og_group in root.iter(f"{{{NS}}}orthologGroup"):
        # Only process top-level groups (direct children of <groups>)
        parent = og_group.getparent()
        if parent is not None and etree.QName(parent).localname == "groups":
            walk_group(og_group)

    # Output
    print("protein_id\tspecies\thog_id\tog")
    for gene_id, hog_id in sorted(gene_hog.items(), key=lambda x: x[0]):
        if gene_id in gene_lookup:
            prot_id, species = gene_lookup[gene_id]
            print(f"{prot_id}\t{species}\t{hog_id}\t{og}")
        else:
            print(f"{gene_id}\tunknown\t{hog_id}\t{og}", file=sys.stderr)


if __name__ == "__main__":
    main()
