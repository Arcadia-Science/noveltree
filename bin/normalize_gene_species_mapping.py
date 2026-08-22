#!/usr/bin/env python3
"""Normalize a staged GeneRax gene-to-species mapping against a species tree."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--species-tree", required=True, type=Path)
    parser.add_argument("--orthogroup", default="unknown")
    return parser.parse_args()


def _skip_comment(newick, index):
    depth = 1
    index += 1
    while index < len(newick) and depth:
        if newick[index] == "[":
            depth += 1
        elif newick[index] == "]":
            depth -= 1
        index += 1
    if depth:
        raise ValueError("Unterminated Newick comment")
    return index


def _read_quoted_label(newick, index):
    characters = []
    index += 1
    while index < len(newick):
        character = newick[index]
        if character == "'":
            if index + 1 < len(newick) and newick[index + 1] == "'":
                characters.append("'")
                index += 2
                continue
            return "".join(characters), index + 1
        characters.append(character)
        index += 1
    raise ValueError("Unterminated quoted Newick label")


def read_species_tree_leaves(path):
    """Return leaf labels without requiring a phylogenetics Python package."""
    newick = path.read_text().strip()
    if not newick:
        raise ValueError(f"Empty species tree: {path}")

    leaves = []
    expecting_leaf = True
    index = 0
    delimiters = set("(),:;[]")
    while index < len(newick):
        character = newick[index]
        if character.isspace():
            index += 1
        elif character == "[":
            index = _skip_comment(newick, index)
        elif character == "(":
            expecting_leaf = True
            index += 1
        elif character == ",":
            expecting_leaf = True
            index += 1
        elif character in "):;":
            expecting_leaf = False
            index += 1
        elif character == "'":
            label, index = _read_quoted_label(newick, index)
            if expecting_leaf:
                leaves.append(label)
            expecting_leaf = False
        else:
            end = index
            while (
                end < len(newick)
                and not newick[end].isspace()
                and newick[end] not in delimiters
            ):
                end += 1
            label = newick[index:end]
            if expecting_leaf:
                leaves.append(label)
            expecting_leaf = False
            index = end

    if len(leaves) < 2:
        raise ValueError(f"Species tree contains only {len(leaves)} leaf labels: {path}")
    duplicates = sorted({leaf for leaf in leaves if leaves.count(leaf) > 1})
    if duplicates:
        raise ValueError(
            "Species tree contains duplicate leaf labels: " + ", ".join(duplicates[:10])
        )
    return set(leaves)


def normalize_mapping(mapping_path, species, orthogroup):
    rows = []
    seen_genes = set()
    corrected = 0
    with mapping_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 2:
                raise ValueError(
                    f"Invalid mapping row for {orthogroup} at "
                    f"{mapping_path}:{line_number}; expected two fields"
                )
            gene, mapped_species = fields
            if gene in seen_genes:
                raise ValueError(
                    f"Duplicate gene {gene!r} in mapping for {orthogroup} at "
                    f"{mapping_path}:{line_number}"
                )
            seen_genes.add(gene)

            derived_species, separator, protein = gene.partition("_")
            if not separator or not derived_species or not protein:
                raise ValueError(
                    f"Cannot derive species from gene {gene!r} for {orthogroup} at "
                    f"{mapping_path}:{line_number}"
                )
            if derived_species not in species:
                raise ValueError(
                    f"Gene {gene!r} derives species {derived_species!r}, which is not "
                    f"a leaf in the GeneRax species tree for {orthogroup}"
                )
            corrected += mapped_species != derived_species
            rows.append((gene, derived_species))

    if not rows:
        raise ValueError(f"Empty mapping for {orthogroup}: {mapping_path}")

    # Nextflow stages this as a task-local copy. Atomic replacement repairs the
    # current GeneRax input without modifying the canonical stored mapping.
    temporary = mapping_path.with_name(f"{mapping_path.name}.normalized.tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)
    temporary.replace(mapping_path)
    return len(rows), corrected


def main():
    args = parse_args()
    species = read_species_tree_leaves(args.species_tree)
    rows, corrected = normalize_mapping(args.mapping, species, args.orthogroup)
    print(
        f"GeneRax mapping normalization for {args.orthogroup}: "
        f"{rows} rows validated; {corrected} corrected"
    )


if __name__ == "__main__":
    main()
