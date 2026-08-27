#!/usr/bin/env python3
"""Build GeneRax mappings for retained families from canonical proteome maps."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-map-glob", required=True)
    parser.add_argument("--species-tree-dir", required=True, type=Path)
    parser.add_argument("--gene-tree-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_canonical_maps(pattern):
    paths = sorted(Path(".").glob(pattern))
    if not paths:
        raise ValueError(f"No canonical protein mappings matched {pattern!r}")
    lookup = {}
    for path in paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"species", "canonical_protein_id"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(f"Invalid canonical mapping header: {path}")
            for row in reader:
                protein = row["canonical_protein_id"]
                species = row["species"]
                if protein in lookup:
                    raise ValueError(f"Duplicate canonical protein ID across mappings: {protein}")
                lookup[protein] = species
    return lookup


def fasta_ids(path):
    identifiers = []
    seen = set()
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            identifier = line[1:].strip().split(maxsplit=1)[0]
            if not identifier:
                raise ValueError(f"Empty FASTA ID at {path}:{line_number}")
            if identifier in seen:
                raise ValueError(f"Duplicate FASTA ID in {path}: {identifier}")
            seen.add(identifier)
            identifiers.append(identifier)
    if not identifiers:
        raise ValueError(f"Empty retained family FASTA: {path}")
    return identifiers


def write_family_maps(lookup, directories, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    seen_families = set()
    for directory in directories:
        for fasta in sorted(directory.glob("*.fa")):
            family = fasta.stem
            if family in seen_families:
                raise ValueError(f"Retained family appears in multiple routes: {family}")
            seen_families.add(family)
            rows = []
            for protein in fasta_ids(fasta):
                try:
                    species = lookup[protein]
                except KeyError as error:
                    raise ValueError(
                        f"Retained protein {protein!r} in {fasta} is absent from canonical mappings"
                    ) from error
                rows.append((protein, species))
            with (output_dir / f"{family}_map.link").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerows(rows)


def main():
    args = parse_args()
    lookup = read_canonical_maps(args.canonical_map_glob)
    write_family_maps(
        lookup,
        [args.species_tree_dir, args.gene_tree_dir],
        args.output_dir,
    )


if __name__ == "__main__":
    main()
