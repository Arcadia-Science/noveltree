#!/usr/bin/env python3
"""Write one compact metadata row for every retained orthogroup FASTA."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-tree-dir", required=True, type=Path)
    parser.add_argument("--gene-tree-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def fasta_stats(path):
    n_sequences = 0
    maximum_length = 0
    current_length = 0
    species = set()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if n_sequences:
                    maximum_length = max(maximum_length, current_length)
                identifier = line[1:].strip().split(maxsplit=1)[0]
                if not identifier:
                    raise ValueError(f"Empty FASTA identifier in {path}:{line_number}")
                n_sequences += 1
                current_length = 0
                species.add(identifier.rsplit("_", 1)[0])
            else:
                if not n_sequences and line.strip():
                    raise ValueError(f"Sequence data before first header in {path}:{line_number}")
                current_length += len(line.strip())

    if n_sequences:
        maximum_length = max(maximum_length, current_length)
    else:
        raise ValueError(f"Retained orthogroup FASTA is empty: {path}")

    return n_sequences, maximum_length, len(species)


def rows_for(directory, family_set):
    for path in sorted(directory.glob("*.fa")):
        n_sequences, maximum_length, n_species = fasta_stats(path)
        yield {
            "orthogroup": path.stem,
            "family_set": family_set,
            "n_seq": n_sequences,
            "max_len": maximum_length,
            "n_species": n_species,
        }


def main():
    args = parse_args()
    fieldnames = ["orthogroup", "family_set", "n_seq", "max_len", "n_species"]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_for(args.species_tree_dir, "species_tree"))
        writer.writerows(rows_for(args.gene_tree_dir, "gene_tree"))


if __name__ == "__main__":
    main()
