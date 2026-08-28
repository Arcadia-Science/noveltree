#!/usr/bin/env python3
"""Subset an authoritative GeneRax mapping to the records in a FASTA."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_mapping(path):
    mapping = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 2:
                raise ValueError(f"Invalid mapping row at {path}:{line_number}")
            protein, species = fields
            if protein in mapping:
                raise ValueError(f"Duplicate protein in mapping: {protein}")
            mapping[protein] = species
    if not mapping:
        raise ValueError(f"Empty mapping: {path}")
    return mapping


def fasta_ids(path):
    identifiers = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].strip().split(maxsplit=1)[0])
    return identifiers


def subset(mapping_path, fasta_path, output_path):
    mapping = read_mapping(mapping_path)
    identifiers = fasta_ids(fasta_path)
    if not identifiers:
        raise ValueError(f"No sequences in FASTA: {fasta_path}")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate sequence IDs in FASTA: {fasta_path}")
    missing = [identifier for identifier in identifiers if identifier not in mapping]
    if missing:
        raise ValueError(
            f"{len(missing)} FASTA IDs are absent from {mapping_path}: "
            + ", ".join(missing[:10])
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows((identifier, mapping[identifier]) for identifier in identifiers)


def main():
    args = parse_args()
    subset(args.mapping, args.fasta, args.output)


if __name__ == "__main__":
    main()
