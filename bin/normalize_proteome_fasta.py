#!/usr/bin/env python3
"""Canonicalize one proteome FASTA and record an authoritative ID mapping."""

import argparse
import csv
import gzip
import json
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    return parser.parse_args()


def open_text(path):
    with path.open("rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return path.open()


def canonical_protein_id(species, raw_identifier):
    prefix = f"{species}_"
    protein = (
        raw_identifier[len(prefix) :]
        if raw_identifier.startswith(prefix)
        else raw_identifier
    )
    # Preserve the accession itself for the two explicitly supported UniProt
    # header conventions while retaining the untouched source ID in the map.
    if protein.startswith(("sp|", "tr|")):
        fields = protein.split("|")
        if len(fields) >= 2 and fields[1]:
            protein = fields[1]
    elif protein.startswith(f"{species}:"):
        protein = protein.split(":", 1)[1]

    # Newick delimiters, whitespace, and underscores inside the protein part
    # are unsafe for the one-underscore canonical species/protein convention.
    protein = re.sub(r"[^A-Za-z0-9.\-]+", "-", protein).strip("-")
    if not protein:
        raise ValueError(f"Empty protein identifier after removing {prefix!r}")
    return f"{species}_{protein}"


def normalize(input_path, species, output_path, mapping_path, qc_path):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.\-]*", species):
        raise ValueError(
            "Canonical species name must contain only letters, digits, dots, "
            f"and hyphens: {species!r}"
        )

    records = []
    seen_original = set()
    seen_canonical = set()
    current = None
    sequence = []
    replacements_u = 0
    replacements_o = 0

    def finish_record():
        nonlocal replacements_u, replacements_o
        if current is None:
            return
        joined = "".join(sequence).replace(" ", "").upper()
        replacements_u += joined.count("U")
        replacements_o += joined.count("O")
        joined = joined.replace("U", "X").replace("O", "X")
        records.append((current[0], current[1], joined))

    with open_text(input_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                finish_record()
                raw = line[1:].strip().split(maxsplit=1)[0]
                if not raw:
                    raise ValueError(f"Empty FASTA identifier at {input_path}:{line_number}")
                canonical = canonical_protein_id(species, raw)
                if raw in seen_original:
                    raise ValueError(f"Duplicate input protein identifier: {raw}")
                if canonical in seen_canonical:
                    raise ValueError(
                        f"Canonical protein identifier collision for {raw!r}: {canonical!r}"
                    )
                seen_original.add(raw)
                seen_canonical.add(canonical)
                current = (raw, canonical)
                sequence = []
            elif line.strip():
                if current is None:
                    raise ValueError(
                        f"Sequence data precedes the first FASTA header at {input_path}:{line_number}"
                    )
                sequence.append("".join(line.split()))
    finish_record()

    if not records:
        raise ValueError(f"No protein records found in {input_path}")

    with output_path.open("w") as output:
        for _original, canonical, residues in records:
            output.write(f">{canonical}\n")
            for start in range(0, len(residues), 80):
                output.write(residues[start : start + 80] + "\n")

    with mapping_path.open("w", newline="") as mapping:
        writer = csv.DictWriter(
            mapping,
            fieldnames=["species", "original_protein_id", "canonical_protein_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for original, canonical, _residues in records:
            writer.writerow(
                {
                    "species": species,
                    "original_protein_id": original,
                    "canonical_protein_id": canonical,
                }
            )

    qc_path.write_text(
        json.dumps(
            {
                "species": species,
                "n_proteins": len(records),
                "selenocysteine_u_replaced_with_x": replacements_u,
                "pyrrolysine_o_replaced_with_x": replacements_o,
                "canonical_id_collisions": 0,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def main():
    args = parse_args()
    normalize(args.input, args.species, args.output, args.mapping, args.qc)


if __name__ == "__main__":
    main()
