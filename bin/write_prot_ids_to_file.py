#!/usr/bin/env python3

import argparse
import re
import sys

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--species",
        type=str,
        required=True
    )
    parser.add_argument(
        "--fasta",
        type=argparse.FileType("r"),
        required=True
    )
    parser.add_argument(
        "--output",
        type=argparse.FileType("w"),
        required=True
    )
    return parser.parse_args()


def main(species, fasta, output):
    with open(output, "w") as output:
        for record in SeqIO.parse(fasta, "fasta"):
            header = record.id
            # Sometimes the species name has an -sp at the end, other times a -.
            # The next line creates a header with only underscore (the standard
            # we require) for searching for the species and delimiter character.
            search_header = header.replace("-", "_")
            result = re.search(f"{species}(_tr)?", search_header)
            if result is None:
                print(f"Species {species}, not present in header {header}.")
                sys.exit(1)
            delimeter = header[result.end()]
            # For the protein id extraction use the header without any character
            # changes.
            prot_id = header.split(" ")[0].split(delimeter)[1]
            output.write(f"{prot_id}\n")


if __name__ == "__main__":
    args = parse_args()
    main(args.species, args.fasta.name, args.output.name)
