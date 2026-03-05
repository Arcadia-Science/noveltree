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
            search_header = record.id
            result = re.search(f"{species}(_tr)?", search_header)
            if result is None:
                print(f"Species {species}, not present in header {search_header}.")
                sys.exit(1)
            delimeter = search_header[result.end()]
            prot_id = search_header.split(" ")[0].split(delimeter)[1]
            output.write(f"{prot_id}\n")


if __name__ == "__main__":
    args = parse_args()
    main(args.species, args.fasta.name, args.output.name)
