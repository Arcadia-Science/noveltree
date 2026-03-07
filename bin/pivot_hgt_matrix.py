#!/usr/bin/env python3
"""Pivot concatenated long-format HGT counts into a species x species matrix."""

import csv
import sys
from collections import defaultdict


def main(input_file, output_file):
    counts = defaultdict(lambda: defaultdict(int))
    species = set()

    with open(input_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            d, r, c = row["donor"], row["recipient"], int(float(row["count"]))
            counts[d][r] += c
            species.add(d)
            species.add(r)

    species_order = sorted(species)

    with open(output_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([""] + species_order)
        for donor in species_order:
            vals = [counts[donor].get(r, 0) for r in species_order]
            writer.writerow([donor] + vals)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: pivot_hgt_matrix.py <input_long.tsv> <output_matrix.tsv>")
    main(sys.argv[1], sys.argv[2])
