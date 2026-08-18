#!/usr/bin/env python3
"""
Summarize orthogroup distribution across species and filter orthogroups
for species tree vs gene tree inference.

The gene-count table is processed as a stream so memory use stays constant
for datasets containing hundreds of thousands or millions of orthogroups.
"""

import csv
import sys


def maximize_csv_field_size():
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def format_mean(value):
    if value == int(value):
        return str(int(value))
    return f"{value:.6f}"


def main():
    if len(sys.argv) != 7:
        print(
            f"Usage: {sys.argv[0]} <og_gene_counts.tsv> <samplesheet.csv> "
            "<min_num_seqs> <min_num_spp> <min_prop_spp_for_spptree> <max_copy_num>",
            file=sys.stderr,
        )
        sys.exit(1)

    og_counts_path = sys.argv[1]
    samplesheet_path = sys.argv[2]
    num_seq_filt = int(sys.argv[3])
    num_spp_filt = int(sys.argv[4])
    prop_spp_spptree_filt = float(sys.argv[5])
    copy_num_filt = float(sys.argv[6])
    maximize_csv_field_size()

    with open(samplesheet_path, newline="") as handle:
        num_species = sum(1 for _row in csv.DictReader(handle))
    num_spp_spptree_filt = round(num_species * prop_spp_spptree_filt)

    fieldnames = ["orthogroup", "num_spp", "total_copy_num", "mean_copy_num"]
    total_ogs = 0
    total_spptree = 0
    total_genetree = 0

    with open(og_counts_path, newline="") as counts_handle, open(
        "all_ogs_counts.csv", "w", newline=""
    ) as all_handle, open(
        "spptree_core_ogs_counts.csv", "w", newline=""
    ) as spptree_handle, open(
        "genetree_core_ogs_counts.csv", "w", newline=""
    ) as genetree_handle:
        reader = csv.reader(counts_handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            raise ValueError(f"Empty orthogroup gene-count table: {og_counts_path}")

        clean_header = [column.split(".", 1)[0] for column in header]
        try:
            total_column = clean_header.index("Total")
        except ValueError as error:
            raise ValueError(
                f"No Total column in orthogroup gene-count table: {og_counts_path}"
            ) from error
        species_columns = [
            index for index in range(1, len(clean_header)) if index != total_column
        ]

        all_writer = csv.DictWriter(all_handle, fieldnames=fieldnames)
        spptree_writer = csv.DictWriter(spptree_handle, fieldnames=fieldnames)
        genetree_writer = csv.DictWriter(genetree_handle, fieldnames=fieldnames)
        for writer in (all_writer, spptree_writer, genetree_writer):
            writer.writeheader()

        for row_number, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"Row {row_number} of {og_counts_path} has {len(row)} fields; "
                    f"expected {len(header)}"
                )
            try:
                total_copy_num = int(row[total_column])
                counts = [int(row[index]) for index in species_columns]
            except ValueError as error:
                raise ValueError(
                    f"Invalid integer count at row {row_number} of {og_counts_path}"
                ) from error

            present_counts = [count for count in counts if count > 0]
            num_spp = len(present_counts)
            mean_copy_num = (
                sum(present_counts) / len(present_counts) if present_counts else 0.0
            )
            result = {
                "orthogroup": row[0],
                "num_spp": num_spp,
                "total_copy_num": total_copy_num,
                "mean_copy_num": format_mean(mean_copy_num),
            }

            is_spptree = (
                total_copy_num >= num_seq_filt
                and mean_copy_num <= copy_num_filt
                and num_spp >= num_spp_filt
                and num_spp >= num_spp_spptree_filt
            )
            is_genetree = (
                total_copy_num >= num_seq_filt
                and num_spp >= num_spp_filt
                and not is_spptree
            )

            all_writer.writerow(result)
            total_ogs += 1
            if is_spptree:
                spptree_writer.writerow(result)
                total_spptree += 1
            elif is_genetree:
                genetree_writer.writerow(result)
                total_genetree += 1

    print(f"Total OGs: {total_ogs}", file=sys.stderr)
    print(f"Species tree core OGs: {total_spptree}", file=sys.stderr)
    print(f"Gene tree core OGs: {total_genetree}", file=sys.stderr)


if __name__ == "__main__":
    main()
