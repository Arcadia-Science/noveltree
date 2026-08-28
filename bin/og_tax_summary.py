#!/usr/bin/env python3
"""
Summarize orthogroup distribution across species and filter orthogroups
for species tree vs gene tree inference.

The gene-count table is processed as a stream so memory use stays constant
for datasets containing hundreds of thousands or millions of orthogroups.
"""

import csv
import math
import sys
from collections import Counter


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
    if len(sys.argv) != 11:
        print(
            f"Usage: {sys.argv[0]} <og_gene_counts.tsv> <samplesheet.csv> "
            "<min_num_seqs> <min_num_spp> <min_prop_spp_for_spptree> "
            "<max_copy_num> <speciesrax_min_occupancy> "
            "<speciesrax_max_mean_copies> <speciesrax_max_copies_per_species> "
            "<speciesrax_max_total_leaves_factor>",
            file=sys.stderr,
        )
        sys.exit(1)

    og_counts_path = sys.argv[1]
    samplesheet_path = sys.argv[2]
    num_seq_filt = int(sys.argv[3])
    num_spp_filt = int(sys.argv[4])
    prop_spp_spptree_filt = float(sys.argv[5])
    copy_num_filt = float(sys.argv[6])
    speciesrax_min_occupancy = float(sys.argv[7])
    speciesrax_max_mean_copies = float(sys.argv[8])
    speciesrax_max_copies_per_species = int(sys.argv[9])
    speciesrax_max_total_leaves_factor = float(sys.argv[10])
    if not 0 < speciesrax_min_occupancy <= 1:
        raise ValueError("SpeciesRax minimum occupancy must be in (0, 1]")
    if speciesrax_max_mean_copies < 1:
        raise ValueError("SpeciesRax maximum mean copies must be at least 1")
    if speciesrax_max_copies_per_species < 1:
        raise ValueError("SpeciesRax maximum copies per species must be at least 1")
    if speciesrax_max_total_leaves_factor < 1:
        raise ValueError("SpeciesRax maximum total-leaves factor must be at least 1")
    maximize_csv_field_size()

    with open(samplesheet_path, newline="") as handle:
        num_species = sum(1 for _row in csv.DictReader(handle))
    num_spp_spptree_filt = math.ceil(num_species * prop_spp_spptree_filt)
    speciesrax_min_species = math.ceil(num_species * speciesrax_min_occupancy)
    speciesrax_max_total_leaves = math.floor(
        num_species * speciesrax_max_total_leaves_factor
    )

    fieldnames = ["orthogroup", "num_spp", "total_copy_num", "mean_copy_num"]
    total_ogs = 0
    total_spptree = 0
    total_genetree = 0
    family_report_rows = []
    selected_species_families = Counter()
    selected_species_copies = Counter()
    candidate_species_families = Counter()

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

            is_base_spptree_candidate = (
                total_copy_num >= num_seq_filt
                and mean_copy_num <= copy_num_filt
                and num_spp >= num_spp_filt
                and num_spp >= num_spp_spptree_filt
            )
            reasons = []
            if num_spp < num_spp_spptree_filt:
                reasons.append("base_species_occupancy_below_minimum")
            if mean_copy_num > copy_num_filt:
                reasons.append("base_mean_copies_above_maximum")
            if num_spp < speciesrax_min_species:
                reasons.append("species_occupancy_below_minimum")
            if mean_copy_num > speciesrax_max_mean_copies:
                reasons.append("mean_copies_above_maximum")
            max_copies = max(present_counts, default=0)
            if max_copies > speciesrax_max_copies_per_species:
                reasons.append("species_copy_count_above_maximum")
            if total_copy_num > speciesrax_max_total_leaves:
                reasons.append("total_leaves_above_maximum")
            is_spptree = is_base_spptree_candidate and not reasons
            is_genetree = (
                total_copy_num >= num_seq_filt
                and num_spp >= num_spp_filt
                and not is_spptree
            )

            if total_copy_num >= num_seq_filt and num_spp >= num_spp_filt:
                for column_index, count in zip(species_columns, counts):
                    if count > 0:
                        species = clean_header[column_index]
                        candidate_species_families[species] += 1
                        if is_spptree:
                            selected_species_families[species] += 1
                            selected_species_copies[species] += count
                family_report_rows.append(
                    {
                        **result,
                        "total_species": num_species,
                        "species_occupancy": f"{num_spp / num_species:.6f}",
                        "max_copies_any_species": max_copies,
                        "max_total_leaves": speciesrax_max_total_leaves,
                        "selected": str(is_spptree).lower(),
                        "exclusion_reasons": ";".join(reasons),
                    }
                )

            all_writer.writerow(result)
            total_ogs += 1
            if is_spptree:
                spptree_writer.writerow(result)
                total_spptree += 1
            elif is_genetree:
                genetree_writer.writerow(result)
                total_genetree += 1

    report_fields = [
        "orthogroup",
        "num_spp",
        "total_copy_num",
        "mean_copy_num",
        "total_species",
        "species_occupancy",
        "max_copies_any_species",
        "max_total_leaves",
        "selected",
        "exclusion_reasons",
    ]
    with open("speciesrax_family_selection.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=report_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(family_report_rows)

    with open(
        "speciesrax_selected_species_coverage.tsv", "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "species",
                "input_families",
                "selected_families",
                "selected_copies",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for species in clean_header[1:]:
            if species == "Total":
                continue
            writer.writerow(
                {
                    "species": species,
                    "input_families": candidate_species_families[species],
                    "selected_families": selected_species_families[species],
                    "selected_copies": selected_species_copies[species],
                }
            )

    if not total_spptree:
        raise ValueError("Upstream SpeciesRax filtering retained no gene families")
    missing_species = sorted(
        species
        for species in clean_header[1:]
        if species != "Total" and not selected_species_families[species]
    )
    if missing_species:
        raise ValueError(
            "Upstream SpeciesRax filtering removed all families for "
            f"{len(missing_species)} species: " + ", ".join(missing_species[:10])
        )

    print(f"Total OGs: {total_ogs}", file=sys.stderr)
    print(f"Species tree core OGs: {total_spptree}", file=sys.stderr)
    print(f"Gene tree core OGs: {total_genetree}", file=sys.stderr)


if __name__ == "__main__":
    main()
