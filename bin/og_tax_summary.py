#!/usr/bin/env python3
"""
Summarize orthogroup distribution across species and filter orthogroups
for species tree vs gene tree inference.
"""

import csv
import math
import sys


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
    copy_num_filt1 = float(sys.argv[6])

    # Read samplesheet to count species
    num_species = 0
    with open(samplesheet_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            num_species += 1

    num_spp_spptree_filt = round(num_species * prop_spp_spptree_filt)

    # Read orthogroup gene counts
    with open(og_counts_path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        # Clean column names (strip everything after first dot, matching R's gsub)
        clean_header = []
        for col in header:
            dot_pos = col.find(".")
            clean_header.append(col[:dot_pos] if dot_pos >= 0 else col)

        # Find column indices
        og_col = 0  # "Orthogroup"
        total_col = clean_header.index("Total")
        spp_cols = [i for i in range(1, len(clean_header)) if i != total_col]

        results = []
        for row in reader:
            og_name = row[og_col]
            total_copy_num = int(row[total_col])

            # Per-species counts
            counts = {}
            for i in spp_cols:
                spp_name = clean_header[i]
                counts[spp_name] = int(row[i])

            # Number of species present (count > 0)
            num_spp = sum(1 for c in counts.values() if c > 0)

            # Mean copy number across species that are present
            present_counts = [c for c in counts.values() if c > 0]
            mean_copy_num = sum(present_counts) / len(present_counts) if present_counts else 0.0

            results.append({
                "orthogroup": og_name,
                "num_spp": num_spp,
                "total_copy_num": total_copy_num,
                "mean_copy_num": mean_copy_num,
            })

    # Filter: species tree core OGs
    spptree_core = [
        r for r in results
        if r["total_copy_num"] >= num_seq_filt
        and r["mean_copy_num"] <= copy_num_filt1
        and r["num_spp"] >= num_spp_filt
        and r["num_spp"] >= num_spp_spptree_filt
    ]
    spptree_og_names = {r["orthogroup"] for r in spptree_core}

    # Filter: gene tree core OGs (excludes species tree core to avoid redundancy)
    genetree_core = [
        r for r in results
        if r["total_copy_num"] >= num_seq_filt
        and r["num_spp"] >= num_spp_filt
        and r["orthogroup"] not in spptree_og_names
    ]

    # Write outputs
    fieldnames = ["orthogroup", "num_spp", "total_copy_num", "mean_copy_num"]

    def write_csv(filename, data):
        with open(filename, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in data:
                # Format mean_copy_num to match R's default CSV output
                row_out = dict(row)
                row_out["mean_copy_num"] = f"{row['mean_copy_num']:.6f}" if row["mean_copy_num"] != int(row["mean_copy_num"]) else str(int(row["mean_copy_num"]))
                writer.writerow(row_out)

    write_csv("all_ogs_counts.csv", results)
    write_csv("spptree_core_ogs_counts.csv", spptree_core)
    write_csv("genetree_core_ogs_counts.csv", genetree_core)

    print(f"Total OGs: {len(results)}", file=sys.stderr)
    print(f"Species tree core OGs: {len(spptree_core)}", file=sys.stderr)
    print(f"Gene tree core OGs: {len(genetree_core)}", file=sys.stderr)


if __name__ == "__main__":
    main()
