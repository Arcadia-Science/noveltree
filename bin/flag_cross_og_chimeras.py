#!/usr/bin/env python3
"""
Detect and remove cross-OG chimeric proteins.

A chimera is a protein whose best DIAMOND hit to a different orthogroup is
>50% of its best hit to its own orthogroup. This signals two distinct
homologous domains fused together.

Reads OrthoFinder's Orthogroups.tsv and the DIAMOND BLAST results already
present in the OrthoFinder working directory.
"""

import argparse
import csv
import glob
import os
import re
import shutil
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Flag cross-OG chimeric proteins")
    parser.add_argument("--orthogroups", required=True, help="Path to Orthogroups.tsv")
    parser.add_argument("--blast_dir", required=True, help="Directory containing DIAMOND BLAST result files")
    parser.add_argument("--og_seqs_dir", required=True, help="Path to Orthogroup_Sequences/ directory")
    parser.add_argument("--report", required=True, help="Output chimera report TSV")
    parser.add_argument(
        "--retained-og-files",
        nargs="+",
        help=(
            "CSV/TSV files containing orthogroups retained by the preliminary "
            "taxon filter. Only proteins assigned to these orthogroups are "
            "evaluated as possible chimeras."
        ),
    )
    parser.add_argument(
        "--gene-counts",
        help="Original OrthoFinder Orthogroups.GeneCount.tsv",
    )
    parser.add_argument(
        "--updated-gene-counts",
        help="Write gene counts after subtracting removed chimeric proteins",
    )
    parser.add_argument(
        "--updated-orthogroups",
        help="Write Orthogroups.tsv after removing chimeric proteins",
    )
    parser.add_argument("--bitscore_min", type=float, default=100, help="Minimum bitscore for a hit (default: 100)")
    parser.add_argument("--evalue_max", type=float, default=1e-10, help="Maximum evalue for a hit (default: 1e-10)")
    parser.add_argument("--ratio_threshold", type=float, default=0.5, help="Min ratio of best non-self OG hit to self OG hit (default: 0.5)")
    return parser.parse_args()


def maximize_csv_field_size():
    """Allow OrthoFinder cells containing arbitrarily large orthogroups."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def read_retained_ogs(paths):
    """Read an orthogroup column from one or more CSV/TSV filter outputs."""
    retained = set()
    if not paths:
        return None

    for path in paths:
        with open(path, newline="") as handle:
            first_line = handle.readline()
            if not first_line:
                continue
            delimiter = "\t" if "\t" in first_line else ","
            handle.seek(0)
            reader = csv.DictReader(handle, delimiter=delimiter)
            if not reader.fieldnames:
                continue
            og_column = next(
                (
                    name
                    for name in reader.fieldnames
                    if name.strip().lower() in {"orthogroup", "og", "og_id"}
                ),
                reader.fieldnames[0],
            )
            for row in reader:
                og = (row.get(og_column) or "").strip()
                if og:
                    retained.add(og)

    return retained


def build_protein_og_map(orthogroups_path, retained_ogs=None):
    """
    Parse Orthogroups.tsv.

    All proteins are mapped to their assigned OG so retained queries can be
    compared with hits in any OG. Species locations and query evaluation are
    retained only for proteins in preliminary retained OGs.
    """
    protein_to_og = {}
    retained_protein_species = {}
    with open(orthogroups_path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError(f"Invalid or empty Orthogroups.tsv: {orthogroups_path}")
        for row in reader:
            if not row:
                continue
            og = row[0]
            is_retained = retained_ogs is None or og in retained_ogs
            for column, cell in enumerate(row[1:], start=1):
                species = header[column] if column < len(header) else None
                for protein in cell.split(","):
                    protein = protein.strip()
                    if protein:
                        previous = protein_to_og.get(protein)
                        if previous is not None and previous != og:
                            raise ValueError(
                                f"Protein {protein!r} is assigned to both {previous} and {og}"
                            )
                        protein_to_og[protein] = og
                        if is_retained:
                            if species is None:
                                raise ValueError(
                                    f"Orthogroup {og} has more fields than the header"
                                )
                            retained_protein_species[protein] = species
    return protein_to_og, retained_protein_species


def iter_blast_score_chunks(
    blast_dir,
    protein_to_og,
    retained_protein_species,
    bitscore_min,
    evalue_max,
):
    """Yield score tables one query species at a time to bound peak memory."""
    blast_files = sorted(glob.glob(os.path.join(blast_dir, "Blast[0-9]*_[0-9]*.txt")))
    if not blast_files:
        raise FileNotFoundError(f"No OrthoFinder Blast*_*.txt files found in {blast_dir}")

    files_by_query_species = defaultdict(list)
    for blast_file in blast_files:
        match = re.fullmatch(r"Blast(\d+)_(\d+)\.txt", os.path.basename(blast_file))
        if match is None:
            raise ValueError(f"Unexpected OrthoFinder BLAST filename: {blast_file}")
        files_by_query_species[int(match.group(1))].append(blast_file)

    files_scanned = 0
    for chunk_number, query_species in enumerate(
        sorted(files_by_query_species), start=1
    ):
        protein_og_scores = defaultdict(dict)
        for blast_file in files_by_query_species[query_species]:
            with open(blast_file) as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 12:
                        continue
                    query = parts[0]
                    if query not in retained_protein_species:
                        continue
                    subject = parts[1]
                    try:
                        evalue = float(parts[10])
                        bitscore = float(parts[11])
                    except ValueError:
                        continue

                    if bitscore < bitscore_min or evalue > evalue_max:
                        continue
                    if query == subject:
                        continue

                    query_og = protein_to_og[query]
                    subject_og = protein_to_og.get(subject)
                    if subject_og is None:
                        continue

                    if bitscore > protein_og_scores[query].get(subject_og, 0.0):
                        protein_og_scores[query][subject_og] = bitscore
            files_scanned += 1

        print(
            f"  scanned query species {chunk_number}/{len(files_by_query_species)} "
            f"({files_scanned}/{len(blast_files)} BLAST files)",
            file=sys.stderr,
        )
        yield query_species, protein_og_scores


def identify_chimeras(protein_og_scores, protein_to_og, ratio_threshold):
    """
    Identify chimeric proteins: those with strong hits to non-self OGs.

    Returns list of (protein, self_og, hit_ogs_str, ratios_str) tuples.
    """
    chimeras = []
    for protein, og_scores in protein_og_scores.items():
        self_og = protein_to_og.get(protein)
        if self_og is None:
            continue

        self_score = og_scores.get(self_og, 0)
        if self_score == 0:
            continue

        hit_ogs = []
        ratios = []
        for og, score in sorted(og_scores.items(), key=lambda x: -x[1]):
            if og == self_og:
                continue
            ratio = score / self_score
            if ratio >= ratio_threshold:
                hit_ogs.append(og)
                ratios.append(f"{ratio:.3f}")

        if hit_ogs:
            chimeras.append((protein, self_og, ";".join(hit_ogs), ";".join(ratios)))

    return chimeras


def remove_from_og_fastas(chimeras, og_seqs_dir):
    """Remove chimeric proteins from only their retained OG FASTA files."""
    chimeras_by_og = defaultdict(set)
    for protein, og, _hit_ogs, _ratios in chimeras:
        chimeras_by_og[og].add(protein)
    removed = 0

    for og, chimera_set in chimeras_by_og.items():
        fasta_file = os.path.join(og_seqs_dir, f"{og}.fa")
        if not os.path.isfile(fasta_file):
            raise FileNotFoundError(f"Missing retained OG FASTA: {fasta_file}")
        lines = []
        skip = False
        modified = False
        with open(fasta_file) as f:
            for line in f:
                if line.startswith(">"):
                    seqid = line[1:].strip().split()[0]
                    if seqid in chimera_set:
                        skip = True
                        modified = True
                        removed += 1
                    else:
                        skip = False
                if not skip:
                    lines.append(line)

        if modified:
            with open(fasta_file, "w") as f:
                f.writelines(lines)

    expected = sum(len(proteins) for proteins in chimeras_by_og.values())
    if removed != expected:
        raise ValueError(
            f"Removed {removed} of {expected} proteins flagged as chimeras from OG FASTAs"
        )
    return removed


def write_updated_gene_counts(
    source_path,
    destination_path,
    chimeras,
    retained_protein_species,
):
    """Stream gene counts while subtracting proteins removed as chimeras."""
    removals = defaultdict(int)
    for protein, og, _hit_ogs, _ratios in chimeras:
        species = retained_protein_species.get(protein)
        if species is None:
            raise ValueError(f"Species is unknown for retained protein {protein}")
        removals[(og, species)] += 1

    if not removals:
        shutil.copyfile(source_path, destination_path)
        return

    with open(source_path, newline="") as source, open(
        destination_path, "w", newline=""
    ) as destination:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        header = next(reader, None)
        if not header or "Total" not in header:
            raise ValueError(f"Invalid gene-count table: {source_path}")
        total_column = header.index("Total")
        unknown_species = {species for _og, species in removals} - set(header)
        if unknown_species:
            raise ValueError(
                "Chimera species absent from gene-count table: "
                + ", ".join(sorted(unknown_species))
            )
        writer.writerow(header)
        applied_removals = set()

        for row_number, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"Row {row_number} of {source_path} has {len(row)} fields; "
                    f"expected {len(header)}"
                )
            og = row[0]
            removed_total = 0
            for column in range(1, min(len(row), len(header))):
                species = header[column]
                if column == total_column:
                    continue
                count_removed = removals.get((og, species), 0)
                if count_removed:
                    original = int(row[column])
                    if count_removed > original:
                        raise ValueError(
                            f"Cannot remove {count_removed} proteins from {og}/{species}; "
                            f"gene-count table contains {original}"
                        )
                    row[column] = str(original - count_removed)
                    removed_total += count_removed
                    applied_removals.add((og, species))
            if removed_total:
                original_total = int(row[total_column])
                if removed_total > original_total:
                    raise ValueError(
                        f"Cannot remove {removed_total} proteins from {og}; "
                        f"total gene count is {original_total}"
                    )
                row[total_column] = str(original_total - removed_total)
            writer.writerow(row)

        missing_removals = set(removals) - applied_removals
        if missing_removals:
            preview = ", ".join(
                f"{og}/{species}" for og, species in sorted(missing_removals)[:10]
            )
            raise ValueError(
                f"Could not apply {len(missing_removals)} chimera count adjustments: {preview}"
            )


def write_updated_orthogroups(source_path, destination_path, chimeras):
    """Stream Orthogroups.tsv while removing proteins flagged as chimeras."""
    removals_by_og = defaultdict(set)
    for protein, og, _hit_ogs, _ratios in chimeras:
        removals_by_og[og].add(protein)

    if not removals_by_og:
        shutil.copyfile(source_path, destination_path)
        return

    removed = 0
    with open(source_path, newline="") as source, open(
        destination_path, "w", newline=""
    ) as destination:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        header = next(reader, None)
        if not header:
            raise ValueError(f"Invalid or empty Orthogroups.tsv: {source_path}")
        writer.writerow(header)

        for row_number, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) > len(header):
                raise ValueError(
                    f"Row {row_number} of {source_path} has {len(row)} fields; "
                    f"expected at most {len(header)}"
                )
            if len(row) < len(header):
                row.extend([""] * (len(header) - len(row)))

            og = row[0]
            og_removals = removals_by_og.get(og)
            if og_removals:
                for column in range(1, len(row)):
                    proteins = [
                        protein.strip()
                        for protein in row[column].split(",")
                        if protein.strip()
                    ]
                    retained = [
                        protein for protein in proteins if protein not in og_removals
                    ]
                    removed += len(proteins) - len(retained)
                    row[column] = ", ".join(retained)
            writer.writerow(row)

    expected = sum(len(proteins) for proteins in removals_by_og.values())
    if removed != expected:
        raise ValueError(
            f"Removed {removed} of {expected} proteins flagged as chimeras "
            "from Orthogroups.tsv"
        )


def main():
    args = parse_args()
    maximize_csv_field_size()

    if bool(args.gene_counts) != bool(args.updated_gene_counts):
        raise ValueError(
            "--gene-counts and --updated-gene-counts must be supplied together"
        )

    retained_ogs = read_retained_ogs(args.retained_og_files)
    if retained_ogs is not None:
        print(
            f"Loaded {len(retained_ogs)} preliminary retained orthogroups",
            file=sys.stderr,
        )

    print("Building protein → OG mapping...", file=sys.stderr)
    protein_to_og, retained_protein_species = build_protein_og_map(
        args.orthogroups, retained_ogs
    )
    print(f"  {len(protein_to_og)} proteins mapped", file=sys.stderr)
    print(
        f"  {len(retained_protein_species)} proteins retained for chimera evaluation",
        file=sys.stderr,
    )

    if retained_protein_species:
        print("Scanning BLAST results...", file=sys.stderr)
        chimeras = []
        qualifying_proteins = 0
        for _query_species, protein_og_scores in iter_blast_score_chunks(
            args.blast_dir,
            protein_to_og,
            retained_protein_species,
            args.bitscore_min,
            args.evalue_max,
        ):
            qualifying_proteins += len(protein_og_scores)
            chimeras.extend(
                identify_chimeras(
                    protein_og_scores, protein_to_og, args.ratio_threshold
                )
            )
    else:
        print("No retained proteins; skipping BLAST scan", file=sys.stderr)
        chimeras = []
        qualifying_proteins = 0
    print(f"  {qualifying_proteins} proteins with qualifying hits", file=sys.stderr)

    print("Chimera identification complete", file=sys.stderr)
    print(f"  {len(chimeras)} chimeric proteins flagged", file=sys.stderr)

    # Write report
    with open(args.report, "w") as f:
        f.write("protein\tassigned_og\thit_ogs\tbitscore_ratios\n")
        for protein, self_og, hit_ogs, ratios in chimeras:
            f.write(f"{protein}\t{self_og}\t{hit_ogs}\t{ratios}\n")

    # Remove chimeras from OG FASTA files
    if chimeras:
        removed = remove_from_og_fastas(chimeras, args.og_seqs_dir)
        print(f"  Removed {removed} chimeric sequences from OG FASTAs", file=sys.stderr)

    if args.gene_counts:
        print("Updating post-chimera gene counts...", file=sys.stderr)
        write_updated_gene_counts(
            args.gene_counts,
            args.updated_gene_counts,
            chimeras,
            retained_protein_species,
        )

    if args.updated_orthogroups:
        print("Updating post-chimera orthogroup membership...", file=sys.stderr)
        write_updated_orthogroups(
            args.orthogroups,
            args.updated_orthogroups,
            chimeras,
        )

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
