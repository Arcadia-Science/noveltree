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
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Flag cross-OG chimeric proteins")
    parser.add_argument("--orthogroups", required=True, help="Path to Orthogroups.tsv")
    parser.add_argument("--blast_dir", required=True, help="Directory containing DIAMOND BLAST result files")
    parser.add_argument("--og_seqs_dir", required=True, help="Path to Orthogroup_Sequences/ directory")
    parser.add_argument("--report", required=True, help="Output chimera report TSV")
    parser.add_argument("--bitscore_min", type=float, default=100, help="Minimum bitscore for a hit (default: 100)")
    parser.add_argument("--evalue_max", type=float, default=1e-10, help="Maximum evalue for a hit (default: 1e-10)")
    parser.add_argument("--ratio_threshold", type=float, default=0.5, help="Min ratio of best non-self OG hit to self OG hit (default: 0.5)")
    return parser.parse_args()


def build_protein_og_map(orthogroups_path):
    """Parse Orthogroups.tsv to build protein → OG mapping."""
    protein_to_og = {}
    with open(orthogroups_path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            og = row[0]
            for cell in row[1:]:
                for protein in cell.split(", "):
                    protein = protein.strip()
                    if protein:
                        protein_to_og[protein] = og
    return protein_to_og


def scan_blast_results(blast_dir, protein_to_og, bitscore_min, evalue_max):
    """
    Scan DIAMOND BLAST results for cross-OG hits.

    Returns dict: protein → {og → best_bitscore} for proteins with hits
    to multiple OGs.
    """
    protein_og_scores = defaultdict(lambda: defaultdict(float))

    # OrthoFinder BLAST files are named like Blast0_1.txt
    blast_files = glob.glob(os.path.join(blast_dir, "Blast*_*.txt"))
    if not blast_files:
        # Try alternative naming
        blast_files = glob.glob(os.path.join(blast_dir, "*.txt"))

    for bf in blast_files:
        with open(bf) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 12:
                    continue
                query = parts[0]
                subject = parts[1]
                evalue = float(parts[10])
                bitscore = float(parts[11])

                if bitscore < bitscore_min or evalue > evalue_max:
                    continue

                if query == subject:
                    continue

                query_og = protein_to_og.get(query)
                subject_og = protein_to_og.get(subject)

                if query_og is None or subject_og is None:
                    continue

                # Track best bitscore per OG for this query
                if bitscore > protein_og_scores[query][subject_og]:
                    protein_og_scores[query][subject_og] = bitscore

    return protein_og_scores


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


def remove_from_og_fastas(chimera_proteins, og_seqs_dir):
    """Remove chimeric proteins from Orthogroup_Sequences FASTA files."""
    chimera_set = set(chimera_proteins)
    removed = 0

    for fasta_file in glob.glob(os.path.join(og_seqs_dir, "*.fa")):
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

    return removed


def main():
    args = parse_args()

    print("Building protein → OG mapping...", file=sys.stderr)
    protein_to_og = build_protein_og_map(args.orthogroups)
    print(f"  {len(protein_to_og)} proteins mapped", file=sys.stderr)

    print("Scanning BLAST results...", file=sys.stderr)
    protein_og_scores = scan_blast_results(
        args.blast_dir, protein_to_og, args.bitscore_min, args.evalue_max
    )
    print(f"  {len(protein_og_scores)} proteins with qualifying hits", file=sys.stderr)

    print("Identifying chimeras...", file=sys.stderr)
    chimeras = identify_chimeras(protein_og_scores, protein_to_og, args.ratio_threshold)
    print(f"  {len(chimeras)} chimeric proteins flagged", file=sys.stderr)

    # Write report
    with open(args.report, "w") as f:
        f.write("protein\tassigned_og\thit_ogs\tbitscore_ratios\n")
        for protein, self_og, hit_ogs, ratios in chimeras:
            f.write(f"{protein}\t{self_og}\t{hit_ogs}\t{ratios}\n")

    # Remove chimeras from OG FASTA files
    if chimeras:
        chimera_proteins = [c[0] for c in chimeras]
        removed = remove_from_og_fastas(chimera_proteins, args.og_seqs_dir)
        print(f"  Removed {removed} chimeric sequences from OG FASTAs", file=sys.stderr)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
