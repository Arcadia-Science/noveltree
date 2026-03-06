#!/usr/bin/env python3
"""
Filter FASTA to keep the longest isoform per gene.

Supports four header patterns:
1. NCBI RefSeq: gene=GENENAME in header → group by gene name
2. NCBI eukaryotic: [gene=LOC...] [protein_id=XP_...] → group by gene
3. Trinity assemblies: TRINITY_DN*_c*_g*_i* → group by gene (before _i)
4. TransDecoder output: GENE.X~~Y.pN → group by gene (before .p)

If no pattern matches, each sequence is treated as its own group (all kept).
"""

import re
import sys


def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def extract_gene_id(header):
    """Extract a gene grouping ID from a FASTA header."""
    seqid = header.split()[0]
    full_header = header

    # Pattern 1: gene=GENENAME (NCBI RefSeq style)
    m = re.search(r"\bgene=(\S+)", full_header)
    if m:
        return m.group(1)

    # Pattern 2: [gene=LOC...] (NCBI eukaryotic)
    m = re.search(r"\[gene=([^\]]+)\]", full_header)
    if m:
        return m.group(1)

    # Pattern 3: Trinity — TRINITY_DN{d}_c{d}_g{d}_i{d}
    m = re.match(r"(TRINITY_DN\d+_c\d+_g\d+)_i\d+", seqid)
    if m:
        return m.group(1)

    # Pattern 4: TransDecoder — GENE.X~~Y.pN
    m = re.match(r"(.+)\.p\d+$", seqid)
    if m:
        return m.group(1)

    # No pattern matched — unique group per sequence
    return seqid


def filter_isoforms(input_path, output_path):
    """Keep only the longest sequence per gene group."""
    gene_best = {}  # gene_id → (header, sequence, length)

    for header, seq in parse_fasta(input_path):
        gene_id = extract_gene_id(header)
        seq_len = len(seq)
        if gene_id not in gene_best or seq_len > gene_best[gene_id][2]:
            gene_best[gene_id] = (header, seq, seq_len)

    with open(output_path, "w") as out:
        for gene_id in gene_best:
            header, seq, _ = gene_best[gene_id]
            out.write(f">{header}\n{seq}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.fasta output.fasta", file=sys.stderr)
        sys.exit(1)
    filter_isoforms(sys.argv[1], sys.argv[2])
