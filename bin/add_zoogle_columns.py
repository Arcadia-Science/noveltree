#!/usr/bin/env python3
"""
Add zoogle-specific columns to the aggregated protein comparison table.

This script uses the exact same processing functions from the 2025-zoogle repository
to ensure consistency and reproducibility. It adds:
- HGNC gene symbols
- Within-gene-family percentiles
- Organism ranks
- Portfolio ranks
- Orthogroup homolog ranks and counts
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Import the processing functions from the zoogle_table module
# These are the exact same functions used in the 2025-zoogle portal
sys.path.insert(0, str(Path(__file__).parent))
from zoogle_table import (
    _merge_hgnc_gene_symbols,
    _append_organism_percentile_and_rank,
    _append_portfolio_rank,
    _append_organism_homolog_rank_and_count,
    ORGANISM_NAMES_TO_EXCLUDE,
)


def main():
    parser = argparse.ArgumentParser(
        description="Add zoogle-specific columns to aggregated protein comparison table"
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input aggregated protein comparison table (TSV or TSV.GZ)"
    )
    parser.add_argument(
        "--hgnc",
        required=True,
        type=Path,
        help="Processed HGNC dataset (TSV)"
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output zoogle table with additional columns (TSV or TSV.GZ)"
    )

    args = parser.parse_args()

    print(f"Reading aggregated table from {args.input}...", file=sys.stderr)

    # Read the input file (handle gzipped or plain text)
    if str(args.input).endswith('.gz'):
        os_dataset = pd.read_csv(args.input, sep="\t", compression='gzip')
    else:
        os_dataset = pd.read_csv(args.input, sep="\t")

    print(f"Loaded {len(os_dataset)} protein comparisons from {os_dataset.gene_family.nunique()} gene families", file=sys.stderr)

    # Ensure trait_dist is numeric
    os_dataset["trait_dist"] = pd.to_numeric(os_dataset["trait_dist"], errors="raise")

    # Drop unnecessary columns if they exist (including disease columns which are NA in Noveltree)
    # In RAAS, disease columns are populated even for non-disease proteins (ref_protein's disease info is duplicated)
    # In Noveltree with clinvar=NULL, these are all NA, so we drop them before the dropna check
    columns_to_drop = ["phylo_dist", "rank_trait_dist", "associated_gene", "disease_mim", "disease_name", "concept_id", "source_name", "source_id"]
    columns_to_drop = [col for col in columns_to_drop if col in os_dataset.columns]
    if columns_to_drop:
        os_dataset.drop(axis=1, labels=columns_to_drop, inplace=True)

    # Drop rows with missing values (after removing disease columns)
    num_rows_before = os_dataset.shape[0]
    os_dataset.dropna(axis=0, how="any", inplace=True)
    num_rows_after = os_dataset.shape[0]
    if num_rows_before > num_rows_after:
        print(f"Dropped {num_rows_before - num_rows_after} rows with missing values", file=sys.stderr)

    # Drop excluded organisms if any
    os_dataset = os_dataset.loc[
        ~os_dataset.nonref_species.str.lower().isin(
            [name.lower() for name in ORGANISM_NAMES_TO_EXCLUDE]
        )
    ].copy()

    # Apply zoogle processing functions (using exact same code as 2025-zoogle)
    print("Calculating organism percentiles and ranks...", file=sys.stderr)
    os_dataset = _append_organism_percentile_and_rank(os_dataset)

    print("Calculating portfolio ranks...", file=sys.stderr)
    os_dataset = _append_portfolio_rank(os_dataset)

    print("Calculating homolog ranks and counts...", file=sys.stderr)
    os_dataset = _append_organism_homolog_rank_and_count(os_dataset)

    # Merge HGNC gene symbols (using exact same code as 2025-zoogle)
    print("Merging HGNC gene symbols...", file=sys.stderr)
    num_rows_before = os_dataset.shape[0]
    os_dataset = _merge_hgnc_gene_symbols(os_dataset, args.hgnc)
    num_rows_after = os_dataset.shape[0]
    if num_rows_before != num_rows_after:
        print(f"Warning: Row count changed during HGNC merge: {num_rows_before} -> {num_rows_after}", file=sys.stderr)

    # Write output (compressed if output filename ends with .gz)
    print(f"Writing processed table to {args.output}...", file=sys.stderr)
    if str(args.output).endswith('.gz'):
        os_dataset.to_csv(args.output, sep="\t", index=False, compression='gzip')
    else:
        os_dataset.to_csv(args.output, sep="\t", index=False)

    print(f"Done! Processed table contains {len(os_dataset)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()