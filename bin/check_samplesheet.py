#!/usr/bin/env python3

import os
import sys
import re
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat NovelTree samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    Validates and normalizes a NovelTree samplesheet.

    Required columns (positional, must be first 3):
        species, input_data, input_type

    Optional columns (any order after required, any subset):
        has_uniprot_ids, include_in_mcl_test, transdecoder,
        filter_isoforms, reference_proteome, busco_shallow, busco_broad

    All boolean columns use yes/no. Defaults: 'no' for booleans, 'NA' for BUSCO.
    Output CSV always contains all columns in a fixed order.
    """

    # Valid extensions for protein FASTAs
    PROTEIN_EXTENSIONS = (".fasta", ".fa", ".fasta.gz", ".fa.gz")
    NUCLEOTIDE_EXTENSIONS = (".fna", ".fna.gz")
    ALL_EXTENSIONS = PROTEIN_EXTENSIONS + NUCLEOTIDE_EXTENSIONS

    # Required columns (positional)
    REQUIRED = ["species", "input_data", "input_type"]

    # Optional columns with defaults
    OPTIONAL = {
        "has_uniprot_ids":    "no",
        "include_in_mcl_test": "no",
        "transdecoder":       "no",
        "filter_isoforms":    "no",
        "reference_proteome": "no",
        "busco_shallow":      "NA",
        "busco_broad":        "NA",
    }
    # Columns that accept free text (not validated as yes/no)
    FREE_TEXT = {"busco_shallow", "busco_broad"}
    VALID_YES_NO = {"yes", "no"}

    # Output column order (fixed, always written)
    OUTPUT_COLS = ["species", "input_data", "input_type"] + list(OPTIONAL.keys())

    species_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        header = [x.strip().strip('"') for x in fin.readline().strip().split(",")]

        # Validate required columns are first
        if header[: len(REQUIRED)] != REQUIRED:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)}")
            print(f"  Required (first {len(REQUIRED)} columns): {','.join(REQUIRED)}")
            sys.exit(1)

        # Detect which optional columns are present by name
        extra_header = header[len(REQUIRED):]
        has_col = {}
        for col in OPTIONAL:
            has_col[col] = col in extra_header

        for line in fin:
            if not line.strip():
                continue
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            if len(lspl) < len(REQUIRED):
                print_error(
                    f"Invalid number of columns (minimum = {len(REQUIRED)})!",
                    "Line", line,
                )

            # Unpack required columns
            species, input_data, input_type = lspl[: len(REQUIRED)]

            # Normalize species name to hyphens
            original_species = species
            species = species.replace(" ", "-").replace("_", "-")
            if species != original_species:
                print(f"WARNING: Species name normalized: {original_species} -> {species}")
            if not species:
                print_error("Species name is empty!", "Line", line)
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.\-]*", species):
                print_error(
                    "species must contain only letters, digits, dots, and hyphens "
                    "after normalization",
                    "Line",
                    line,
                )

            # Validate input_type
            if input_type.lower() not in ("proteins", "transcriptome"):
                print_error(
                    f"'input_type' must be 'proteins' or 'transcriptome', got '{input_type}'!",
                    "Line", line,
                )
            input_type = input_type.lower()

            # Parse optional columns
            opt = {}
            for col, default in OPTIONAL.items():
                if has_col[col]:
                    col_idx = len(REQUIRED) + extra_header.index(col)
                    val = lspl[col_idx].strip() if col_idx < len(lspl) and lspl[col_idx].strip() else default
                else:
                    val = default
                # Validate yes/no columns
                if col not in FREE_TEXT:
                    val = val.lower()
                    if val not in VALID_YES_NO:
                        print_error(
                            f"'{col}' must be 'yes' or 'no', got '{val}'!",
                            "Line", line,
                        )
                opt[col] = val

            # Cross-field validations
            if opt["reference_proteome"] == "yes" and opt["has_uniprot_ids"] != "yes":
                print_error(
                    "reference_proteome=yes requires has_uniprot_ids=yes!",
                    "Line", line,
                )
            if opt["include_in_mcl_test"] == "yes" and opt["has_uniprot_ids"] != "yes":
                print_error(
                    "include_in_mcl_test=yes requires has_uniprot_ids=yes!",
                    "Line", line,
                )

            # Validate input_data (file/URL/accession)
            if input_data:
                if " " in input_data:
                    print_error("input_data contains spaces!", "Line", line)
                is_url = any(input_data.startswith(p) for p in ("http://", "https://", "ftp://", "s3://"))
                is_ncbi = bool(re.match(r'^GC[AF]_\d+(\.\d+)?$', input_data))
                is_uniprot = bool(re.match(r'^UP\d{9,}$', input_data))
                is_remote = is_url or is_ncbi or is_uniprot
                check_name = input_data.split("?")[0].split("/")[-1] if is_url else input_data
                valid_exts = ALL_EXTENSIONS if opt["transdecoder"] == "yes" else PROTEIN_EXTENSIONS
                if not any(check_name.endswith(ext) for ext in valid_exts) and not is_remote:
                    print_error(
                        f"input_data does not have a valid extension ('{', '.join(valid_exts)}')!",
                        "Line", line,
                    )

            # Build row data
            row_data = [input_data, input_type] + [opt[col] for col in OPTIONAL]

            if species not in species_mapping_dict:
                species_mapping_dict[species] = [row_data]
            else:
                if row_data in species_mapping_dict[species]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    species_mapping_dict[species].append(row_data)

    # Write normalized output
    if len(species_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(OUTPUT_COLS) + "\n")
            for species in sorted(species_mapping_dict.keys()):
                for row_data in species_mapping_dict[species]:
                    fout.write(",".join([species] + row_data) + "\n")
    else:
        print_error("No entries to process!", f"Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
