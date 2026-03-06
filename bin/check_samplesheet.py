#!/usr/bin/env python3

import os
import sys
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
    This function checks that the samplesheet follows the following structure:
    species,file,taxonomy,shallow_db,broad_db,mode,uniprot,mcl_test[,transdecoder,isoform,reference]

    The last three columns are optional. If absent, they default to 'no'.
    """

    # Valid extensions for protein FASTAs
    PROTEIN_EXTENSIONS = (".fasta", ".fa", ".fasta.gz", ".fa.gz")
    # Additional extensions allowed when transdecoder=yes (nucleotide input)
    NUCLEOTIDE_EXTENSIONS = (".fna", ".fna.gz")
    ALL_EXTENSIONS = PROTEIN_EXTENSIONS + NUCLEOTIDE_EXTENSIONS

    # Optional columns and their defaults
    OPTIONAL_COLS = ["transdecoder", "isoform", "reference"]
    VALID_YES_NO = {"yes", "no"}

    species_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        ## Check header
        MIN_COLS = 8
        HEADER = ["species", "file", "taxonomy", "shallow_db", "broad_db", "mode", "uniprot", "mcl_test"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            sys.exit(1)

        # Detect which optional columns are present
        extra_header = header[len(HEADER):]
        has_optional = {}
        for col in OPTIONAL_COLS:
            has_optional[col] = col in extra_header

        ## Check sample entries
        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                ## Check valid number of columns per row
                if len(lspl) < len(HEADER):
                    print_error(
                        f"Invalid number of columns (minimum = {len(HEADER)})!",
                        "Line",
                        line,
                    )

                num_cols = len([x for x in lspl[:MIN_COLS] if x])
                if num_cols < MIN_COLS:
                    print_error(
                        f"Invalid number of populated columns (minimum = {MIN_COLS})!",
                        "Line",
                        line,
                    )

                ## Check sample name entries
                species, file, taxonomy, shallow_db, broad_db, mode, uniprot, mcl_test = lspl[: len(HEADER)]
                # Normalize species names to use hyphens (e.g. "Homo sapiens" or
                # "Homo_sapiens" both become "Homo-sapiens"). This ensures
                # unambiguous parsing of tip labels (Species-name_ProteinID) since
                # the only underscore becomes the species-protein delimiter.
                original_species = species
                species = species.replace(" ", "-").replace("_", "-")
                if species != original_species:
                    print(f"WARNING: Species name normalized to hyphens: {original_species} -> {species}")
                if not species:
                    print_error("Sample entry has not been specified!", "Line", line)

                # Parse optional columns (default to 'no')
                optional_values = {}
                for col in OPTIONAL_COLS:
                    if has_optional[col]:
                        col_idx = len(HEADER) + extra_header.index(col)
                        val = lspl[col_idx].lower() if col_idx < len(lspl) and lspl[col_idx] else "no"
                    else:
                        val = "no"
                    if val not in VALID_YES_NO:
                        print_error(
                            f"'{col}' column must be 'yes' or 'no', got '{val}'!",
                            "Line",
                            line,
                        )
                    optional_values[col] = val

                transdecoder = optional_values["transdecoder"]
                isoform = optional_values["isoform"]
                reference = optional_values["reference"]

                # Validate: reference=yes requires uniprot=true
                if reference == "yes" and uniprot.lower() != "true":
                    print_error(
                        "reference=yes requires uniprot=true (reference proteomes come from UniProt)!",
                        "Line",
                        line,
                    )

                ## Check fasta file extension
                for fasta in [file]:
                    if fasta:
                        if fasta.find(" ") != -1:
                            print_error("fasta file contains spaces!", "Line", line)
                        # For URLs, validate extension from basename (strip query params)
                        check_name = fasta
                        if any(fasta.startswith(p) for p in ("http://", "https://", "ftp://", "s3://")):
                            check_name = fasta.split("?")[0].split("/")[-1]
                        # Determine valid extensions based on transdecoder flag
                        valid_exts = ALL_EXTENSIONS if transdecoder == "yes" else PROTEIN_EXTENSIONS
                        if not any(check_name.endswith(ext) for ext in valid_exts):
                            ext_str = "', '".join(valid_exts)
                            print_error(
                                f"File does not have a valid extension ('{ext_str}')!",
                                "Line",
                                line,
                            )

                ## populate sample data
                species_info = [
                    file,
                    taxonomy,
                    shallow_db,
                    broad_db,
                    mode,
                    uniprot,
                    mcl_test,
                    transdecoder,
                    isoform,
                    reference,
                ]

                ## Create species mapping dictionary
                if species not in species_mapping_dict:
                    species_mapping_dict[species] = [species_info]
                else:
                    if species_info in species_mapping_dict[species]:
                        print_error("Samplesheet contains duplicate rows!", "Line", line)
                    else:
                        species_mapping_dict[species].append(species_info)

    ## Write validated samplesheet with appropriate columns
    if len(species_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(["species", "file", "taxonomy", "shallow_db", "broad_db", "mode", "uniprot", "mcl_test",
                          "transdecoder", "isoform", "reference"])
                + "\n"
            )
            for species in sorted(species_mapping_dict.keys()):
                for idx, val in enumerate(species_mapping_dict[species]):
                    fout.write(",".join([f"{species}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
