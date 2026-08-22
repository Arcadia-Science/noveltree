#!/usr/bin/env python3
"""Select a bounded multicopy family set for SpeciesRax inference.

This is intentionally a late, SpeciesRax-local guard. It limits pathological
family sizes without turning SpeciesRax into a single-copy analysis. Once the
cutoffs have been validated across production datasets, the same classification
should move into upstream orthogroup routing so excluded families are never
staged as species-tree inputs.
"""

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


MANIFEST_FIELDS = ["orthogroup", "gene_tree", "mapping"]
REPORT_FIELDS = [
    "orthogroup",
    "gene_tree",
    "mapping",
    "total_species",
    "species_present",
    "species_occupancy",
    "total_leaves",
    "mean_copies_present_species",
    "max_copies_any_species",
    "max_total_leaves",
    "selected",
    "exclusion_reasons",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest-glob",
        default="speciesrax_inputs_*.tsv",
        help="Glob selecting extracted SpeciesRax input manifests",
    )
    parser.add_argument(
        "--validation-report",
        required=True,
        type=Path,
        help="Tree-validation report produced by prepare_speciesrax_gene_trees.py",
    )
    parser.add_argument(
        "--expected-species-file",
        required=True,
        type=Path,
        help="One complete input-dataset species name per line",
    )
    parser.add_argument("--min-species-occupancy", required=True, type=float)
    parser.add_argument("--max-mean-copies", required=True, type=float)
    parser.add_argument("--max-copies-per-species", required=True, type=int)
    parser.add_argument("--max-total-leaves-factor", required=True, type=float)
    parser.add_argument("--selected-manifest", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--species-coverage-report", required=True, type=Path)
    return parser.parse_args()


def validate_cutoffs(args):
    if not 0 < args.min_species_occupancy <= 1:
        raise ValueError("min species occupancy must be in (0, 1]")
    if args.max_mean_copies < 1:
        raise ValueError("max mean copies must be at least 1")
    if args.max_copies_per_species < 1:
        raise ValueError("max copies per species must be at least 1")
    if args.max_total_leaves_factor < 1:
        raise ValueError("max total leaves factor must be at least 1")


def read_manifests(pattern):
    manifests = sorted(Path(".").glob(pattern))
    if not manifests:
        raise ValueError(f"No SpeciesRax manifests matched: {pattern}")

    rows = []
    seen = set()
    for manifest in manifests:
        with manifest.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or not set(MANIFEST_FIELDS).issubset(
                reader.fieldnames
            ):
                raise ValueError(f"Invalid SpeciesRax manifest header: {manifest}")
            for row in reader:
                orthogroup = row["orthogroup"]
                if orthogroup in seen:
                    raise ValueError(f"Duplicate SpeciesRax family: {orthogroup}")
                seen.add(orthogroup)
                rows.append({field: row[field] for field in MANIFEST_FIELDS})
    return rows


def read_validation_report(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"orthogroup", "leaves"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(f"Invalid SpeciesRax validation report: {path}")
        leaves = {}
        for row in reader:
            orthogroup = row["orthogroup"]
            if orthogroup in leaves:
                raise ValueError(
                    f"Duplicate family in SpeciesRax validation report: {orthogroup}"
                )
            leaves[orthogroup] = int(row["leaves"])
    return leaves


def read_expected_species(path):
    species = []
    seen = set()
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            name = line.strip()
            if not name:
                continue
            if name in seen:
                raise ValueError(
                    f"Duplicate expected species {name!r} at {path}:{line_number}"
                )
            seen.add(name)
            species.append(name)
    if len(species) < 3:
        raise ValueError(f"Expected species list contains only {len(species)} species")
    return species


def read_mapping(path, orthogroup, expected_species):
    mapping_path = Path(path)
    if not mapping_path.is_file():
        raise ValueError(f"Missing mapping for {orthogroup}: {mapping_path}")

    gene_to_species = {}
    normalized_rows = []
    corrected_rows = 0
    with mapping_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 2:
                raise ValueError(
                    f"Invalid mapping row for {orthogroup} at "
                    f"{mapping_path}:{line_number}"
                )
            gene, mapped_species = fields
            if gene in gene_to_species:
                raise ValueError(
                    f"Duplicate gene {gene!r} in mapping for {orthogroup}"
                )
            species, separator, _protein = gene.partition("_")
            if not separator or species not in expected_species:
                raise ValueError(
                    f"Cannot derive a run-level species from gene {gene!r} for "
                    f"{orthogroup} at {mapping_path}:{line_number}"
                )
            if mapped_species != species:
                corrected_rows += 1
            gene_to_species[gene] = species
            normalized_rows.append((gene, species))
    if not gene_to_species:
        raise ValueError(f"Empty mapping for {orthogroup}: {mapping_path}")

    # Mapping files are staged copies inside the SpeciesRax task. Normalize a
    # malformed cached mapping here so GeneRax consumes the same validated
    # Genus-species assignment used by the family-size filters. Stored upstream
    # mapping outputs are never modified.
    if corrected_rows:
        temporary = mapping_path.with_name(f"{mapping_path.name}.normalized.tmp")
        with temporary.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(normalized_rows)
        temporary.replace(mapping_path)

    return Counter(gene_to_species.values()), corrected_rows


def write_tsv(path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    validate_cutoffs(args)
    manifest_rows = read_manifests(args.manifest_glob)
    validated_leaves = read_validation_report(args.validation_report)
    expected_species = read_expected_species(args.expected_species_file)
    expected_species_set = set(expected_species)

    families = []
    all_species = set()
    corrected_mapping_rows = 0
    corrected_mapping_files = 0
    for row in manifest_rows:
        orthogroup = row["orthogroup"]
        if orthogroup not in validated_leaves:
            raise ValueError(f"Family lacks tree validation: {orthogroup}")
        copies, corrected_rows = read_mapping(
            row["mapping"], orthogroup, expected_species_set
        )
        corrected_mapping_rows += corrected_rows
        corrected_mapping_files += corrected_rows > 0
        mapping_leaves = sum(copies.values())
        if mapping_leaves != validated_leaves[orthogroup]:
            raise ValueError(
                f"Tree/mapping leaf-count mismatch for {orthogroup}: "
                f"tree={validated_leaves[orthogroup]}, mapping={mapping_leaves}"
            )
        all_species.update(copies)
        families.append((row, copies, mapping_leaves))

    if set(validated_leaves) != {
        row["orthogroup"] for row, _copies, _mapping_leaves in families
    }:
        raise ValueError("Tree-validation report and SpeciesRax manifests differ")

    unexpected_species = sorted(all_species - expected_species_set)
    missing_input_species = sorted(expected_species_set - all_species)
    if unexpected_species or missing_input_species:
        details = []
        if unexpected_species:
            details.append(
                "unexpected mapping species: " + ", ".join(unexpected_species[:10])
            )
        if missing_input_species:
            details.append(
                "input species absent from all staged families: "
                + ", ".join(missing_input_species[:10])
            )
        raise ValueError(
            "SpeciesRax mappings do not match the complete input species set; "
            + "; ".join(details)
        )
    total_species = len(expected_species)
    max_total_leaves = math.floor(args.max_total_leaves_factor * total_species)

    selected_manifest = []
    report_rows = []
    selected_species_families = Counter()
    selected_species_copies = Counter()
    reason_counts = Counter()

    for row, copies, total_leaves in families:
        species_present = len(copies)
        occupancy = species_present / total_species
        mean_copies = total_leaves / species_present
        max_copies = max(copies.values())
        reasons = []
        if occupancy < args.min_species_occupancy:
            reasons.append("species_occupancy_below_minimum")
        if mean_copies > args.max_mean_copies:
            reasons.append("mean_copies_above_maximum")
        if max_copies > args.max_copies_per_species:
            reasons.append("species_copy_count_above_maximum")
        if total_leaves > max_total_leaves:
            reasons.append("total_leaves_above_maximum")

        selected = not reasons
        if selected:
            selected_manifest.append(row)
            for species, copy_count in copies.items():
                selected_species_families[species] += 1
                selected_species_copies[species] += copy_count
        else:
            reason_counts.update(reasons)

        report_rows.append(
            {
                **row,
                "total_species": total_species,
                "species_present": species_present,
                "species_occupancy": f"{occupancy:.6f}",
                "total_leaves": total_leaves,
                "mean_copies_present_species": f"{mean_copies:.6f}",
                "max_copies_any_species": max_copies,
                "max_total_leaves": max_total_leaves,
                "selected": str(selected).lower(),
                "exclusion_reasons": ";".join(reasons),
            }
        )

    if not selected_manifest:
        raise ValueError("SpeciesRax family filtering retained no families")

    missing_species = sorted(expected_species_set - set(selected_species_families))
    if missing_species:
        preview = ", ".join(missing_species[:10])
        raise ValueError(
            f"SpeciesRax filtering removed all families for {len(missing_species)} "
            f"species: {preview}"
        )

    write_tsv(args.selected_manifest, MANIFEST_FIELDS, selected_manifest)
    write_tsv(args.report, REPORT_FIELDS, report_rows)
    coverage_rows = [
        {
            "species": species,
            "input_families": sum(species in copies for _, copies, _ in families),
            "selected_families": selected_species_families[species],
            "selected_copies": selected_species_copies[species],
        }
        for species in sorted(expected_species)
    ]
    write_tsv(
        args.species_coverage_report,
        ["species", "input_families", "selected_families", "selected_copies"],
        coverage_rows,
    )

    selected_leaves = sum(
        int(row["total_leaves"])
        for row in report_rows
        if row["selected"] == "true"
    )
    reason_summary = ", ".join(
        f"{reason}={count}" for reason, count in sorted(reason_counts.items())
    )
    print(
        "SpeciesRax family selection: "
        f"{len(selected_manifest)}/{len(families)} retained; "
        f"{selected_leaves} leaves; {total_species} species; "
        f"max leaves/family={max_total_leaves}"
    )
    if reason_summary:
        print(f"SpeciesRax family exclusions: {reason_summary}")
    if corrected_mapping_rows:
        print(
            "SpeciesRax mapping normalization: "
            f"{corrected_mapping_rows} rows corrected across "
            f"{corrected_mapping_files} task-local mapping files"
        )


if __name__ == "__main__":
    main()
