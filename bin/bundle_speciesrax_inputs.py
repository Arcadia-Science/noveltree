#!/usr/bin/env python3
"""Validate and bundle a bounded shard of SpeciesRax family inputs."""

import argparse
import csv
import tarfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("."))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest-name", required=True)
    return parser.parse_args()


def orthogroup(path):
    if "_" not in path.name:
        raise ValueError(f"Cannot derive orthogroup from input filename: {path.name}")
    return path.name.split("_", 1)[0]


def index_unique(paths, kind):
    result = {}
    for path in sorted(paths):
        og = orthogroup(path)
        if og in result:
            raise ValueError(
                f"Multiple {kind} files for {og}: {result[og].name}, {path.name}"
            )
        result[og] = path
    return result


def collect_families(input_dir):
    trees = index_unique(input_dir.glob("*.newick"), "gene tree")
    mappings = index_unique(input_dir.glob("*_map.link"), "mapping")
    if not trees:
        raise ValueError(f"No SpeciesRax gene trees found in {input_dir}")

    all_ogs = set(trees) | set(mappings)
    incomplete = [
        og
        for og in sorted(all_ogs)
        if og not in trees or og not in mappings
    ]
    if incomplete:
        preview = ", ".join(incomplete[:10])
        raise ValueError(
            f"{len(incomplete)} SpeciesRax families lack exactly one gene tree "
            f"and mapping; first: {preview}"
        )
    return [(og, trees[og], mappings[og]) for og in sorted(all_ogs)]


def main():
    args = parse_args()
    families = collect_families(args.input_dir)
    manifest_path = args.input_dir / args.manifest_name
    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orthogroup", "gene_tree", "mapping"])
        for og, tree, mapping in families:
            writer.writerow([og, tree.name, mapping.name])

    with tarfile.open(args.output, "w") as archive:
        archive.add(manifest_path, arcname=manifest_path.name, recursive=False)
        for _og, tree, mapping in families:
            for path in (tree, mapping):
                archive.add(path, arcname=path.name, recursive=False)
    print(f"Bundled {len(families)} SpeciesRax families into {args.output}")


if __name__ == "__main__":
    main()
