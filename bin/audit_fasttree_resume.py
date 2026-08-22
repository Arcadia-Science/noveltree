#!/usr/bin/env python3
"""Audit stored FastTree outputs before resuming a production run.

The audit is read-only. It identifies the families routed directly to
FASTTREE_TIER2 from the OrthoFinder metadata, compares their expected output
names with the stable tree store, and can enforce an exact known resume state.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import urlparse


@dataclass(frozen=True)
class StoredObject:
    name: str
    size: int
    etag: str
    last_modified: str


@dataclass(frozen=True)
class AuditResult:
    expected: tuple[str, ...]
    records: tuple[StoredObject, ...]
    missing: tuple[str, ...]
    fingerprint: str

    @property
    def stored(self) -> tuple[str, ...]:
        return tuple(record.name for record in self.records)


def _aws(*args: str) -> str:
    completed = subprocess.run(
        ("aws", *args),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"AWS CLI failed: {message}")
    return completed.stdout


def _read_text(path: str) -> str:
    if path.startswith("s3://"):
        return _aws("s3", "cp", path, "-", "--no-progress")
    return Path(path).read_text()


def expected_fasttree_outputs(
    metadata_path: str,
    *,
    family_set: str = "gene_tree",
    min_sequences: int = 1000,
) -> tuple[str, ...]:
    reader = csv.DictReader(io.StringIO(_read_text(metadata_path)), delimiter="\t")
    required = {"orthogroup", "family_set", "n_seq"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        missing = sorted(required.difference(reader.fieldnames or ()))
        raise ValueError(f"Metadata is missing required column(s): {', '.join(missing)}")

    outputs = {
        f"{row['orthogroup']}_famsa_clipkit_ft.newick"
        for row in reader
        if row["family_set"] == family_set and int(row["n_seq"]) > min_sequences
    }
    return tuple(sorted(outputs))


def _parse_s3_uri(uri: str) -> tuple[str, str]:
    parsed = urlparse(uri)
    if parsed.scheme != "s3" or not parsed.netloc:
        raise ValueError(f"Invalid S3 URI: {uri}")
    prefix = parsed.path.lstrip("/").rstrip("/") + "/"
    return parsed.netloc, prefix


def stored_objects(tree_store: str) -> dict[str, StoredObject]:
    if tree_store.startswith("s3://"):
        bucket, prefix = _parse_s3_uri(tree_store)
        payload = json.loads(
            _aws(
                "s3api",
                "list-objects-v2",
                "--bucket",
                bucket,
                "--prefix",
                prefix,
                "--page-size",
                "1000",
                "--output",
                "json",
            )
        )
        records = payload.get("Contents", ())
        return {
            Path(record["Key"]).name: StoredObject(
                name=Path(record["Key"]).name,
                size=int(record["Size"]),
                etag=str(record.get("ETag", "")).strip('"'),
                last_modified=str(record.get("LastModified", "")),
            )
            for record in records
            if int(record["Size"]) > 0
        }

    root = Path(tree_store)
    if not root.is_dir():
        raise ValueError(f"Tree store is not a directory: {tree_store}")
    return {
        path.name: StoredObject(
            name=path.name,
            size=path.stat().st_size,
            etag="",
            last_modified=str(path.stat().st_mtime_ns),
        )
        for path in root.iterdir()
        if path.is_file() and path.stat().st_size > 0
    }


def audit_fasttree_resume(
    metadata_path: str,
    tree_store: str,
    *,
    min_sequences: int = 1000,
) -> AuditResult:
    expected = expected_fasttree_outputs(
        metadata_path,
        min_sequences=min_sequences,
    )
    objects = stored_objects(tree_store)
    stored = tuple(name for name in expected if name in objects)
    missing = tuple(name for name in expected if name not in objects)

    digest = hashlib.sha256()
    for name in stored:
        record = objects[name]
        digest.update(
            f"{name}\t{record.size}\t{record.etag}\t{record.last_modified}\n".encode()
        )
    return AuditResult(
        expected,
        tuple(objects[name] for name in stored),
        missing,
        digest.hexdigest(),
    )


def write_snapshot(path: str, result: AuditResult) -> None:
    payload = {
        "schema_version": 1,
        "fingerprint": result.fingerprint,
        "objects": [
            {
                "name": record.name,
                "size": record.size,
                "etag": record.etag,
                "last_modified": record.last_modified,
            }
            for record in result.records
        ],
    }
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def verify_snapshot(path: str, result: AuditResult) -> tuple[str, ...]:
    payload = json.loads(Path(path).read_text())
    if payload.get("schema_version") != 1 or not isinstance(payload.get("objects"), list):
        raise ValueError(f"Invalid FastTree resume snapshot: {path}")

    current = {record.name: record for record in result.records}
    changed: list[str] = []
    for saved in payload["objects"]:
        name = saved.get("name", "")
        observed = current.get(name)
        if observed is None or (
            observed.size != saved.get("size")
            or observed.etag != saved.get("etag")
            or observed.last_modified != saved.get("last_modified")
        ):
            changed.append(name)
    return tuple(sorted(changed))


def _orthogroups(names: Iterable[str]) -> tuple[str, ...]:
    return tuple(name.split("_", 1)[0] for name in names)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--tree-store", required=True)
    parser.add_argument("--min-sequences", type=int, default=1000)
    parser.add_argument("--expect-candidates", type=int)
    parser.add_argument("--expect-stored", type=int)
    snapshot = parser.add_mutually_exclusive_group()
    snapshot.add_argument(
        "--write-snapshot",
        metavar="JSON",
        help="Write the stored-object fingerprint manifest for post-run verification",
    )
    snapshot.add_argument(
        "--verify-snapshot",
        metavar="JSON",
        help="Fail if any object recorded in a previous snapshot changed or disappeared",
    )
    parser.add_argument(
        "--expect-missing",
        nargs="*",
        default=None,
        metavar="ORTHOGROUP",
        help="Exact orthogroup IDs expected to be missing",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = audit_fasttree_resume(
            args.metadata,
            args.tree_store,
            min_sequences=args.min_sequences,
        )
    except (OSError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    print(f"Expected FASTTREE_TIER2 outputs: {len(result.expected)}")
    print(f"Stored nonempty outputs: {len(result.stored)}")
    print(f"Missing outputs: {len(result.missing)}")
    print(f"Stored-object fingerprint: {result.fingerprint}")
    if result.missing:
        print("Missing orthogroups: " + ", ".join(_orthogroups(result.missing)))

    errors: list[str] = []
    if args.expect_candidates is not None and len(result.expected) != args.expect_candidates:
        errors.append(
            f"expected {args.expect_candidates} candidates, observed {len(result.expected)}"
        )
    if args.expect_stored is not None and len(result.stored) != args.expect_stored:
        errors.append(f"expected {args.expect_stored} stored outputs, observed {len(result.stored)}")
    if args.expect_missing is not None:
        expected_missing = tuple(sorted(args.expect_missing))
        observed_missing = tuple(sorted(_orthogroups(result.missing)))
        if observed_missing != expected_missing:
            errors.append(
                "missing set differs: expected "
                + ",".join(expected_missing)
                + "; observed "
                + ",".join(observed_missing)
            )

    if args.verify_snapshot:
        try:
            changed = verify_snapshot(args.verify_snapshot, result)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            errors.append(str(error))
        else:
            if changed:
                errors.append("stored output(s) changed or disappeared: " + ",".join(changed))

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print("Resume preflight FAILED; do not launch Nextflow.", file=sys.stderr)
        return 1

    if args.write_snapshot:
        try:
            write_snapshot(args.write_snapshot, result)
        except OSError as error:
            print(f"ERROR: Could not write snapshot: {error}", file=sys.stderr)
            return 2
        print(f"Wrote stored-output snapshot: {args.write_snapshot}")

    print("Resume preflight passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
