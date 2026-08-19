#!/usr/bin/env python3
"""Create, restore, and remove durable OrthoFinder MCL checkpoints.

The checkpoint contains the raw orthogroup membership tables and per-orthogroup
FASTAs. Those are sufficient for NovelTree postprocessing without duplicating
OrthoFinder's much larger transient graph/working directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


CHECKPOINT_FORMAT_VERSION = 1
NOT_RESTORABLE = 10
BUFFER_SIZE = 8 * 1024 * 1024


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(BUFFER_SIZE):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint_inputs(args: argparse.Namespace) -> int:
    paths = sorted((Path(value) for value in args.inputs), key=lambda path: path.name)
    metadata_paths = sorted(
        (Path(value) for value in args.metadata_only_inputs), key=lambda path: path.name
    )
    if not paths:
        raise ValueError("At least one checkpoint input is required")

    basenames: set[str] = set()
    digest = hashlib.sha256()
    metadata = {
        "checkpoint_format": CHECKPOINT_FORMAT_VERSION,
        "orthofinder_version": args.orthofinder_version,
        "inflation": args.inflation,
        "orthofinder_options": args.orthofinder_options,
        "extra_args": args.extra_args,
    }
    digest.update(json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode())
    digest.update(b"\0")

    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Checkpoint fingerprint input is missing: {path}")
        if path.name in basenames:
            raise ValueError(f"Duplicate checkpoint input basename: {path.name}")
        basenames.add(path.name)
        record = f"{path.name}\0{path.stat().st_size}\0{sha256_file(path)}\n"
        digest.update(record.encode())

    # BLAST bundle contents are deterministically generated from the fully
    # hashed proteomes and sequence maps. Recording their names and sizes
    # verifies the expected bundle set without rereading the enormous all-v-all
    # search results solely to construct the checkpoint key.
    for path in metadata_paths:
        if not path.is_file():
            raise FileNotFoundError(f"Checkpoint metadata input is missing: {path}")
        if path.name in basenames:
            raise ValueError(f"Duplicate checkpoint input basename: {path.name}")
        basenames.add(path.name)
        digest.update(f"metadata\0{path.name}\0{path.stat().st_size}\n".encode())

    print(digest.hexdigest())
    return 0


def checkpoint_names(fingerprint: str) -> tuple[str, str]:
    prefix = f"orthofinder_mcl_{fingerprint}"
    return f"{prefix}.tar.gz", f"{prefix}.json"


def checkpoint_location(root: str, name: str) -> str:
    return f"{root.rstrip('/')}/{name}"


def is_s3(path: str) -> bool:
    return path.startswith("s3://")


def aws_cli() -> str:
    candidates = [
        os.environ.get("AWS_CLI"),
        shutil.which("aws"),
        "/home/ec2-user/miniconda/bin/aws",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    raise RuntimeError(
        "An AWS CLI executable is required for an s3:// OrthoFinder checkpoint"
    )


def run_aws(
    arguments: list[str], *, allow_missing: bool = False, attempts: int = 8
) -> bool:
    command = [aws_cli(), *arguments, "--only-show-errors"]
    environment = os.environ.copy()
    environment.setdefault("AWS_RETRY_MODE", "adaptive")
    environment.setdefault("AWS_MAX_ATTEMPTS", "10")

    for attempt in range(1, attempts + 1):
        result = subprocess.run(
            command,
            env=environment,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode == 0:
            return True
        error = result.stderr.strip()
        missing = any(
            token in error
            for token in ("404", "NoSuchKey", "does not exist", "not found")
        )
        if allow_missing and missing:
            return False
        if attempt == attempts:
            raise RuntimeError(
                f"AWS CLI failed after {attempts} attempts: {' '.join(command)}\n{error}"
            )
        delay = min(60, 2 ** (attempt - 1))
        print(
            f"AWS checkpoint transfer failed (attempt {attempt}/{attempts}); "
            f"retrying in {delay}s: {error}",
            file=sys.stderr,
        )
        time.sleep(delay)
    return False


def validate_results(results_dir: Path) -> None:
    orthogroups = results_dir / "Orthogroups"
    membership = orthogroups / "Orthogroups.tsv"
    counts = orthogroups / "Orthogroups.GeneCount.tsv"
    sequences = results_dir / "Orthogroup_Sequences"
    for path in (membership, counts):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Required OrthoFinder checkpoint file is missing: {path}")
    if not sequences.is_dir() or next(sequences.glob("*.fa"), None) is None:
        raise ValueError(
            f"OrthoFinder checkpoint has no orthogroup FASTAs: {sequences}"
        )


def create_archive(results_dir: Path, archive: Path) -> None:
    validate_results(results_dir)
    cwd = Path.cwd().resolve()
    results_dir = results_dir.resolve()
    try:
        relative_results = results_dir.relative_to(cwd)
    except ValueError as error:
        raise ValueError(
            f"Results directory must be below the task directory: {results_dir}"
        ) from error

    members = [
        str(relative_results / "Orthogroups"),
        str(relative_results / "Orthogroup_Sequences"),
    ]
    with archive.open("wb") as archive_handle:
        tar_process = subprocess.Popen(
            ["tar", "-C", str(cwd), "-cf", "-", *members],
            stdout=subprocess.PIPE,
        )
        assert tar_process.stdout is not None
        gzip_process = subprocess.Popen(
            ["gzip", "-1"], stdin=tar_process.stdout, stdout=archive_handle
        )
        tar_process.stdout.close()
        gzip_status = gzip_process.wait()
        tar_status = tar_process.wait()
    if tar_status != 0 or gzip_status != 0:
        archive.unlink(missing_ok=True)
        raise RuntimeError(
            f"Checkpoint archive creation failed (tar={tar_status}, gzip={gzip_status})"
        )


def write_receipt(
    receipt: Path, root: str, fingerprint: str, archive_name: str, manifest_name: str
) -> None:
    payload = {
        "checkpoint_format": CHECKPOINT_FORMAT_VERSION,
        "checkpoint_root": root,
        "fingerprint": fingerprint,
        "archive": checkpoint_location(root, archive_name),
        "manifest": checkpoint_location(root, manifest_name),
    }
    receipt.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def save_checkpoint(args: argparse.Namespace) -> int:
    results_dir = Path(args.results_dir)
    receipt = Path(args.receipt)
    archive_name, manifest_name = checkpoint_names(args.fingerprint)
    archive = Path(f".{archive_name}")
    manifest = Path(f".{manifest_name}")

    create_archive(results_dir, archive)
    payload = {
        "checkpoint_format": CHECKPOINT_FORMAT_VERSION,
        "fingerprint": args.fingerprint,
        "archive_name": archive_name,
        "archive_sha256": sha256_file(archive),
        "contents": ["Orthogroups", "Orthogroup_Sequences"],
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    try:
        if is_s3(args.root):
            # The manifest is the completion marker and must be uploaded last.
            run_aws(
                ["s3", "cp", str(archive), checkpoint_location(args.root, archive_name)]
            )
            run_aws(
                ["s3", "cp", str(manifest), checkpoint_location(args.root, manifest_name)]
            )
        else:
            root = Path(args.root)
            root.mkdir(parents=True, exist_ok=True)
            archive_target = root / archive_name
            manifest_target = root / manifest_name
            archive_temporary = root / f".{archive_name}.{os.getpid()}.tmp"
            manifest_temporary = root / f".{manifest_name}.{os.getpid()}.tmp"
            shutil.copyfile(archive, archive_temporary)
            os.replace(archive_temporary, archive_target)
            shutil.copyfile(manifest, manifest_temporary)
            os.replace(manifest_temporary, manifest_target)
        write_receipt(
            receipt, args.root, args.fingerprint, archive_name, manifest_name
        )
    finally:
        archive.unlink(missing_ok=True)
        manifest.unlink(missing_ok=True)

    print(
        f"Saved OrthoFinder MCL checkpoint: "
        f"{checkpoint_location(args.root, archive_name)}"
    )
    return 0


def obtain_checkpoint_file(source: str, destination: Path, *, marker: bool) -> bool:
    if is_s3(source):
        return run_aws(
            ["s3", "cp", source, str(destination)], allow_missing=marker
        )
    source_path = Path(source)
    if not source_path.is_file():
        return False
    shutil.copyfile(source_path, destination)
    return True


def remove_partial_restore(results_dir: Path) -> None:
    if results_dir.exists():
        shutil.rmtree(results_dir)


def restore_checkpoint(args: argparse.Namespace) -> int:
    results_dir = Path(args.results_dir)
    receipt = Path(args.receipt)
    archive_name, manifest_name = checkpoint_names(args.fingerprint)
    archive_source = checkpoint_location(args.root, archive_name)
    manifest_source = checkpoint_location(args.root, manifest_name)
    archive = Path(f".{archive_name}.restore")
    manifest = Path(f".{manifest_name}.restore")

    try:
        manifest_available = obtain_checkpoint_file(
            manifest_source, manifest, marker=True
        )
        payload = None
        if manifest_available:
            payload = json.loads(manifest.read_text())
            if (
                payload.get("checkpoint_format") != CHECKPOINT_FORMAT_VERSION
                or payload.get("fingerprint") != args.fingerprint
                or payload.get("archive_name") != archive_name
            ):
                print("OrthoFinder checkpoint manifest is incompatible", file=sys.stderr)
                return NOT_RESTORABLE
            obtain_checkpoint_file(archive_source, archive, marker=False)
            if sha256_file(archive) != payload.get("archive_sha256"):
                print("OrthoFinder checkpoint archive checksum failed", file=sys.stderr)
                return NOT_RESTORABLE
        else:
            # S3 exposes a multipart upload only after it is complete. If the
            # archive upload succeeded but the small completion manifest did
            # not, validate and recover the content-addressed archive rather
            # than spending another day rebuilding it.
            if not obtain_checkpoint_file(archive_source, archive, marker=True):
                print("No complete OrthoFinder MCL checkpoint was found")
                return NOT_RESTORABLE
            print(
                "Checkpoint manifest is absent; validating the completed archive",
                file=sys.stderr,
            )

        remove_partial_restore(results_dir)
        extraction = subprocess.run(["tar", "-xzf", str(archive)])
        if extraction.returncode != 0:
            remove_partial_restore(results_dir)
            print("OrthoFinder checkpoint archive could not be extracted", file=sys.stderr)
            return NOT_RESTORABLE
        try:
            validate_results(results_dir)
        except ValueError as error:
            remove_partial_restore(results_dir)
            print(str(error), file=sys.stderr)
            return NOT_RESTORABLE
        write_receipt(
            receipt, args.root, args.fingerprint, archive_name, manifest_name
        )
    finally:
        archive.unlink(missing_ok=True)
        manifest.unlink(missing_ok=True)

    print(f"Restored OrthoFinder MCL checkpoint: {archive_source}")
    return 0


def remove_location(location: str) -> None:
    if is_s3(location):
        run_aws(["s3", "rm", location])
    else:
        Path(location).unlink(missing_ok=True)


def cleanup_checkpoint(args: argparse.Namespace) -> int:
    payload = json.loads(Path(args.receipt).read_text())
    if payload.get("checkpoint_format") != CHECKPOINT_FORMAT_VERSION:
        raise ValueError("Unsupported OrthoFinder checkpoint receipt")
    # Delete the completion marker first so a partial cleanup cannot be restored.
    remove_location(payload["manifest"])
    remove_location(payload["archive"])
    print(f"Removed completed OrthoFinder MCL checkpoint {payload['fingerprint']}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    fingerprint = subparsers.add_parser("fingerprint")
    fingerprint.add_argument("--orthofinder-version", required=True)
    fingerprint.add_argument("--inflation", required=True)
    fingerprint.add_argument("--orthofinder-options", required=True)
    fingerprint.add_argument("--extra-args", default="")
    fingerprint.add_argument("--metadata-only-inputs", nargs="*", default=[])
    fingerprint.add_argument("inputs", nargs="+")
    fingerprint.set_defaults(function=fingerprint_inputs)

    for command, function in (("save", save_checkpoint), ("restore", restore_checkpoint)):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--root", required=True)
        subparser.add_argument("--fingerprint", required=True)
        subparser.add_argument("--results-dir", required=True)
        subparser.add_argument("--receipt", required=True)
        subparser.set_defaults(function=function)

    cleanup = subparsers.add_parser("cleanup")
    cleanup.add_argument("--receipt", required=True)
    cleanup.set_defaults(function=cleanup_checkpoint)
    return parser


def main() -> int:
    parser = build_parser()
    arguments = parser.parse_args()
    try:
        return arguments.function(arguments)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
