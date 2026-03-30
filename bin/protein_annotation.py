#!/usr/bin/env python3
"""Retrieve InterPro annotations from UniProt for a list of protein accessions.

Uses the UniProt ID Mapping API (POST-based, handles up to 100K IDs per job).
This is the officially recommended approach for bulk accession lookups. It
gracefully handles obsolete/merged accessions by simply returning results for
valid ones.
"""
import argparse
import sys
import time

import requests

# UniProt ID Mapping endpoints
IDMAPPING_RUN = "https://rest.uniprot.org/idmapping/run"
IDMAPPING_STATUS = "https://rest.uniprot.org/idmapping/status/{job_id}"
IDMAPPING_STREAM = "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job_id}"

FIELDS = "accession,organism_name,organism_id,xref_interpro"
OUTPUT_COLUMNS = ["organism_name", "organism_id", "accession", "xref_interpro"]

MAX_IDS_PER_JOB = 50000  # UniProt limit is 100K but use 50K for safety
POLL_INTERVAL = 3  # seconds between status checks
MAX_POLL_TIME = 600  # 10 minutes max wait per job
MAX_RETRIES = 2
MIN_PROP_RETRIEVED = 0.5

# Map from UniProt TSV header names to our internal names
HEADER_MAP = {
    "Entry": "accession",
    "From": "_from",  # ID mapping adds a "From" column
    "Organism": "organism_name",
    "Organism (ID)": "organism_id",
    "InterPro": "xref_interpro",
}


def submit_id_mapping(accessions):
    """Submit an ID mapping job. Returns the job ID or None on failure."""
    data = {
        "from": "UniProtKB_AC-ID",
        "to": "UniProtKB",
        "ids": ",".join(accessions),
    }
    try:
        resp = requests.post(IDMAPPING_RUN, data=data, timeout=60)
        resp.raise_for_status()
        return resp.json()["jobId"]
    except (requests.RequestException, KeyError) as e:
        print(f"  WARNING: Failed to submit ID mapping job: {e}", file=sys.stderr)
        return None


def poll_job(job_id):
    """Poll until the job is complete. Returns True if ready, False on timeout/error."""
    start = time.time()
    while time.time() - start < MAX_POLL_TIME:
        try:
            resp = requests.get(
                IDMAPPING_STATUS.format(job_id=job_id), timeout=30
            )
            resp.raise_for_status()
            result = resp.json()
            if "jobStatus" in result:
                status = result["jobStatus"]
                if status == "FINISHED":
                    return True
                if status in ("ERROR", "FAILED"):
                    print(f"  WARNING: Job {job_id} failed: {result}", file=sys.stderr)
                    return False
            elif "results" in result or "failedIds" in result:
                # Job is complete (some endpoints return results directly)
                return True
        except requests.RequestException as e:
            print(f"  Poll error for {job_id}: {e}", file=sys.stderr)
        time.sleep(POLL_INTERVAL)

    print(f"  WARNING: Job {job_id} timed out after {MAX_POLL_TIME}s", file=sys.stderr)
    return False


def fetch_results(job_id):
    """Fetch TSV results for a completed job. Returns list of row dicts."""
    params = {
        "format": "tsv",
        "fields": FIELDS,
    }
    try:
        resp = requests.get(
            IDMAPPING_STREAM.format(job_id=job_id),
            params=params,
            timeout=300,
        )
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"  WARNING: Failed to fetch results for {job_id}: {e}", file=sys.stderr)
        return []

    rows = []
    lines = resp.text.strip().split("\n")
    if len(lines) < 2:
        return rows

    header = lines[0].split("\t")
    col_indices = {}
    for i, col in enumerate(header):
        if col in HEADER_MAP:
            col_indices[HEADER_MAP[col]] = i

    # Check we have the columns we need (excluding _from which is just the mapping source)
    needed = set(OUTPUT_COLUMNS)
    have = set(col_indices.keys()) - {"_from"}
    missing = needed - have
    if missing:
        print(f"  WARNING: Missing columns in response: {missing}", file=sys.stderr)
        return []

    for line in lines[1:]:
        fields = line.split("\t")
        row = {}
        for name, idx in col_indices.items():
            if name == "_from":
                continue
            row[name] = fields[idx] if idx < len(fields) else ""
        rows.append(row)

    return rows


def get_annotations(accessions):
    """Retrieve InterPro annotations for all accessions via ID Mapping API."""
    all_rows = []
    n_chunks = (len(accessions) + MAX_IDS_PER_JOB - 1) // MAX_IDS_PER_JOB

    for chunk_i in range(0, len(accessions), MAX_IDS_PER_JOB):
        chunk = accessions[chunk_i : chunk_i + MAX_IDS_PER_JOB]
        chunk_num = chunk_i // MAX_IDS_PER_JOB + 1
        print(f"  Chunk {chunk_num}/{n_chunks} ({len(chunk)} accessions)")

        rows = None
        for attempt in range(MAX_RETRIES + 1):
            job_id = submit_id_mapping(chunk)
            if job_id is None:
                if attempt < MAX_RETRIES:
                    wait = 2 ** (attempt + 1)
                    print(f"  Retry {attempt + 1}/{MAX_RETRIES} in {wait}s",
                          file=sys.stderr)
                    time.sleep(wait)
                continue

            print(f"    Job submitted: {job_id}")
            if poll_job(job_id):
                rows = fetch_results(job_id)
                if rows is not None:
                    break
            elif attempt < MAX_RETRIES:
                wait = 2 ** (attempt + 1)
                print(f"  Retry {attempt + 1}/{MAX_RETRIES} in {wait}s",
                      file=sys.stderr)
                time.sleep(wait)

        if rows:
            all_rows.extend(rows)
        else:
            print(f"  WARNING: Failed chunk {chunk_num} after {MAX_RETRIES} retries",
                  file=sys.stderr)

    return all_rows


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve UniProt InterPro annotations for protein accessions."
    )
    parser.add_argument("spp", help="Species name")
    parser.add_argument(
        "ids", help="File with UniProt protein accessions (one per line)"
    )
    args = parser.parse_args()

    with open(args.ids) as f:
        accessions = [line.strip() for line in f if line.strip()]

    if not accessions:
        print(f"WARNING: No accessions found in {args.ids}", file=sys.stderr)
        output_file = f"{args.spp}_cogeqc_annotations.tsv"
        with open(output_file, "w") as f:
            f.write("\t".join(OUTPUT_COLUMNS) + "\n")
        return

    print(f"Retrieving InterPro annotations for {args.spp} "
          f"({len(accessions)} accessions)")
    rows = get_annotations(accessions)
    print(f"  Retrieved {len(rows)} annotations")

    if len(rows) < MIN_PROP_RETRIEVED * len(accessions):
        print(
            f"  WARNING: Only {len(rows)}/{len(accessions)} accessions "
            f"retrieved ({len(rows)/len(accessions)*100:.1f}%). "
            f"Check accession format.",
            file=sys.stderr,
        )

    output_file = f"{args.spp}_cogeqc_annotations.tsv"
    with open(output_file, "w") as f:
        f.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            values = [row.get(col, "") for col in OUTPUT_COLUMNS]
            f.write("\t".join(values) + "\n")

    print(f"  Written to {output_file}")


if __name__ == "__main__":
    main()
