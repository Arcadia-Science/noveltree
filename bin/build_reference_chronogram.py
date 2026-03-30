#!/usr/bin/env python3
"""
Build a reference chronogram from TimeTree.org divergence times.

Queries TimeTree.org for pairwise divergence times among species,
constructs a UPGMA tree via scipy, and outputs a Newick ultrametric tree.

Species-level queries are attempted first, with genus-level fallback
when species pairs are not found.
"""

import argparse
import json
import logging
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# Rate limit: minimum seconds between NCBI requests
NCBI_RATE_LIMIT = 0.4
# Retry settings
MAX_RETRIES = 3
RETRY_BASE_DELAY = 1.0  # seconds, doubles each retry
# TimeTree API base URL (temple.edu is the actual API host)
TIMETREE_API = "http://timetree.temple.edu/api/pairwise"
# Lock to serialize NCBI Entrez calls (they share a global rate limit)
_ncbi_lock = threading.Lock()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a reference chronogram from TimeTree.org"
    )
    parser.add_argument(
        "--species-names",
        required=True,
        help="File with one species per line (Genus_species format)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output Newick file path",
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Email for NCBI Entrez queries",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel workers for TimeTree queries (default: 4)",
    )
    parser.add_argument(
        "--cache-dir",
        default=".",
        help="Directory for cache files (default: current directory)",
    )
    return parser.parse_args()


def load_cache(path):
    """Load a JSON cache file, returning empty dict if not found."""
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_cache(data, path):
    """Save data to a JSON cache file."""
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def ncbi_taxid_lookup(species_name, email, taxid_cache, cache_path):
    """
    Look up NCBI taxonomy ID for a species name via Entrez.

    Thread-safe: uses a global lock to serialize NCBI requests and
    respect rate limits. Retries with exponential backoff on transient errors.

    Parameters
    ----------
    species_name : str
        Species name in 'Genus species' format (spaces).
    email : str
        Email for NCBI Entrez.
    taxid_cache : dict
        Cache mapping species names to taxids.
    cache_path : str
        Path to save updated cache.

    Returns
    -------
    str or None
        NCBI taxonomy ID, or None if not found.
    """
    if species_name in taxid_cache:
        return taxid_cache[species_name]

    from Bio import Entrez

    Entrez.email = email

    for attempt in range(MAX_RETRIES):
        with _ncbi_lock:
            time.sleep(NCBI_RATE_LIMIT)
            try:
                handle = Entrez.esearch(db="taxonomy", term=species_name)
                record = Entrez.read(handle)
                handle.close()
            except Exception as e:
                delay = RETRY_BASE_DELAY * (2 ** attempt)
                if attempt < MAX_RETRIES - 1:
                    log.warning(
                        "NCBI lookup failed for '%s' (attempt %d/%d): %s. "
                        "Retrying in %.1fs...",
                        species_name, attempt + 1, MAX_RETRIES, e, delay,
                    )
                    time.sleep(delay)
                    continue
                else:
                    log.warning(
                        "NCBI lookup failed for '%s' after %d attempts: %s",
                        species_name, MAX_RETRIES, e,
                    )
                    return None

        id_list = record.get("IdList", [])
        if id_list:
            taxid = id_list[0]
            taxid_cache[species_name] = taxid
            save_cache(taxid_cache, cache_path)
            log.info("Taxid for '%s': %s", species_name, taxid)
            return taxid
        else:
            log.warning("No taxid found for '%s'", species_name)
            taxid_cache[species_name] = None
            save_cache(taxid_cache, cache_path)
            return None

    return None


def query_timetree(taxid_a, taxid_b):
    """
    Query TimeTree API for pairwise divergence time.

    The default (no suffix) endpoint returns a CSV with columns:
      taxon_a_id,taxon_b_id,scientific_name_a,scientific_name_b,
      all_total,precomputed_age,precomputed_ci_low,precomputed_ci_high,adjusted_age

    Retries with exponential backoff on transient HTTP errors.
    Returns divergence time in Mya (float), or None if not found.
    """
    url = f"{TIMETREE_API}/{taxid_a}/{taxid_b}"
    for attempt in range(MAX_RETRIES):
        req = Request(url)
        try:
            with urlopen(req, timeout=30) as resp:
                text = resp.read().decode().strip()
                if not text:
                    return None
                # Parse CSV: header line + data line
                lines = text.split("\n")
                if len(lines) < 2:
                    return None
                header = lines[0].split(",")
                values = lines[1].split(",")
                row = dict(zip(header, values))
                age = row.get("precomputed_age", "").strip()
                if age:
                    return float(age)
                return None
        except HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < MAX_RETRIES - 1:
                delay = RETRY_BASE_DELAY * (2 ** attempt)
                log.debug(
                    "TimeTree HTTP %d for %s vs %s, retrying in %.1fs...",
                    e.code, taxid_a, taxid_b, delay,
                )
                time.sleep(delay)
                continue
            log.debug("TimeTree query failed for %s vs %s: %s", taxid_a, taxid_b, e)
            return None
        except (URLError, ValueError, ConnectionError, TimeoutError) as e:
            if attempt < MAX_RETRIES - 1:
                delay = RETRY_BASE_DELAY * (2 ** attempt)
                log.debug(
                    "TimeTree error for %s vs %s: %s, retrying in %.1fs...",
                    taxid_a, taxid_b, e, delay,
                )
                time.sleep(delay)
                continue
            log.debug("TimeTree query failed for %s vs %s: %s", taxid_a, taxid_b, e)
            return None
    return None


def query_timetree_pair(sp_a, sp_b, taxid_map, genus_taxid_map, timetree_cache):
    """
    Query TimeTree for a species pair with genus-level fallback.

    Parameters
    ----------
    sp_a, sp_b : str
        Species names in 'Genus species' format.
    taxid_map : dict
        Species name -> NCBI taxid.
    genus_taxid_map : dict
        Genus name -> NCBI taxid (pre-computed).
    timetree_cache : dict
        Cache of previous TimeTree results.

    Returns
    -------
    tuple
        (sp_a, sp_b, divergence_time_mya_or_None)
    """
    cache_key = f"{sp_a}|{sp_b}"
    reverse_key = f"{sp_b}|{sp_a}"

    if cache_key in timetree_cache:
        return (sp_a, sp_b, timetree_cache[cache_key])
    if reverse_key in timetree_cache:
        return (sp_a, sp_b, timetree_cache[reverse_key])

    # Species-level query
    tid_a = taxid_map.get(sp_a)
    tid_b = taxid_map.get(sp_b)

    result = None
    if tid_a and tid_b:
        result = query_timetree(tid_a, tid_b)

    # Genus-level fallback (taxids pre-computed before parallel loop)
    if result is None:
        genus_a = sp_a.split()[0]
        genus_b = sp_b.split()[0]
        gtid_a = genus_taxid_map.get(genus_a)
        gtid_b = genus_taxid_map.get(genus_b)
        if gtid_a and gtid_b and gtid_a != gtid_b:
            result = query_timetree(gtid_a, gtid_b)
            if result is not None:
                log.info(
                    "Genus fallback succeeded for %s vs %s: %.1f Mya",
                    sp_a, sp_b, result,
                )

    timetree_cache[cache_key] = result
    return (sp_a, sp_b, result)


def build_distance_matrix(species_list, pair_times):
    """
    Build a complete distance matrix, dropping species with most missing pairs
    until all remaining pairs have values.

    Parameters
    ----------
    species_list : list of str
        Species names.
    pair_times : dict
        Mapping of "sp_a|sp_b" -> divergence_time.

    Returns
    -------
    tuple
        (remaining_species, distance_matrix as np.ndarray, dropped_species)
    """
    remaining = list(species_list)
    dropped = []

    while True:
        n = len(remaining)
        if n < 4:
            log.error(
                "Fewer than 4 species remain (%d) after dropping species with "
                "missing pairwise times. Cannot build chronogram.", n
            )
            sys.exit(1)

        # Count missing pairs per species
        missing_counts = {sp: 0 for sp in remaining}
        has_missing = False
        for i, sp_a in enumerate(remaining):
            for j in range(i + 1, n):
                sp_b = remaining[j]
                key = f"{sp_a}|{sp_b}"
                rkey = f"{sp_b}|{sp_a}"
                val = pair_times.get(key, pair_times.get(rkey))
                if val is None:
                    missing_counts[sp_a] += 1
                    missing_counts[sp_b] += 1
                    has_missing = True

        if not has_missing:
            break

        # Drop species with most missing pairs
        worst = max(missing_counts, key=missing_counts.get)
        log.info(
            "Dropping '%s' (%d missing pairs) to complete distance matrix",
            worst, missing_counts[worst],
        )
        remaining.remove(worst)
        dropped.append(worst)

    # Build the complete matrix
    n = len(remaining)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sp_a, sp_b = remaining[i], remaining[j]
            key = f"{sp_a}|{sp_b}"
            rkey = f"{sp_b}|{sp_a}"
            val = pair_times.get(key, pair_times.get(rkey))
            dist_matrix[i, j] = val
            dist_matrix[j, i] = val

    return remaining, dist_matrix, dropped


def linkage_to_newick(Z, labels):
    """
    Convert scipy linkage matrix to Newick string.

    Branch lengths are set so that node heights equal divergence times (Mya).
    Since UPGMA halves distances, we multiply by 2.

    Parameters
    ----------
    Z : np.ndarray
        Linkage matrix from scipy.cluster.hierarchy.linkage.
    labels : list of str
        Tip labels corresponding to original observations.

    Returns
    -------
    str
        Newick string with branch lengths.
    """
    tree_root = to_tree(Z)

    def _to_newick(node, parent_height=None):
        # UPGMA distances are halved, so node height = dist / 2
        # We multiply by 2 so heights represent actual Mya
        height = node.dist  # already the UPGMA merge height

        if node.is_leaf():
            label = labels[node.id]
            bl = parent_height - 0.0 if parent_height is not None else 0.0
            return f"{label}:{bl:.6f}"

        left = _to_newick(node.get_left(), height)
        right = _to_newick(node.get_right(), height)

        if parent_height is not None:
            bl = parent_height - height
            return f"({left},{right}):{bl:.6f}"
        else:
            return f"({left},{right})"

    return _to_newick(tree_root) + ";"


def main():
    args = parse_args()

    # Read species names
    with open(args.species_names) as f:
        raw_names = [line.strip() for line in f if line.strip()]

    # Normalize: convert hyphens to underscores, then to spaces for API queries
    species_names = [name.replace("-", " ").replace("_", " ") for name in raw_names]
    # Keep mapping back to underscore format for output
    underscore_names = [name.replace(" ", "_") for name in species_names]

    log.info("Loaded %d species", len(species_names))
    if len(species_names) < 4:
        log.error("Need at least 4 species to build a chronogram, got %d", len(species_names))
        sys.exit(1)

    # Setup cache paths
    taxid_cache_path = os.path.join(args.cache_dir, "taxid_cache.json")
    timetree_cache_path = os.path.join(args.cache_dir, "timetree_cache.json")

    taxid_cache = load_cache(taxid_cache_path)
    timetree_cache = load_cache(timetree_cache_path)

    # Step 1: NCBI taxid lookup for all species
    log.info("Looking up NCBI taxonomy IDs...")
    taxid_map = {}
    for sp in species_names:
        tid = ncbi_taxid_lookup(sp, args.email, taxid_cache, taxid_cache_path)
        if tid:
            taxid_map[sp] = tid
        else:
            log.warning("Could not find taxid for '%s'", sp)

    log.info("Found taxids for %d/%d species", len(taxid_map), len(species_names))

    # Step 2: Pre-compute genus-level taxids for fallback queries
    unique_genera = sorted(set(sp.split()[0] for sp in species_names))
    log.info("Looking up %d unique genus taxids for fallback...", len(unique_genera))
    genus_taxid_map = {}
    for genus in unique_genera:
        genus_taxid_map[genus] = ncbi_taxid_lookup(
            genus, args.email, taxid_cache, taxid_cache_path
        )
    save_cache(taxid_cache, taxid_cache_path)

    # Step 3: Query TimeTree for all pairwise combinations
    pairs = list(combinations(species_names, 2))
    log.info("Querying TimeTree.org for %d pairwise divergence times...", len(pairs))

    def _query_pair(pair):
        return query_timetree_pair(
            pair[0], pair[1], taxid_map, genus_taxid_map,
            timetree_cache,
        )

    found = 0
    missing = 0
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(_query_pair, p): p for p in pairs}
        for future in as_completed(futures):
            sp_a, sp_b, result = future.result()
            if result is not None:
                found += 1
            else:
                missing += 1

    # Save caches
    save_cache(timetree_cache, timetree_cache_path)
    save_cache(taxid_cache, taxid_cache_path)

    log.info("TimeTree results: %d found, %d missing", found, missing)

    # Step 3: Build distance matrix with greedy drop
    remaining_species, dist_matrix, dropped = build_distance_matrix(
        species_names, timetree_cache
    )
    log.info(
        "%d species retained, %d dropped", len(remaining_species), len(dropped)
    )

    # Map remaining species back to underscore format
    sp_to_underscore = dict(zip(species_names, underscore_names))
    remaining_underscore = [sp_to_underscore[sp] for sp in remaining_species]

    # Step 4: UPGMA via scipy
    # Convert full distance matrix to condensed form for scipy
    n = len(remaining_species)
    condensed = []
    for i in range(n):
        for j in range(i + 1, n):
            condensed.append(dist_matrix[i, j])
    condensed = np.array(condensed)

    Z = linkage(condensed, method="average")

    # Convert to Newick
    newick = linkage_to_newick(Z, remaining_underscore)

    # Write output
    with open(args.output, "w") as f:
        f.write(newick + "\n")
    log.info("Wrote chronogram to %s", args.output)

    # Write dropped species file
    dropped_path = os.path.join(args.cache_dir, "dropped_species.txt")
    dropped_underscore = [sp_to_underscore[sp] for sp in dropped]
    with open(dropped_path, "w") as f:
        if dropped_underscore:
            for sp in dropped_underscore:
                f.write(sp + "\n")
        else:
            f.write("# No species were dropped\n")
    log.info("Wrote dropped species list to %s", dropped_path)


if __name__ == "__main__":
    main()
