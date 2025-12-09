#!/usr/bin/env python3
import argparse
import os
import tempfile

# Fix for bioservices in containers - set writable config/cache directories
# bioservices uses appdirs which looks for XDG_CONFIG_HOME and XDG_CACHE_HOME
tmpdir = tempfile.gettempdir()
os.environ['HOME'] = tmpdir
os.environ['XDG_CONFIG_HOME'] = tmpdir
os.environ['XDG_CACHE_HOME'] = tmpdir

from bioservices import UniProt
from bioservices import EUtils
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
import concurrent.futures
import pandas as pd

# Columns (query fields) that we will use when accessing UniProt protein annotations
COLUMNS = ['organism_name', 'organism_id', 'accession', 'xref_interpro', 'xref_oma']

    # Now define the different annotation sets - we will pull out tables for each.
ANNOTATION_SETS = {
    'cogeqc': ['organism_name', 'organism_id', 'accession', 'xref_interpro', 'xref_oma']
}

MIN_PROP_RETRIEVED = 0.5  # Minimum proportion of accessions that must be retrieved
MAX_RETRY_DEPTH = 5  # Maximum recursion depth for retries

# Suffixes for each annotation filename
ANNOT_NAMES = ['_cogeqc_annotations.tsv']

def fetch_batch(batch_accessions, organism_name, columns, depth=0):
    """Fetch annotations for a batch, recursively splitting on errors."""
    # Safety check: prevent excessive recursion
    if depth > MAX_RETRY_DEPTH:
        print(f"Warning: Max recursion depth reached for batch of size {len(batch_accessions)}")
        return []

    uniprot = UniProt()
    query = " OR ".join([f"accession:{acc}" for acc in batch_accessions])

    # Try to fetch results
    result = None
    try:
        result = uniprot.search(query, frmt="tsv", columns=",".join(columns), limit=None)
    except (KeyError, AttributeError):
        pass  # result stays None, will be handled below

    # Check if result is valid (covers both exception case and invalid return)
    if not result or not isinstance(result, str):
        # Invalid result - split and retry if possible
        if len(batch_accessions) > 1:
            mid = len(batch_accessions) // 2
            first_half = fetch_batch(batch_accessions[:mid], organism_name, columns, depth + 1)
            second_half = fetch_batch(batch_accessions[mid:], organism_name, columns, depth + 1)
            return first_half + second_half
        else:
            return []

    # Parse valid result
    annotations = []
    lines = result.strip().split('\n')[1:]
    for line in lines:
        fields = line.split('\t')
        annotations.append(fields)
    return annotations

def get_annotations(organism_name, input_file, columns, num_workers=None):
    annotations = []

    with open(input_file, 'r') as file:
        accessions = file.read().splitlines()

    # Split the accessions into batches of 50
    batch_size = 50
    batches = [accessions[i:i + batch_size] for i in range(0, len(accessions), batch_size)]
    
    if num_workers is None:
        num_workers = min(4, len(batches))
    
    # Print out some info:
    print(f"Pulling down annotations for {organism_name}")
    
    # Create a lock for thread-safe list updates
    annotations_lock = Lock()

    # Run the loop in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_results = [executor.submit(fetch_batch, batch_accessions, organism_name, columns) for batch_accessions in batches]
        for i, future in enumerate(concurrent.futures.as_completed(future_results)):
            batch_annotations = future.result()
            with annotations_lock:
                annotations.extend(batch_annotations)
            print(f"Completed batch {i + 1} of {len(batches)}")

    if len(annotations) < MIN_PROP_RETRIEVED * len(accessions):
        raise RuntimeError("Less than 50% of accessions were retrieved. Possible error in UniProt query.")

    return annotations

def main():
    parser = argparse.ArgumentParser(description='Retrieve UniProt annotations for a given species and list of protein accessions.')
    parser.add_argument('spp', help='Species to be queried (in Genus_species format)')
    parser.add_argument('ids', help='File containing the UniProt protein accessions')

    args = parser.parse_args()

    organism_name = args.spp
    input_file = args.ids

    # Begin by pulling down all annotations that we may want
    annotations = get_annotations(organism_name, input_file, COLUMNS)
    annots_df = pd.DataFrame(annotations, columns=COLUMNS)
    
    # Create a dictionary of data frames for each annotation set. 
    annot_df_dict = {key: annots_df[cols] for key, cols in ANNOTATION_SETS.items()}
    
    # Save each DataFrame subset to a TSV file with the appropriate prefix and suffix
    for idx, (key, df) in enumerate(annot_df_dict.items()):
        output_filename = f"{organism_name}{ANNOT_NAMES[idx]}"
        df.to_csv(output_filename, sep='\t', index=False)

if __name__ == "__main__":
    main()
    
