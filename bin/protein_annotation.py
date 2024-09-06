#!/usr/bin/env python3
import argparse
from bioservices import UniProt
from bioservices import EUtils
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
import concurrent.futures
import pandas as pd
import os

# Columns (query fields) that we will use when accessing UniProt protein annotations
COLUMNS = ['organism_name', 'organism_id', 'accession', 'xref_interpro', 'xref_oma']

    # Now define the different annotation sets - we will pull out tables for each.
ANNOTATION_SETS = {
    'cogeqc': ['organism_name', 'organism_id', 'accession', 'xref_interpro', 'xref_oma']
}

# Suffixes for each annotation filename
ANNOT_NAMES = ['_cogeqc_annotations.tsv']

def fetch_batch(batch_accessions, organism_name, columns):
    uniprot = UniProt()
    query = " OR ".join([f"accession:{acc}" for acc in batch_accessions])
    result = uniprot.search(query, frmt="tsv", columns=",".join(columns), limit=None)

    annotations = []
    if result and isinstance(result, str):
        lines = result.strip().split('\n')[1:]  # Ignore the header row
        for line in lines:
            fields = line.split('\t')
            annotations.append(fields)
    return annotations

def get_annotations(organism_name, input_file, columns, num_workers=None):
    uniprot = UniProt()
    annotations = []

    with open(input_file, 'r') as file:
        accessions = file.read().splitlines()

    # Split the accessions into batches of 100
    batch_size = 100
    batches = [accessions[i:i + batch_size] for i in range(0, len(accessions), batch_size)]
    
    if num_workers is None:
        num_workers = min(os.cpu_count(), len(batches))
    
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
    
