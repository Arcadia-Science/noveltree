#!/usr/bin/env python3
"""Translate GeneRax gene tree leaf labels to OrthoFinder internal IDs.

Reads SequenceIDs.txt to build a name->id mapping, then for each
*_reconciled_gft.newick file, tokenizes the Newick string and replaces
leaf labels with OrthoFinder's internal numeric IDs (e.g. 27_153).
Writes output to {results_dir}/Gene_Trees/Trees_ids/{OG}_tree_id.txt.
"""

import os
import sys
import re
import glob


def main():
    seq_ids_file = sys.argv[1]
    results_dir = sys.argv[2]

    # Build reverse lookup: gene_name -> internal_id
    name_to_id = {}
    with open(seq_ids_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            internal_id, gene_name = line.split(": ", 1)
            name_to_id[gene_name] = internal_id

    # Newick delimiters that separate tokens
    splitter = re.compile(r'([(),;:\[\]])')

    trees_dir = os.path.join(results_dir, "Gene_Trees", "Trees_ids")
    os.makedirs(trees_dir, exist_ok=True)

    n_trees = 0
    for tree_file in glob.glob("*_reconciled_gft.newick"):
        og = tree_file.replace("_reconciled_gft.newick", "")
        with open(tree_file) as f:
            newick = f.read().strip()

        # Split into tokens, replace leaf labels
        tokens = splitter.split(newick)
        translated = []
        for token in tokens:
            if token in name_to_id:
                translated.append(name_to_id[token])
            else:
                translated.append(token)

        out_path = os.path.join(trees_dir, og + "_tree_id.txt")
        with open(out_path, "w") as f:
            f.write("".join(translated) + "\n")
        n_trees += 1

    print("Converted {} gene trees to OrthoFinder format".format(n_trees), file=sys.stderr)


if __name__ == "__main__":
    main()
