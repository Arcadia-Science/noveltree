#!/usr/bin/env Rscript
# Load the required library
library(ape)

# Function to resolve polytomies and save the tree to a file
process_tree <- function(input_file, output_file) {
  # Read the Newick formatted tree from the input file
  tree <- read.tree(input_file)

  # Resolve polytomies using multi2di from the ape package
  resolved_tree <- multi2di(tree)

  # Write the resolved tree to the output file in Newick format
  write.tree(resolved_tree, file=output_file)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage: Rscript resolve_polytomies.R <input_file> <output_file>\n")
} else {
  # Assign input and output file names from command line arguments
  input_file <- args[1]
  output_file <- args[2]

  # Call the function with the input and output file names
  process_tree(input_file, output_file)
}
