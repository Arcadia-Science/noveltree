#!/usr/bin/env Rscript

# Time-calibrate a consensus species tree using a reference time-calibrated tree
# This script uses geiger::congruify.phylo() to transfer time-calibration from
# a reference tree to a target consensus tree.

# Load required libraries
library(ape)
library(geiger)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(paste(
    "Usage: time_calibrate_species_tree.R <consensus_tree> <reference_tree> <output_tree> <calibration_method>",
    "\n  consensus_tree: Path to consensus species tree (Newick format)",
    "\n  reference_tree: Path to time-calibrated reference tree (Newick format)",
    "\n  output_tree: Path for output time-calibrated tree",
    "\n  calibration_method: Method for time calibration (treePL or PATHd8)",
    sep = "\n"
  ))
}

consensus_tree_path <- args[1]
reference_tree_path <- args[2]
output_tree_path <- args[3]
calibration_method <- args[4]

# Validate calibration method
if (!calibration_method %in% c("treePL", "PATHd8")) {
  stop("calibration_method must be either 'treePL' or 'PATHd8'")
}

cat("Reading consensus species tree from:", consensus_tree_path, "\n")
consensus_tree <- read.tree(consensus_tree_path)

cat("Reading reference time tree from:", reference_tree_path, "\n")
reference_tree <- read.tree(reference_tree_path)

# Check that reference tree is ultrametric (time-calibrated)
if (!is.ultrametric(reference_tree)) {
  stop("Reference tree must be ultrametric (time-calibrated)")
}

cat("Reference tree has", length(reference_tree$tip.label), "species\n")
cat("Consensus tree has", length(consensus_tree$tip.label), "species\n")

# Find overlapping species
shared_species <- intersect(consensus_tree$tip.label, reference_tree$tip.label)
cat("Number of shared species:", length(shared_species), "\n")

if (length(shared_species) < 2) {
  stop("At least 2 shared species are required for time calibration")
}

cat("Shared species:", paste(shared_species, collapse = ", "), "\n")

# Create taxonomy matrix for congruify
# This maps species names from consensus tree to reference tree
# For now, we assume direct mapping (species names match)
taxonomy <- matrix(consensus_tree$tip.label,
                   dimnames = list(consensus_tree$tip.label, NULL),
                   ncol = 1)

cat("\nApplying time-calibration using", calibration_method, "\n")

# Set up environment for PATHd8 if needed
if (calibration_method == "PATHd8") {
  pathd8_path <- paste0(
    Sys.getenv("CONDA_PREFIX"), "/bin:",
    Sys.getenv("CONDA_PREFIX"), "/lib64:",
    Sys.getenv("LD_LIBRARY_PATH")
  )
  Sys.setenv(LD_LIBRARY_PATH = pathd8_path)
}

# Apply congruify to transfer time-calibration
# This will:
# 1. Prune both trees to shared species
# 2. Transfer time-calibration from reference to consensus
# 3. Use treePL or PATHd8 to estimate times for species not in reference
calibrated_tree <- suppressWarnings(
  geiger::congruify.phylo(
    reference = reference_tree,
    target = consensus_tree,
    scale = calibration_method,
    taxonomy = taxonomy
  )$phy
)

cat("Time-calibration complete\n")
cat("Calibrated tree has", length(calibrated_tree$tip.label), "species\n")

# Check if the resulting tree is ultrametric
if (is.ultrametric(calibrated_tree)) {
  cat("Output tree is ultrametric (time-calibrated)\n")
} else {
  warning("Output tree is not ultrametric - may have numerical precision issues")
}

# Write the calibrated tree
cat("Writing time-calibrated tree to:", output_tree_path, "\n")
write.tree(calibrated_tree, file = output_tree_path)

cat("Done!\n")
