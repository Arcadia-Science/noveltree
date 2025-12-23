#!/usr/bin/env Rscript

# Time-calibrate a species tree using a reference time-calibrated tree
# This script uses geiger::congruify.phylo() to transfer time-calibration from
# a reference tree to a target species tree.

# Load required libraries
library(ape)
library(geiger)
library(phytools)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(paste(
    "Usage: time_calibrate_species_tree.R <species_tree> <reference_tree> <output_tree> <calibration_method>",
    "\n  species_tree: Path to species tree (Newick format)",
    "\n  reference_tree: Path to time-calibrated reference tree (Newick format)",
    "\n  output_tree: Path for output time-calibrated tree",
    "\n  calibration_method: Method for time calibration (treePL or PATHd8)",
    sep = "\n"
  ))
}

species_tree_path <- args[1]
reference_tree_path <- args[2]
output_tree_path <- args[3]
calibration_method <- args[4]

# Validate calibration method
if (!calibration_method %in% c("treePL", "PATHd8")) {
  stop("calibration_method must be either 'treePL' or 'PATHd8'")
}

cat("Reading species tree from:", species_tree_path, "\n")
species_tree <- read.tree(species_tree_path)

cat("Reading reference time tree from:", reference_tree_path, "\n")
reference_tree <- read.tree(reference_tree_path)

# Check that reference tree is ultrametric (time-calibrated)
# Force ultrametricity to handle rounding errors if needed
if (!is.ultrametric(reference_tree)) {
  cat("Reference tree not strictly ultrametric - forcing ultrametricity to handle rounding errors\n")
  reference_tree <- force.ultrametric(reference_tree)
  if (!is.ultrametric(reference_tree)) {
    stop("Reference tree cannot be made ultrametric")
  }
}

cat("Reference tree has", length(reference_tree$tip.label), "species\n")
cat("Species tree has", length(species_tree$tip.label), "species\n")

# Find overlapping species
shared_species <- intersect(species_tree$tip.label, reference_tree$tip.label)
cat("Number of shared species:", length(shared_species), "\n")

if (length(shared_species) < 2) {
  stop("At least 2 shared species are required for time calibration")
}

cat("Shared species:", paste(shared_species, collapse = ", "), "\n")

# Create taxonomy matrix for congruify
# This maps species names from species tree to reference tree
# For now, we assume direct mapping (species names match)
taxonomy <- matrix(species_tree$tip.label,
                   dimnames = list(species_tree$tip.label, NULL),
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
# 2. Transfer time-calibration from reference to species tree
# 3. Use treePL or PATHd8 to estimate times for species not in reference
calibrated_tree <- suppressWarnings(
  geiger::congruify.phylo(
    reference = reference_tree,
    target = species_tree,
    scale = calibration_method,
    taxonomy = taxonomy
  )$phy
)

cat("Time-calibration complete\n")
cat("Calibrated tree has", length(calibrated_tree$tip.label), "species\n")

# Check if the resulting tree is ultrametric
# Force ultrametricity to handle rounding errors if needed
if (!is.ultrametric(calibrated_tree)) {
  cat("Output tree not strictly ultrametric - forcing ultrametricity to handle rounding errors\n")
  calibrated_tree <- force.ultrametric(calibrated_tree)
  if (!is.ultrametric(calibrated_tree)) {
    warning("Output tree cannot be made ultrametric - may have structural issues")
  }
} else {
  cat("Output tree is ultrametric (time-calibrated)\n")
}

# Write the calibrated tree
cat("Writing time-calibrated tree to:", output_tree_path, "\n")
write.tree(calibrated_tree, file = output_tree_path)

cat("Done!\n")
