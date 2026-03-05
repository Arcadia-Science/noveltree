#!/usr/bin/env Rscript

# Time-calibrate a species tree using a reference time-calibrated tree.
# Extracts calibration points from the reference chronogram, then runs treePL
# with cross-validation to find optimal smoothing, followed by a final ADOLC pass.

# Load required libraries
suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phangorn)
})

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

# Normalize reference tree tip labels: convert underscores to hyphens to match
# the pipeline's internal species name convention (e.g. Homo_sapiens -> Homo-sapiens).
# The species tree already uses hyphens from the samplesheet normalization.
reference_tree$tip.label <- gsub("_", "-", reference_tree$tip.label)
cat("Normalized reference tree tip labels (underscores -> hyphens)\n")

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

# ============================================================================
# Step 1: Extract calibrations from reference chronogram
# ============================================================================

cat("\nExtracting calibrations from reference chronogram...\n")

ref_depths <- node.depth.edgelength(reference_tree)
ref_root_depth <- max(ref_depths)

n_ref_tips <- length(reference_tree$tip.label)
n_ref_internal <- reference_tree$Nnode

calibrations <- data.frame(
  spp_mrca = integer(0),
  age_mya = numeric(0),
  tipA = character(0),
  tipB = character(0),
  n_desc_tips = integer(0),
  stringsAsFactors = FALSE
)

for (i in seq_len(n_ref_internal)) {
  ref_node <- n_ref_tips + i

  # Get descendant tip labels in the reference tree
  desc_tips <- reference_tree$tip.label[unlist(
    Descendants(reference_tree, ref_node, type = "tips")
  )]

  # Intersect with species tree tip labels
  shared_desc <- intersect(desc_tips, species_tree$tip.label)
  if (length(shared_desc) < 2) next

  # Find MRCA in species tree
  spp_mrca <- getMRCA(species_tree, shared_desc)
  if (is.null(spp_mrca)) next

  # Reference node age = root depth - node depth
  age_mya <- ref_root_depth - ref_depths[ref_node]
  if (is.na(age_mya) || age_mya <= 0) next

  # Pick two representative tips (must exist in species tree)
  tipA <- shared_desc[1]
  tipB <- shared_desc[length(shared_desc)]

  calibrations <- rbind(calibrations, data.frame(
    spp_mrca = spp_mrca,
    age_mya = age_mya,
    tipA = tipA,
    tipB = tipB,
    n_desc_tips = length(shared_desc),
    stringsAsFactors = FALSE
  ))
}

cat("Raw calibrations extracted:", nrow(calibrations), "\n")

# Deduplicate: one calibration per species-tree MRCA node (keep entry with most descendant tips)
if (nrow(calibrations) > 0) {
  calibrations <- calibrations[order(-calibrations$n_desc_tips), ]
  calibrations <- calibrations[!duplicated(calibrations$spp_mrca), ]
}

cat("Calibrations after dedup:", nrow(calibrations), "\n")

# Remove parent-child age conflicts:
# Sort by node depth (root-to-tip), drop calibrations where child age > parent age
if (nrow(calibrations) > 1) {
  spp_depths <- node.depth.edgelength(species_tree)
  calibrations$node_depth <- spp_depths[calibrations$spp_mrca]
  calibrations <- calibrations[order(calibrations$node_depth), ]

  dropped <- c()
  for (i in 2:nrow(calibrations)) {
    current_node <- calibrations$spp_mrca[i]
    current_age <- calibrations$age_mya[i]

    for (j in 1:(i - 1)) {
      if (j %in% dropped) next
      ancestor_node <- calibrations$spp_mrca[j]
      ancestor_age <- calibrations$age_mya[j]

      if (ancestor_node %in% Ancestors(species_tree, current_node, type = "all")) {
        if (current_age > ancestor_age) {
          cat("  Dropping conflicting calibration: node", current_node,
              "(age", current_age, ") > ancestor node", ancestor_node,
              "(age", ancestor_age, ")\n")
          dropped <- c(dropped, i)
          break
        }
      }
    }
  }
  if (length(dropped) > 0) {
    calibrations <- calibrations[-dropped, ]
  }
  calibrations$node_depth <- NULL
}

cat("Calibrations after conflict resolution:", nrow(calibrations), "\n")

if (nrow(calibrations) < 2) {
  stop("Need at least 2 calibration points from the reference tree, got ", nrow(calibrations))
}

# ============================================================================
# Step 2: Build calibration lines for treePL config
# ============================================================================

# Use fixed-point calibrations (min = max = age) since these are known TimeTree ages
cal_lines <- c()
for (i in seq_len(nrow(calibrations))) {
  cal_name <- paste0("cal", i)
  cal_lines <- c(cal_lines,
    paste("mrca =", cal_name, calibrations$tipA[i], calibrations$tipB[i]),
    paste("min =", cal_name, calibrations$age_mya[i]),
    paste("max =", cal_name, calibrations$age_mya[i])
  )
}

cat("Using", nrow(calibrations), "fixed-point calibrations\n")

# ============================================================================
# Step 3: Write species tree for treePL
# ============================================================================

# Fix negative or zero branch lengths (treePL cannot handle them)
if (any(species_tree$edge.length <= 0)) {
  min_pos <- min(species_tree$edge.length[species_tree$edge.length > 0])
  n_fixed <- sum(species_tree$edge.length <= 0)
  species_tree$edge.length[species_tree$edge.length <= 0] <- min_pos * 0.01
  cat("Fixed", n_fixed, "non-positive branch lengths (set to", min_pos * 0.01, ")\n")
}

tree_file <- "spp_tree_for_treepl.newick"
write.tree(species_tree, file = tree_file)

# ============================================================================
# Step 4: treePL CV pass
# ============================================================================

cat("\nApplying time-calibration using", calibration_method, "\n")

if (calibration_method == "treePL") {

  cv_config_file <- "spp_tree_treepl_cv.config"
  cv_out_file <- "spp_tree_cv_out.newick"
  cv_config <- c(
    paste("treefile =", tree_file),
    "numsites = 10000",
    "smooth = 100",
    cal_lines,
    paste("outfile =", cv_out_file),
    "opt = 1",
    "optad = 1",
    "cvstart = 1000",
    "cvstop = 0.1",
    "cviter = 3",
    "cv"
  )
  writeLines(cv_config, cv_config_file)

  cat("Running treePL cross-validation...\n")
  cv_output <- tryCatch({
    system2("/opt/conda/bin/treePL", args = cv_config_file,
            stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    cat("treePL CV pass failed:", e$message, "\n")
    NULL
  })

  # Parse optimal smoothing from CV output
  optimal_smooth <- 100  # default fallback
  if (!is.null(cv_output)) {
    smooth_line <- grep("Optimal smoothing value", cv_output, value = TRUE)
    if (length(smooth_line) > 0) {
      smooth_val <- as.numeric(sub(".*: *", "", smooth_line[1]))
      if (!is.na(smooth_val) && smooth_val > 0) {
        optimal_smooth <- smooth_val
        cat("Optimal smoothing from CV:", optimal_smooth, "\n")
      } else {
        cat("Could not parse optimal smoothing, using default: 100\n")
      }
    } else {
      cat("No optimal smoothing line found in CV output, using default: 100\n")
    }
  } else {
    cat("CV pass failed, using default smoothing: 100\n")
  }

  # ============================================================================
  # Step 5: treePL final pass with optimal smoothing
  # ============================================================================

  final_config_file <- "spp_tree_treepl_final.config"
  final_out_file <- "spp_tree_treepl_dated.newick"
  final_config <- c(
    paste("treefile =", tree_file),
    "numsites = 10000",
    paste("smooth =", optimal_smooth),
    cal_lines,
    paste("outfile =", final_out_file),
    "opt = 1",
    "optad = 1"
  )
  writeLines(final_config, final_config_file)

  cat("Running treePL final pass (smooth =", optimal_smooth, ")...\n")
  treepl_result <- tryCatch({
    system2("/opt/conda/bin/treePL", args = final_config_file,
            stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    cat("treePL final pass failed:", e$message, "\n")
    NULL
  })

  if (!file.exists(final_out_file) || file.info(final_out_file)$size == 0) {
    cat("treePL output:\n")
    if (!is.null(treepl_result)) cat(paste(treepl_result, collapse = "\n"), "\n")
    stop("treePL final pass failed - output file missing or empty. Species tree dating is critical.")
  }

  calibrated_tree <- read.tree(final_out_file)
  cat("treePL dating succeeded\n")

} else if (calibration_method == "PATHd8") {

  # PATHd8 fallback path
  pathd8_path <- paste0(
    Sys.getenv("CONDA_PREFIX"), "/bin:",
    Sys.getenv("CONDA_PREFIX"), "/lib64:",
    Sys.getenv("LD_LIBRARY_PATH")
  )
  Sys.setenv(LD_LIBRARY_PATH = pathd8_path)

  pathd8_input <- "spp_tree_pathd8_input.txt"
  tree_string <- write.tree(species_tree)

  pathd8_lines <- c(tree_string)
  pathd8_lines <- c(pathd8_lines, "Sequence length = 10000;")
  for (i in seq_len(nrow(calibrations))) {
    pathd8_lines <- c(pathd8_lines,
      paste0("mrca: ", calibrations$tipA[i], ", ", calibrations$tipB[i],
             ", fixage=", calibrations$age_mya[i], ";")
    )
  }
  writeLines(pathd8_lines, pathd8_input)

  pathd8_output <- "spp_tree_pathd8_output.txt"
  cat("Running PATHd8...\n")
  pathd8_result <- tryCatch({
    system2("/usr/local/bin/PATHd8", args = c("-i", pathd8_input, "-o", pathd8_output),
            stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    cat("PATHd8 failed:", e$message, "\n")
    NULL
  })

  if (!file.exists(pathd8_output)) {
    stop("PATHd8 failed - output file not found. Species tree dating is critical.")
  }

  pathd8_out_lines <- readLines(pathd8_output, warn = FALSE)
  d8_line <- grep("^d8 tree", pathd8_out_lines, value = TRUE)
  if (length(d8_line) == 0) {
    stop("Could not parse PATHd8 output - no d8 tree line found.")
  }
  d8_tree_str <- sub("^d8 tree *: *", "", d8_line[1])
  calibrated_tree <- read.tree(text = d8_tree_str)
  cat("PATHd8 dating succeeded\n")
}

# ============================================================================
# Post-processing
# ============================================================================

cat("Time-calibration complete\n")
cat("Calibrated tree has", length(calibrated_tree$tip.label), "species\n")

# Force ultrametricity to handle rounding errors if needed
if (!is.ultrametric(calibrated_tree)) {
  if (is.ultrametric(calibrated_tree, tol = 1)) {
    cat("Output tree not strictly ultrametric - forcing ultrametricity to handle rounding errors\n")
    calibrated_tree <- force.ultrametric(calibrated_tree, method = "extend")
    if (!is.ultrametric(calibrated_tree)) {
      warning("Output tree cannot be made ultrametric - may have structural issues")
    }
  } else {
    warning("Output tree is far from ultrametric - may have structural issues")
  }
} else {
  cat("Output tree is ultrametric (time-calibrated)\n")
}

# Write the calibrated tree
cat("Writing time-calibrated tree to:", output_tree_path, "\n")
write.tree(calibrated_tree, file = output_tree_path)

cat("Done!\n")
