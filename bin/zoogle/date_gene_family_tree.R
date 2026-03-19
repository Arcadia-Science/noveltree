#!/usr/bin/env Rscript

# date_gene_family_tree.R
# Time-calibrate gene family trees using reconciliation-filtered calibrations.
# Only speciation (S, SL) nodes from GeneRax NHX trees are used as calibration
# points. Duplication (D) and transfer (T, TL) nodes are excluded.
#
# Usage:
#   Rscript date_gene_family_tree.R <reconciled_tree> <species_tree> \
#     <alignment> <og_name> <max_treepl_tips> \
#     <out_dated_tree> <out_calibrations_csv>

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phangorn)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript date_gene_family_tree.R <reconciled_tree> ",
       "<species_tree> <alignment> <og_name> <max_treepl_tips> ",
       "<out_dated_tree> <out_calibrations_csv>")
}

tree_path       <- args[1]
spp_tree_path   <- args[2]
alignment_path  <- args[3]
og_name         <- args[4]
max_treepl_tips <- as.integer(args[5])
out_tree_path   <- args[6]
out_csv_path    <- args[7]

cat("=== date_gene_family_tree.R ===\n")
cat("OG:", og_name, "\n")
cat("Reconciled tree:", tree_path, "\n")
cat("Species tree:", spp_tree_path, "\n")
cat("Alignment:", alignment_path, "\n")
cat("Max treePL tips:", max_treepl_tips, "\n")

# ============================================================================
# Step 0: Read inputs
# ============================================================================

# Read the reconciled gene family tree — has both ML branch lengths and
# S/D/T internal node labels from GeneRax reconciliation
gf_tree <- read.tree(tree_path)

# Read the time-calibrated species tree
spp_tree <- read.tree(spp_tree_path)

# ============================================================================
# Step 1: Parse event annotations from reconciled tree node labels
# ============================================================================

# The reconciled tree includes internal node labels indicating event types:
# S (speciation), D (duplication), T@donor@recipient (transfer).
# ape::read.tree() parses these directly as $node.label.

# Extract per-internal-node event types from node labels
# Node labels may include transfer info like "T@donor@recipient"
node_events <- gf_tree$node.label
node_events <- sub("@.*", "", node_events)  # strip transfer details

n_tips_ev <- length(gf_tree$tip.label)
n_internal <- gf_tree$Nnode

# Classify: S -> speciation (usable for calibration); D, T -> excluded
is_speciation <- node_events == "S"
n_spec <- sum(is_speciation, na.rm = TRUE)
n_dup <- sum(node_events == "D", na.rm = TRUE)
n_trans <- sum(node_events == "T", na.rm = TRUE)
cat("  Events found: S =", n_spec, ", D =", n_dup,
    ", T =", n_trans, "\n")

# ============================================================================
# Step 2: Extract speciation-only calibrations
# ============================================================================

# Collect descendant tips reachable only through S (speciation) and D
# (duplication) events.  Traversal stops at T (transfer) nodes so that
# species introduced by horizontal transfer do not inflate the MRCA age
# used for calibration.
get_native_tips <- function(tree, start_node, node_events, n_tips) {
  tips  <- character(0)
  queue <- start_node
  while (length(queue) > 0) {
    node  <- queue[1]
    queue <- queue[-1]
    if (node <= n_tips) {
      tips <- c(tips, tree$tip.label[node])
    } else {
      ev <- node_events[node - n_tips]
      # Stop traversal at transfer nodes (but always enter start_node itself)
      if (!is.na(ev) && ev == "T" && node != start_node) next
      children <- tree$edge[tree$edge[, 1] == node, 2]
      queue <- c(queue, children)
    }
  }
  tips
}

extract_speciation_calibrations <- function(gf_tree, spp_tree,
                                            is_speciation, node_events) {
  n_tips <- length(gf_tree$tip.label)

  # Get species tree node depths for age calculation
  spp_depths <- node.depth.edgelength(spp_tree)
  spp_root_depth <- max(spp_depths)

  calibrations <- data.frame(
    gf_mrca = integer(0),
    age_mya = numeric(0),
    tipA = character(0),
    tipB = character(0),
    n_desc_tips = integer(0),
    events_node = integer(0),
    stringsAsFactors = FALSE
  )

  spec_nodes <- which(is_speciation)
  if (length(spec_nodes) == 0) return(calibrations)

  for (idx in spec_nodes) {
    ev_node <- n_tips + idx

    # Get descendant tips reachable through S/D only (stop at T nodes)
    native_tips <- get_native_tips(gf_tree, ev_node, node_events, n_tips)
    if (length(native_tips) < 2) next

    # Map to species names (remove protein ID after last underscore)
    desc_species <- unique(sub("_[^_]+$", "", native_tips))

    # Find species shared with the species tree
    shared_spp <- intersect(desc_species, spp_tree$tip.label)
    if (length(shared_spp) < 2) next

    # Get MRCA age from species tree
    spp_mrca <- getMRCA(spp_tree, shared_spp)
    if (is.null(spp_mrca)) next
    age_mya <- spp_root_depth - spp_depths[spp_mrca]

    if (is.na(age_mya) || age_mya <= 0) next

    # The calibration only makes sense if the S node's two children each
    # contribute at least one native tip (otherwise, all tips are from one
    # daughter lineage and the calibration would duplicate a descendant node).
    children <- gf_tree$edge[gf_tree$edge[, 1] == ev_node, 2]
    left_tips  <- intersect(native_tips, gf_tree$tip.label[unlist(
      Descendants(gf_tree, children[1], type = "tips"))])
    right_tips <- intersect(native_tips, gf_tree$tip.label[unlist(
      Descendants(gf_tree, children[2], type = "tips"))])
    if (length(left_tips) == 0 || length(right_tips) == 0) next

    # gf_mrca must be ev_node since tips span both children
    gf_mrca <- ev_node

    # Pick one representative tip from each child for treePL mrca specification
    tipA <- left_tips[1]
    tipB <- right_tips[1]

    calibrations <- rbind(calibrations, data.frame(
      gf_mrca = gf_mrca,
      age_mya = age_mya,
      tipA = tipA,
      tipB = tipB,
      n_desc_tips = length(native_tips),
      events_node = ev_node,
      stringsAsFactors = FALSE
    ))
  }

  calibrations
}

calibrations <- extract_speciation_calibrations(
  gf_tree, spp_tree, is_speciation, node_events
)
cat("  Raw speciation calibrations:", nrow(calibrations), "\n")

# ============================================================================
# Step 3: Deduplicate and resolve conflicts
# ============================================================================

if (nrow(calibrations) > 0) {
  # One calibration per gene-tree MRCA node (keep entry with most descendant tips)
  calibrations <- calibrations[order(-calibrations$n_desc_tips), ]
  calibrations <- calibrations[!duplicated(calibrations$gf_mrca), ]

  # Remove parent-child age conflicts:
  # Sort by node depth (root-to-tip), drop calibrations where child age > parent age
  gf_depths <- node.depth.edgelength(gf_tree)
  calibrations$node_depth <- gf_depths[calibrations$gf_mrca]
  calibrations <- calibrations[order(calibrations$node_depth), ]

  dropped <- c()
  if (nrow(calibrations) > 1) {
    for (i in 2:nrow(calibrations)) {
      # Check if any ancestor of this node has a calibration with younger age
      current_node <- calibrations$gf_mrca[i]
      current_age <- calibrations$age_mya[i]

      for (j in 1:(i - 1)) {
        if (j %in% dropped) next
        ancestor_node <- calibrations$gf_mrca[j]
        ancestor_age <- calibrations$age_mya[j]

        # Check if j is an ancestor of i
        if (ancestor_node %in% Ancestors(gf_tree, current_node, type = "all")) {
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
  }
  calibrations$node_depth <- NULL
}

cat("  Final calibrations after dedup/conflict resolution:", nrow(calibrations), "\n")

# ============================================================================
# Step 3b: Ancestor-descendant spacing thinning
# ============================================================================
# On any root-to-tip path, adjacent calibrations must be >= 5% of root age
# apart.  When two are closer, the shallower (younger) one is dropped.
# This prevents over-constraining treePL in densely-calibrated clades.

if (nrow(calibrations) > 1) {
  min_age_gap <- max(calibrations$age_mya) * 0.05
  # Sort oldest-first so we preferentially keep deeper calibrations
  calibrations <- calibrations[order(-calibrations$age_mya), ]

  drop_idx <- c()
  for (i in 1:(nrow(calibrations) - 1)) {
    if (i %in% drop_idx) next
    for (j in (i + 1):nrow(calibrations)) {
      if (j %in% drop_idx) next
      node_i <- calibrations$gf_mrca[i]
      node_j <- calibrations$gf_mrca[j]
      # Check if i is ancestor of j (i is older, j is younger)
      if (node_i %in% Ancestors(gf_tree, node_j, type = "all")) {
        age_diff <- calibrations$age_mya[i] - calibrations$age_mya[j]
        if (age_diff < min_age_gap) {
          drop_idx <- c(drop_idx, j)
        }
      }
    }
  }

  if (length(drop_idx) > 0) {
    cat("  Spacing thinning: dropped", length(drop_idx), "of",
        nrow(calibrations), "calibrations (min gap =",
        round(min_age_gap, 2), "Mya)\n")
    calibrations <- calibrations[-drop_idx, ]
  }
}

cat("  Calibrations after spacing thinning:", nrow(calibrations), "\n")

# ============================================================================
# Step 4: Set ±10% calibration brackets
# ============================================================================
# Gene tree calibrations use ±10% brackets around species tree ages,
# giving treePL freedom for lineage-specific rate variation.

if (nrow(calibrations) > 0) {
  calibrations$min_mya <- calibrations$age_mya * 0.90
  calibrations$max_mya <- calibrations$age_mya * 1.10
}

# ============================================================================
# Step 5: Count alignment columns (numsites for treePL/PATHd8)
# ============================================================================

count_alignment_columns <- function(fasta_path) {
  lines <- readLines(fasta_path, warn = FALSE)
  seq_lines <- lines[!grepl("^>", lines)]
  if (length(seq_lines) == 0) return(0)
  nchar(paste(seq_lines[1:min(length(seq_lines), 1)], collapse = ""))
}

numsites <- count_alignment_columns(alignment_path)
cat("  Alignment columns (numsites):", numsites, "\n")

# ============================================================================
# Step 6: Run dating
# ============================================================================

n_cal <- nrow(calibrations)
n_tips <- length(gf_tree$tip.label)
cat("  Number of calibrations:", n_cal, "\n")
cat("  Number of tips:", n_tips, "\n")

if (n_cal >= 2) {
  if (n_tips < max_treepl_tips) {
    # --- treePL ---
    cat("  Using treePL for dating...\n")

    # Write gene tree for treePL
    tree_file <- paste0(og_name, "_for_treepl.newick")
    write.tree(gf_tree, file = tree_file)

    # Build calibration lines
    cal_lines <- c()
    for (i in 1:n_cal) {
      cal_name <- paste0("cal", i)
      cal_lines <- c(cal_lines,
        paste("mrca =", cal_name, calibrations$tipA[i], calibrations$tipB[i]),
        paste("min =", cal_name, calibrations$min_mya[i]),
        paste("max =", cal_name, calibrations$max_mya[i])
      )
    }

    # CV pass to find optimal smoothing
    cv_config_file <- paste0(og_name, "_treepl_cv.config")
    cv_out_file <- paste0(og_name, "_treepl_cv_out.newick")
    cv_config <- c(
      paste("treefile =", tree_file),
      paste("numsites =", numsites),
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

    cv_output <- tryCatch({
      system2("/opt/conda/bin/treePL", args = cv_config_file,
              stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      cat("  treePL CV pass failed:", e$message, "\n")
      NULL
    })

    # Parse optimal smoothing from CV output
    optimal_smooth <- 100  # default
    if (!is.null(cv_output)) {
      smooth_line <- grep("Optimal smoothing value", cv_output, value = TRUE)
      if (length(smooth_line) > 0) {
        smooth_val <- as.numeric(sub(".*: *", "", smooth_line[1]))
        if (!is.na(smooth_val) && smooth_val > 0) {
          optimal_smooth <- smooth_val
          cat("  Optimal smoothing from CV:", optimal_smooth, "\n")
        }
      }
    }

    # Final pass with optimal smoothing
    final_config_file <- paste0(og_name, "_treepl_final.config")
    final_out_file <- paste0(og_name, "_treepl_dated.newick")
    final_config <- c(
      paste("treefile =", tree_file),
      paste("numsites =", numsites),
      paste("smooth =", optimal_smooth),
      cal_lines,
      paste("outfile =", final_out_file),
      "opt = 1",
      "optad = 1"
    )
    writeLines(final_config, final_config_file)

    treepl_result <- tryCatch({
      system2("/opt/conda/bin/treePL", args = final_config_file,
              stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      cat("  treePL final pass failed:", e$message, "\n")
      NULL
    })

    if (file.exists(final_out_file)) {
      dated_tree <- read.tree(final_out_file)
      cat("  treePL dating succeeded\n")
    } else {
      cat("  WARNING: treePL output not found, falling back to undated tree\n")
      dated_tree <- gf_tree
    }

  } else {
    # --- PATHd8 ---
    cat("  Using PATHd8 for dating (tree has", n_tips, "tips)...\n")

    # Write PATHd8 input
    pathd8_input <- paste0(og_name, "_pathd8_input.txt")
    tree_string <- write.tree(gf_tree)

    pathd8_lines <- c(tree_string)
    pathd8_lines <- c(pathd8_lines, paste("Sequence length =", numsites, ";"))
    for (i in 1:n_cal) {
      midpoint <- (calibrations$min_mya[i] + calibrations$max_mya[i]) / 2
      pathd8_lines <- c(pathd8_lines,
        paste0("mrca: ", calibrations$tipA[i], ", ", calibrations$tipB[i],
               ", fixage=", midpoint, ";")
      )
    }
    writeLines(pathd8_lines, pathd8_input)

    pathd8_output <- paste0(og_name, "_pathd8_output.txt")
    pathd8_result <- tryCatch({
      system2("/usr/local/bin/PATHd8", args = c("-i", pathd8_input,
                                                 "-o", pathd8_output),
              stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      cat("  PATHd8 failed:", e$message, "\n")
      NULL
    })

    if (file.exists(pathd8_output)) {
      # Parse PATHd8 output - extract the d8 dated tree
      pathd8_out_lines <- readLines(pathd8_output, warn = FALSE)
      d8_line <- grep("^d8 tree", pathd8_out_lines, value = TRUE)
      if (length(d8_line) > 0) {
        d8_tree_str <- sub("^d8 tree *: *", "", d8_line[1])
        dated_tree <- read.tree(text = d8_tree_str)
        cat("  PATHd8 dating succeeded\n")
      } else {
        cat("  WARNING: Could not parse PATHd8 output, using undated tree\n")
        dated_tree <- gf_tree
      }
    } else {
      cat("  WARNING: PATHd8 output not found, using undated tree\n")
      dated_tree <- gf_tree
    }
  }

} else if (n_cal == 1) {
  # --- Strict molecular clock with single calibration ---
  cat("  Using strict molecular clock (1 calibration)...\n")

  # Fix the calibration node age and use chronos with strict clock
  cal_node <- calibrations$gf_mrca[1]
  cal_age <- calibrations$age_mya[1]

  dated_tree <- tryCatch({
    calib <- makeChronosCalib(gf_tree, node = cal_node,
                               age.min = cal_age * 0.90, age.max = cal_age * 1.10)
    chronos(gf_tree, model = "strict", calibration = calib)
  }, error = function(e) {
    cat("  chronos() failed:", e$message, "\n")
    cat("  Returning undated tree\n")
    gf_tree
  })

} else {
  # --- No calibrations: output undated tree ---
  cat("  WARNING: No speciation calibrations found. Outputting undated tree.\n")
  cat("  Tree will still be usable for GLS but without time calibration.\n")
  dated_tree <- gf_tree
}

# ============================================================================
# Step 7: Post-processing
# ============================================================================

# Force ultrametricity if the tree was dated but isn't perfectly ultrametric
if (n_cal >= 1 && !is.null(dated_tree$edge.length)) {
  if (!is.ultrametric(dated_tree, tol = 0.01) && is.ultrametric(dated_tree, tol = 1)) {
    cat("  Forcing ultrametricity (minor adjustments)...\n")
    dated_tree <- force.ultrametric(dated_tree, method = "extend")
  }
}

# Write outputs
write.tree(dated_tree, file = out_tree_path)
cat("  Dated tree written to:", out_tree_path, "\n")

# Write calibrations CSV
if (nrow(calibrations) > 0) {
  write.csv(calibrations, file = out_csv_path, row.names = FALSE)
} else {
  # Write header-only CSV
  write.csv(data.frame(
    gf_mrca = integer(0), age_mya = numeric(0),
    tipA = character(0), tipB = character(0),
    n_desc_tips = integer(0), events_node = integer(0),
    min_mya = numeric(0), max_mya = numeric(0)
  ), file = out_csv_path, row.names = FALSE)
}
cat("  Calibrations CSV written to:", out_csv_path, "\n")
cat("=== Done ===\n")
