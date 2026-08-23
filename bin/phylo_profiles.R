#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: phylo_profiles.R <species_event_counts> <species_coverage> <orthogroup>")
}

species_event_counts_file <- args[1]
species_coverage_file <- args[2]
og <- args[3]

read_species_coverage <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2 || trimws(lines[1]) != "SPECIES: FAMILY_COVERAGE") {
    stop(paste("Invalid GeneRax per-species coverage file:", path))
  }
  fields <- strsplit(lines[-1], ":", fixed = TRUE)
  if (any(lengths(fields) != 2)) {
    stop(paste("Invalid GeneRax coverage row in:", path))
  }
  species <- trimws(vapply(fields, `[[`, character(1), 1))
  coverage <- suppressWarnings(as.numeric(trimws(vapply(fields, `[[`, character(1), 2))))
  if (length(species) == 0 || any(!nzchar(species)) || anyDuplicated(species)) {
    stop(paste("Coverage file has no species or duplicate species:", path))
  }
  if (any(!is.finite(coverage)) || any(coverage < 0)) {
    stop(paste("Coverage file has invalid coverage values:", path))
  }
  data.frame(species = species, coverage = coverage, stringsAsFactors = FALSE)
}

build_profile <- function(events, coverage, column, og) {
  species <- coverage$species
  positions <- match(species, events$species_label)
  missing_covered <- is.na(positions) & coverage$coverage > 0
  if (any(missing_covered)) {
    missing <- paste(species[missing_covered], collapse = ", ")
    stop(paste("Covered species absent from GeneRax event table:", missing))
  }

  # GeneRax's species event table is sparse: zero-coverage terminal species
  # with no inferred events are omitted entirely. Preserve explicit event rows
  # (including losses in zero-coverage species), and zero-fill only omitted,
  # zero-coverage species.
  values <- numeric(length(species))
  represented <- !is.na(positions)
  represented_values <- events[[column]][positions[represented]]
  if (any(!is.finite(represented_values)) || any(represented_values < 0)) {
    stop(paste("Invalid GeneRax event values in column:", column))
  }
  values[represented] <- represented_values

  result <- data.frame(gene_family = og, check.names = FALSE)
  for (index in seq_along(species)) {
    result[[species[index]]] <- values[index]
  }
  result
}

coverage <- read_species_coverage(species_coverage_file)
events <- read.csv(
  species_event_counts_file,
  check.names = FALSE,
  strip.white = TRUE,
  stringsAsFactors = FALSE
)
required <- c("species_label", "speciations", "duplications", "losses")
missing_columns <- setdiff(required, colnames(events))
if (length(missing_columns) > 0) {
  stop(paste("Missing GeneRax event columns:", paste(missing_columns, collapse = ", ")))
}
if (anyDuplicated(events$species_label)) {
  stop("GeneRax event table contains duplicate species labels")
}
terminal_species <- events$species_label[!startsWith(events$species_label, "Node_")]
unexpected_species <- setdiff(terminal_species, coverage$species)
if (length(unexpected_species) > 0) {
  stop(paste(
    "GeneRax event table contains species absent from coverage:",
    paste(unexpected_species, collapse = ", ")
  ))
}

outputs <- list(
  speciation = build_profile(events, coverage, "speciations", og),
  duplication = build_profile(events, coverage, "duplications", og),
  loss = build_profile(events, coverage, "losses", og)
)
for (kind in names(outputs)) {
  write.table(
    outputs[[kind]],
    file = paste0(og, "_", kind, "_count.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}
