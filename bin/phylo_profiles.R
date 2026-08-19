#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: phylo_profiles.R <species_event_counts> <species_coverage> <orthogroup>")
}

species_event_counts_file <- args[1]
species_coverage_file <- args[2]
og <- args[3]

read_covered_species <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2 || trimws(lines[1]) != "SPECIES: FAMILY_COVERAGE") {
    stop(paste("Invalid GeneRax per-species coverage file:", path))
  }
  species <- trimws(sub(":.*$", "", lines[-1]))
  species <- species[nzchar(species)]
  if (length(species) == 0 || anyDuplicated(species)) {
    stop(paste("Coverage file has no species or duplicate species:", path))
  }
  species
}

build_profile <- function(events, species, column, og) {
  positions <- match(species, events$species_label)
  if (anyNA(positions)) {
    missing <- paste(species[is.na(positions)], collapse = ", ")
    stop(paste("Species absent from GeneRax event table:", missing))
  }
  values <- events[[column]][positions]
  result <- data.frame(gene_family = og, check.names = FALSE)
  for (index in seq_along(species)) {
    result[[species[index]]] <- values[index]
  }
  result
}

species <- read_covered_species(species_coverage_file)
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

outputs <- list(
  speciation = build_profile(events, species, "speciations", og),
  duplication = build_profile(events, species, "duplications", og),
  loss = build_profile(events, species, "losses", og)
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
