#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Get the species ID and the filepath to the accessions from the commandline
species <- args[1]
accession_fpath <- args[2]

# Read in
accessions <- read.table(accession_fpath)[, 1]

# Specify the columns we want returned
cogeqc_annots <-
  c("organism_name", "organism_id", "accession", "xref_interpro")

# Define the query, and then pull them down.
query <-
  list(accession = accessions)

res <-
  queryup::query_uniprot(query, columns = cogeqc_annots, show_progress = FALSE)

# Update the column names
colnames(res) <- cogeqc_annots

# And save out to file
write.table(res, paste0(species, "_cogeqc_annotations.tsv"),
            col.names = TRUE, row.names = FALSE,
            sep = "\t")
