#!/usr/bin/env Rscript
# Merge outputs from batched PHYLO_PROFILES runs
# - Concatenates per-gene-family TSV files
# - Sums transfer matrices element-wise

library(data.table)

message("Starting merge of phylo profiles batches")

# Function to concatenate TSV files (skip header on subsequent files)
concat_tsvs <- function(input_dir, output_file) {
    files <- list.files(input_dir, full.names = TRUE)
    message(paste("  Concatenating", length(files), "files from", input_dir))

    all_data <- rbindlist(lapply(files, fread), fill = TRUE)
    fwrite(all_data, output_file, sep = "\t")
}

# Function to sum transfer matrices
sum_matrices <- function(input_dir, output_file) {
    files <- list.files(input_dir, full.names = TRUE)
    message(paste("  Summing", length(files), "matrices from", input_dir))

    # Read first matrix to get dimensions and row/col names
    first_mat <- as.matrix(read.table(files[1], header = TRUE, row.names = 1,
                                       check.names = FALSE, sep = "\t"))
    result <- first_mat

    # Sum remaining matrices
    if (length(files) > 1) {
        for (f in files[2:length(files)]) {
            mat <- as.matrix(read.table(f, header = TRUE, row.names = 1,
                                        check.names = FALSE, sep = "\t"))
            result <- result + mat
        }
    }

    write.table(result, output_file, sep = "\t", quote = FALSE,
                row.names = TRUE, col.names = NA)
}

# Concatenate per-gene-family files
message("Concatenating per-gene-family files...")
concat_tsvs("duplication_counts", "duplication_count_per_species_per_gene_family.tsv")
concat_tsvs("loss_counts", "loss_count_per_per_species_gene_family.tsv")
concat_tsvs("speciation_counts", "speciation_count_per_species_per_gene_family.tsv")
concat_tsvs("transfer_donor_counts", "transfer_donor_count_per_species_per_gene_family.tsv")
concat_tsvs("transfer_recipient_counts", "transfer_recipient_count_per_species_per_gene_family.tsv")

# Sum transfer matrices
message("Summing transfer matrices...")
sum_matrices("hgt_matrices", "hgt_summed_counts_recip_donor.tsv")

message("Merge complete")
