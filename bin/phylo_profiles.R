#!/usr/bin/env Rscript
library(data.table)
library(plyr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)
event_counts_file <- args[1]
species_event_counts_file <- args[2]
species_coverage_file <- args[3]
og <- args[4]
orthogroup_dir <- args[5]

get_per_spp_og_counts <-
  function(orthogroup_dir){
    og_counts <-
      read.table(paste0(orthogroup_dir, "/Orthogroups/Orthogroups.GeneCount.tsv"),
                 header = T, check.names = F)
    colnames(og_counts) <- gsub("\\..*", "", colnames(og_counts))

    # Calculate the number of species in each gene family
    og_counts$NumSpecies <-
      rowSums(og_counts[,-c(1,ncol(og_counts))] > 0)

    return(og_counts)
  }

get_og_event_counts <-
  function(per_spp_og_counts, event_counts_file, og){
    per_og_event_counts <-
      data.frame(
        gene_family = NA,
        speciation = NA,
        speciation_loss = NA,
        duplication = NA,
        loss = NA,
        number_gene_copies = NA,
        number_species = NA)

    tmp <- read.table(event_counts_file, sep = ":", check.names = F)

    species <-
      colnames(per_spp_og_counts)[-c(1, c((ncol(per_spp_og_counts)-1):ncol(per_spp_og_counts)))]

    counts <-
      per_spp_og_counts[which(per_spp_og_counts$Orthogroup == og),]

    per_og_event_counts[1,] <- c(og, tmp$V2, counts$NumSpecies)

    return(per_og_event_counts)
  }

get_og_events_per_spp <-
  function(per_spp_og_counts, species_event_counts_file, og){
    species <-
      colnames(per_spp_og_counts)[-c(1, c((ncol(per_spp_og_counts)-1):ncol(per_spp_og_counts)))]

    tmp <- read.table(species_event_counts_file, row.names = 1, check.names = F, sep = ",")
    tmp <- data.frame(t(tmp[which(rownames(tmp) %in% species),]), check.names = F)

    gf_col <- data.frame(gene_family = og)

    per_spp_og_speciation <-
      data.frame(matrix(ncol = length(species)+1, nrow = 0))
    colnames(per_spp_og_speciation) <- c("gene_family", species)
    per_spp_og_duplication <-
      data.frame(matrix(ncol = length(species)+1, nrow = 0))
    colnames(per_spp_og_duplication) <- c("gene_family", species)
    per_spp_og_loss <-
      data.frame(matrix(ncol = length(species)+1, nrow = 0))
    colnames(per_spp_og_loss) <- c("gene_family", species)

    # Use drop=FALSE to preserve column names when there's only 1 species
    specs <- cbind(gf_col, tmp[1, , drop=FALSE])
    dups <- cbind(gf_col, tmp[2, , drop=FALSE])
    loss <- cbind(gf_col, tmp[3, , drop=FALSE])

    per_spp_og_speciation <-
      plyr::rbind.fill(per_spp_og_speciation, specs)
    per_spp_og_duplication <-
      plyr::rbind.fill(per_spp_og_duplication, dups)
    per_spp_og_loss <-
      plyr::rbind.fill(per_spp_og_loss, loss)

    return(list(speciations = per_spp_og_speciation,
                duplications = per_spp_og_duplication,
                losses = per_spp_og_loss))
  }

# Load OG counts table
per_spp_og_counts <- get_per_spp_og_counts(orthogroup_dir)

# Get per-OG event counts
per_og_event_res <- get_og_event_counts(per_spp_og_counts, event_counts_file, og)

# Get per-species event counts
per_spp_events <- get_og_events_per_spp(per_spp_og_counts, species_event_counts_file, og)

# Write outputs with OG-specific filenames
write.table(per_spp_events$speciations,
            file = paste0(og, "_speciation_count.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(per_spp_events$duplications,
            file = paste0(og, "_duplication_count.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(per_spp_events$losses,
            file = paste0(og, "_loss_count.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
