#!/usr/bin/env Rscript
library(data.table)
library(plyr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)
event_counts_file <- args[1]
species_event_counts_file <- args[2]
transfer_event_counts_file <- args[3]
species_coverage_file <- args[4]
og <- args[5]
orthogroup_dir <- args[6]

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
        transfer = NA,
        transfer_loss = NA,
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

get_tranfer_donor_recips <-
  function(per_spp_og_counts, transfer_event_counts_file, og){
    species <-
      colnames(per_spp_og_counts)[-c(1,(ncol(per_spp_og_counts)-1):ncol(per_spp_og_counts))]

    per_spp_og_transfer_donor <-
      data.frame(matrix(ncol = length(species)+1, nrow = 0))
    colnames(per_spp_og_transfer_donor) <- c("gene_family", species)
    per_spp_og_transfer_recip <-
      data.frame(matrix(ncol = length(species)+1, nrow = 0))
    colnames(per_spp_og_transfer_recip) <- c("gene_family", species)

    transf_count_mat <-
      matrix(nrow = length(species),
             ncol = length(species),
             dimnames = list(species, species),
             data = 0)

    # Species present in this gene family
    og_spps <- per_spp_og_counts[which(per_spp_og_counts$Orthogroup == og),
                                 -c(1,(ncol(per_spp_og_counts)-1):ncol(per_spp_og_counts))]
    og_spps <- colnames(og_spps[which(og_spps[1,] > 0),])

    # Only process if transfer events were inferred
    if(file.size(transfer_event_counts_file) != 0L){
      tmp <- read.table(transfer_event_counts_file, check.names = F)
      donor_spps <- og_spps[which(og_spps %in% tmp$V1)]
      donors <- summary(as.factor(tmp$V1))
      donors <- data.frame(as.list(donors[which(names(donors) %in% species)]), check.names = F)
      recip_spps <- og_spps[which(og_spps %in% tmp$V2)]
      recips <- summary(as.factor(tmp$V2))
      recips <- data.frame(as.list(recips[which(names(recips) %in% species)]), check.names = F)

      if(nrow(donors) > 0){
        donors$gene_family <- og
        non_donors <- og_spps[-which(og_spps %in% donor_spps)]
        if(length(non_donors) > 0){
          non_donors <-
            data.frame(matrix(ncol = length(non_donors),
                              dimnames = list(NULL, non_donors),
                              data = 0), check.names = F)
          donors <- cbind(donors, non_donors)
        }
        per_spp_og_transfer_donor <-
          plyr::rbind.fill(per_spp_og_transfer_donor, donors)
      }else{
        non_donors <-
          data.frame(matrix(ncol = length(og_spps)+1,
                            dimnames = list(NULL, c("gene_family", og_spps)),
                            data = c(og, rep(0, length(og_spps)))),
                        check.names = F)
        per_spp_og_transfer_donor <-
          plyr::rbind.fill(per_spp_og_transfer_donor, non_donors)
      }
      if(nrow(recips) > 0){
        recips$gene_family <- og
        non_recips <- og_spps[-which(og_spps %in% recip_spps)]
        if(length(non_recips) > 0){
          non_recips <-
            data.frame(matrix(ncol = length(non_recips),
                              dimnames = list(NULL, non_recips),
                              data = 0), check.names = F)
          recips <- cbind(recips, non_recips)
        }
        per_spp_og_transfer_recip <-
          plyr::rbind.fill(per_spp_og_transfer_recip, recips)
      }else{
        non_recips <-
          data.frame(matrix(ncol = length(og_spps)+1,
                            dimnames = list(NULL, c("gene_family", og_spps)),
                            data = c(og, rep(0, length(og_spps)))),
                     check.names = F)
        per_spp_og_transfer_recip <-
          plyr::rbind.fill(per_spp_og_transfer_recip, non_recips)
      }
      # Fill in the donor-recipient matrix
      tmp <- tmp[which(tmp$V1 %in% species & tmp$V2 %in% species),]
      for(x in 1:nrow(tmp)){
        donor <- which(species == tmp$V1[x])
        recip <- which(species == tmp$V2[x])
        transf_count_mat[donor, recip] <-
          transf_count_mat[donor, recip] + 1
      }
    }else{
      per_spp_og_transfer_donor <-
        data.frame(matrix(ncol = length(og_spps)+1,
                          dimnames = list(NULL, c("gene_family", og_spps)),
                          data = c(og, rep(0, length(og_spps)))),
                   check.names = F)
      per_spp_og_transfer_recip <-
        data.frame(matrix(ncol = length(og_spps)+1,
                          dimnames = list(NULL, c("gene_family", og_spps)),
                          data = c(og, rep(0, length(og_spps)))),
                   check.names = F)
    }
    return(list("summed_matrix" = transf_count_mat,
                "gf_transfer_donors" = per_spp_og_transfer_donor,
                "gf_transfer_recips" = per_spp_og_transfer_recip))
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
    per_spp_og_transfer <-
      data.frame(matrix(ncol = length(species)+1, nrow = 0))
    colnames(per_spp_og_transfer) <- c("gene_family", species)
    per_spp_og_loss <-
      data.frame(matrix(ncol = length(species)+1, nrow = 0))
    colnames(per_spp_og_loss) <- c("gene_family", species)

    # Use drop=FALSE to preserve column names when there's only 1 species
    specs <- cbind(gf_col, tmp[1, , drop=FALSE])
    dups <- cbind(gf_col, tmp[2, , drop=FALSE])
    loss <- cbind(gf_col, tmp[3, , drop=FALSE])
    transf <- cbind(gf_col, tmp[4, , drop=FALSE])

    per_spp_og_speciation <-
      plyr::rbind.fill(per_spp_og_speciation, specs)
    per_spp_og_duplication <-
      plyr::rbind.fill(per_spp_og_duplication, dups)
    per_spp_og_loss <-
      plyr::rbind.fill(per_spp_og_loss, loss)
    per_spp_og_transfer <-
      plyr::rbind.fill(per_spp_og_transfer, transf)

    return(list(speciations = per_spp_og_speciation,
                duplications = per_spp_og_duplication,
                transfers = per_spp_og_transfer,
                losses = per_spp_og_loss))
  }

# Load OG counts table
per_spp_og_counts <- get_per_spp_og_counts(orthogroup_dir)

# Get per-OG event counts
per_og_event_res <- get_og_event_counts(per_spp_og_counts, event_counts_file, og)

# Get per-species event counts
per_spp_events <- get_og_events_per_spp(per_spp_og_counts, species_event_counts_file, og)

# Get transfer donor-recipient info
transf_res <- get_tranfer_donor_recips(per_spp_og_counts, transfer_event_counts_file, og)

# Write outputs with OG-specific filenames
write.table(transf_res$summed_matrix,
            file = paste0(og, "_hgt_summed_counts.tsv"),
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(per_spp_events$speciations,
            file = paste0(og, "_speciation_count.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(per_spp_events$duplications,
            file = paste0(og, "_duplication_count.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(per_spp_events$losses,
            file = paste0(og, "_loss_count.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(transf_res$gf_transfer_donors,
            file = paste0(og, "_transfer_donor_count.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(transf_res$gf_transfer_recips,
            file = paste0(og, "_transfer_recipient_count.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
