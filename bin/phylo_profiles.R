#!/usr/bin/env Rscript
library(data.table)
library(parallel)
library(plyr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)
event_counts_files <- args[1]
species_event_counts_files <- args[2]
transfer_event_counts_files <- args[3]
species_coverage_files <- args[4]
ogs <- args[5]
orthogroup_dir <- args[6]

per_og_events <- unlist(strsplit(event_counts_files, " "))
per_spp_og_events <- unlist(strsplit(species_event_counts_files, " "))
spp_tranf_rates_fpaths <- unlist(strsplit(transfer_event_counts_files, " "))
ogs <- unlist(strsplit(ogs, " "))

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
  function(i, per_spp_og_counts = per_spp_og_counts, per_og_events = per_og_events, ogs = ogs){
    # Populate an empty dataframe
    # Create a dataframe to store the per-OG event counts
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

    # Read in the the orthogroup-wide (across spp) event counts
    tmp <- read.table(per_og_events[i], sep = ":", check.names = F)
    gf <- ogs[i]
    gf_col <- data.frame(gene_family = gf)

    # Get the per-species gene-count for this gene family
    species <-
      colnames(per_spp_og_counts)[-c(1, c((ncol(per_spp_og_counts)-1):ncol(per_spp_og_counts)))]

    counts <-
      per_spp_og_counts[which(per_spp_og_counts$Orthogroup == gf),]
    og_spps <- species[which(species %in% colnames(counts))]

    # Now fill
    per_og_event_counts[1,] <- c(gf, tmp$V2, counts$NumSpecies)

    # Return these event results as output
    return(per_og_event_counts)
  }

get_tranfer_donor_recips <-
  function(i, per_spp_og_counts = NULL,
           spp_tranf_rates_fpaths = NULL,
           ogs = NULL){
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

    gf <- ogs[i]
    # And the species that are in this gene family
    og_spps <- per_spp_og_counts[which(per_spp_og_counts$Orthogroup == gf),
                                 -c(1,(ncol(per_spp_og_counts)-1):ncol(per_spp_og_counts))]
    og_spps <- colnames(og_spps[which(og_spps[1,] > 0),])
    # And then the table of recipient-donor events
    # But only if transfer events were inferred
    if(file.size(spp_tranf_rates_fpaths[i]) != 0L){
      tmp <- read.table(spp_tranf_rates_fpaths[i], check.names = F)
      donor_spps <- og_spps[which(og_spps %in% tmp$V1)]
      donors <- summary(as.factor(tmp$V1))
      donors <- data.frame(as.list(donors[which(names(donors) %in% species)]), check.names = F)
      recip_spps <- og_spps[which(og_spps %in% tmp$V2)]
      recips <- summary(as.factor(tmp$V2))
      recips <- data.frame(as.list(recips[which(names(recips) %in% species)]), check.names = F)

      # And only if the transfer events occurred between tips
      # Make sure species included in the gene family have integer
      # counts, all other species NA
      if(nrow(donors) > 0){
        donors$gene_family <- gf
        non_donors <- og_spps[-which(og_spps %in% donor_spps)]
        # If there are species who did not donate gene copies via transfer:
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
                            data = c(gf, rep(0, length(og_spps)))),
                        check.names = F)
        per_spp_og_transfer_donor <-
          plyr::rbind.fill(per_spp_og_transfer_donor, non_donors)
      }
      if(nrow(recips) > 0){
        recips$gene_family <- gf
        non_recips <- og_spps[-which(og_spps %in% recip_spps)]
        # If there are species who did not receive gene copies via transfer:
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
                            data = c(gf, rep(0, length(og_spps)))),
                     check.names = F)
        per_spp_og_transfer_recip <-
          plyr::rbind.fill(per_spp_og_transfer_recip, non_recips)
      }
      # Now fill in the donor-recipient matrix
      tmp <- tmp[which(tmp$V1 %in% species & tmp$V2 %in% species),]
      for(x in 1:nrow(tmp)){
        donor <- which(species == tmp$V1[x])
        recip <- which(species == tmp$V2[x])

        # y-axis: recipient, x-axis: donor
        transf_count_mat[donor, recip] <-
          transf_count_mat[donor, recip] + 1
      }
    }else{
      per_spp_og_transfer_donor <-
        data.frame(matrix(ncol = length(og_spps)+1,
                          dimnames = list(NULL, c("gene_family", og_spps)),
                          data = c(gf, rep(0, length(og_spps)))),
                   check.names = F)
      per_spp_og_transfer_recip <-
        data.frame(matrix(ncol = length(og_spps)+1,
                          dimnames = list(NULL, c("gene_family", og_spps)),
                          data = c(gf, rep(0, length(og_spps)))),
                   check.names = F)
    }
    return(list("summed_matrix" = transf_count_mat,
                "gf_transfer_donors" = per_spp_og_transfer_donor,
                "gf_transfer_recips" = per_spp_og_transfer_recip))
  }

get_og_events_per_spp <-
  function(i, per_spp_og_counts = per_spp_og_counts,
           per_spp_og_events = per_spp_og_events,
           ogs = ogs){
    # Get the names of species
    species <-
      colnames(per_spp_og_counts)[-c(1, c((ncol(per_spp_og_counts)-1):ncol(per_spp_og_counts)))]

    # Read in table of per-species event counts for this gene family
    tmp <- read.table(per_spp_og_events[i], row.names = 1, check.names = F, sep = ",")
    tmp <- data.frame(t(tmp[which(rownames(tmp) %in% species),]), check.names = F)

    # Identify which gene family we"re dealing with
    gf <- ogs[i]
    gf_col <- data.frame(gene_family = gf)

    # Create one empty table for each parameter for the per-species counts per OG
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

    # And populate counts of duplication, transfer and loss
    # Use drop=FALSE to preserve column names when there's only 1 species
    specs <- cbind(gf_col, tmp[1, , drop=FALSE])
    dups <- cbind(gf_col, tmp[2, , drop=FALSE])
    loss <- cbind(gf_col, tmp[3, , drop=FALSE])
    transf <- cbind(gf_col, tmp[4, , drop=FALSE])

    # Now populate, allowing for species to not be present
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

# Now some quick helper functions to pull out speciations, duplications,
# transfers, and losses
get_speciations <-
  function(i, per_spp_events = per_spp_events){
    speciactions <- per_spp_events[[i]]$speciations; return(speciactions)}
get_duplications <-
  function(i, per_spp_events = per_spp_events){
    duplications <- per_spp_events[[i]]$duplications; return(duplications)}
get_losses <-
  function(i, per_spp_events = per_spp_events){
    losses <- per_spp_events[[i]]$losses; return(losses)}

summarize_generax_per_species <-
  function(per_og_events,
           per_spp_og_events,
           spp_tranf_rates_fpaths,
           ogs,
           per_spp_og_counts,
           nparallel = detectCores()-1){
    # Get the counts of each event type (duplications, transfers, losses, etc)
    # per-og, across all species
    message("Extracting event counts for each species per gene family.")
    per_og_event_res <-
      do.call(rbind, mclapply(X = 1:length(per_og_events),
                              get_og_event_counts, per_spp_og_counts = per_spp_og_counts,
                              per_og_events = per_og_events, ogs = ogs,
                              mc.cores = nparallel))

    # Get the counts of events per species, per orthogroup
    # Begin by first summarizing these event counts per orthogroup
    message("Extracting event counts per-species, per-orthogroup.")
    per_spp_events <-
      mclapply(1:length(per_spp_og_events),
               get_og_events_per_spp, per_spp_og_counts = per_spp_og_counts,
               per_spp_og_events = per_spp_og_events, ogs = ogs,
               mc.cores = nparallel)

    # And then pull out each event type individually
    message("Now, pulling out each event type individually.")
    per_spp_og_speciation <-
      do.call(rbind, mclapply(1:length(per_spp_events), get_speciations,
                              per_spp_events = per_spp_events,
                              mc.cores = nparallel))
    per_spp_og_duplication <-
      do.call(rbind, mclapply(1:length(per_spp_events), get_duplications,
                              per_spp_events = per_spp_events,
                              mc.cores = nparallel))
    per_spp_og_loss <-
      do.call(rbind, mclapply(1:length(per_spp_events), get_losses,
                              per_spp_events = per_spp_events,
                              mc.cores = nparallel))

    # Clean up the large interim list
    rm(per_spp_events)

    # Now, focusing on transfers - get a summed matrix of transfers among species,
    # with donors along the x-axis, and recipients along the y.
    # y-axis: recipient, x-axis: donor
    message("Summarizing gene transfer recipient events.")
    transf_res <-
      transpose(mclapply(1:length(spp_tranf_rates_fpaths),
                         get_tranfer_donor_recips,
                         per_spp_og_counts = per_spp_og_counts,
                         spp_tranf_rates_fpaths = spp_tranf_rates_fpaths,
                         ogs = ogs,
                         mc.cores = nparallel))
    message("Summarizing gene transfer events into a matrix of donor-recipient species pairs.")
    transf_count_mat <- Reduce("+", transf_res$summed_matrix)
    message("Pulling out the count of transfer-donor events for each species per gene family")
    transf_donors <- do.call("rbind", transf_res$gf_transfer_donors)
    message("Pulling out the count of transfer-recipient events for each species per gene family")
    transf_recips <- do.call("rbind", transf_res$gf_transfer_recips)

    # Generate a list containing alll required outputs for plotting
    results <- list(
      events_per_og = per_og_event_res,
      lgt_count_mat = transf_count_mat,
      speciations_per_spp =  per_spp_og_speciation,
      duplications_per_spp =  per_spp_og_duplication,
      losses_per_spp = per_spp_og_loss,
      transfer_donor_counts = transf_donors,
      transfer_recip_counts = transf_recips)

    out_dir = "."

    write.table(transf_count_mat,
                file = "hgt_summed_counts_recip_donor.tsv",
                sep = "\t", quote = F, row.names = T, col.names = NA)
    write.table(per_spp_og_speciation,
                file = "speciation_count_per_species_per_gene_family.tsv",
                sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(per_spp_og_duplication,
                file = "duplication_count_per_species_per_gene_family.tsv",
                sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(per_spp_og_loss,
                file = "loss_count_per_per_species_gene_family.tsv",
                sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(transf_donors,
                file = "transfer_donor_count_per_species_per_gene_family.tsv",
                sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(transf_recips,
                file = "transfer_recipient_count_per_species_per_gene_family.tsv",
                sep = "\t", quote = F, row.names = F, col.names = T)
  }

per_spp_og_counts <- get_per_spp_og_counts(orthogroup_dir)

generax_res_per_species <-
  summarize_generax_per_species(
    per_og_events,
    per_spp_og_events,
    spp_tranf_rates_fpaths,
    ogs,
    per_spp_og_counts
  )
