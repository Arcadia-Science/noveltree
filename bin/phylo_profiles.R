#!/usr/bin/env Rscript
message("DEBUG: Script starting"); flush(stderr())
library(data.table)
message("DEBUG: Loaded data.table"); flush(stderr())
library(parallel)
message("DEBUG: Loaded parallel"); flush(stderr())
library(plyr)
message("DEBUG: Loaded plyr"); flush(stderr())
library(purrr)
message("DEBUG: Loaded purrr"); flush(stderr())

args = commandArgs(trailingOnly=TRUE)
message("DEBUG: Got args"); flush(stderr())
event_counts_files <- args[1]
species_event_counts_files <- args[2]
transfer_event_counts_files <- args[3]
species_coverage_files <- args[4]
ogs <- args[5]
orthogroup_dir <- args[6]
message("DEBUG: Parsed args"); flush(stderr())

per_og_events <- unlist(strsplit(event_counts_files, " "))
message(paste("DEBUG: Split event_counts_files, length:", length(per_og_events))); flush(stderr())
per_spp_og_events <- unlist(strsplit(species_event_counts_files, " "))
message(paste("DEBUG: Split species_event_counts_files, length:", length(per_spp_og_events))); flush(stderr())
spp_tranf_rates_fpaths <- unlist(strsplit(transfer_event_counts_files, " "))
message(paste("DEBUG: Split transfer_event_counts_files, length:", length(spp_tranf_rates_fpaths))); flush(stderr())
ogs <- unlist(strsplit(ogs, " "))
message(paste("DEBUG: Split ogs, length:", length(ogs))); flush(stderr())

# Chunk size for memory-efficient processing
CHUNK_SIZE <- 100

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
           nparallel = detectCores()-1,
           chunk_size = CHUNK_SIZE){

    # Get species names for initializing the transfer matrix
    species <-
      colnames(per_spp_og_counts)[-c(1, (ncol(per_spp_og_counts)-1):ncol(per_spp_og_counts))]

    # Initialize the cumulative transfer count matrix (summed incrementally)
    transf_count_mat <- matrix(0, nrow = length(species), ncol = length(species),
                               dimnames = list(species, species))

    # Calculate number of chunks
    n_total <- length(ogs)
    n_chunks <- ceiling(n_total / chunk_size)
    message(paste("Processing", n_total, "orthogroups in", n_chunks, "chunks of", chunk_size))

    # Initialize lists to collect results (will use rbindlist at the end)
    all_og_event_res <- vector("list", n_chunks)
    all_speciations <- vector("list", n_chunks)
    all_duplications <- vector("list", n_chunks)
    all_losses <- vector("list", n_chunks)
    all_transf_donors <- vector("list", n_chunks)
    all_transf_recips <- vector("list", n_chunks)

    # Process in chunks
    for (chunk_idx in 1:n_chunks) {
      start_idx <- (chunk_idx - 1) * chunk_size + 1
      end_idx <- min(chunk_idx * chunk_size, n_total)
      chunk_indices <- start_idx:end_idx

      message(paste("Processing chunk", chunk_idx, "of", n_chunks,
                    "(orthogroups", start_idx, "to", end_idx, ")"))

      # Get the counts of each event type per-og for this chunk
      chunk_og_event_res <-
        rbindlist(mclapply(X = chunk_indices,
                           get_og_event_counts, per_spp_og_counts = per_spp_og_counts,
                           per_og_events = per_og_events, ogs = ogs,
                           mc.cores = nparallel))
      all_og_event_res[[chunk_idx]] <- chunk_og_event_res

      # Get the counts of events per species, per orthogroup for this chunk
      chunk_spp_events <-
        mclapply(chunk_indices,
                 get_og_events_per_spp, per_spp_og_counts = per_spp_og_counts,
                 per_spp_og_events = per_spp_og_events, ogs = ogs,
                 mc.cores = nparallel)

      # Extract each event type and store
      all_speciations[[chunk_idx]] <-
        rbindlist(lapply(chunk_spp_events, function(x) x$speciations), fill = TRUE)
      all_duplications[[chunk_idx]] <-
        rbindlist(lapply(chunk_spp_events, function(x) x$duplications), fill = TRUE)
      all_losses[[chunk_idx]] <-
        rbindlist(lapply(chunk_spp_events, function(x) x$losses), fill = TRUE)

      # Clean up chunk per-species events
      rm(chunk_spp_events)

      # Process transfers for this chunk - sum matrices incrementally
      message(paste("  Processing transfer events for chunk", chunk_idx))
      chunk_transf_res <-
        mclapply(chunk_indices,
                 get_tranfer_donor_recips,
                 per_spp_og_counts = per_spp_og_counts,
                 spp_tranf_rates_fpaths = spp_tranf_rates_fpaths,
                 ogs = ogs,
                 mc.cores = nparallel)

      # Incrementally sum transfer matrices (the key memory optimization)
      for (res in chunk_transf_res) {
        transf_count_mat <- transf_count_mat + res$summed_matrix
      }

      # Collect donor/recipient tables
      all_transf_donors[[chunk_idx]] <-
        rbindlist(lapply(chunk_transf_res, function(x) x$gf_transfer_donors), fill = TRUE)
      all_transf_recips[[chunk_idx]] <-
        rbindlist(lapply(chunk_transf_res, function(x) x$gf_transfer_recips), fill = TRUE)

      # Clean up chunk transfer results
      rm(chunk_transf_res)
      gc()
    }

    # Combine all chunks using rbindlist (memory efficient)
    message("Combining results from all chunks...")
    per_og_event_res <- rbindlist(all_og_event_res, fill = TRUE)
    per_spp_og_speciation <- rbindlist(all_speciations, fill = TRUE)
    per_spp_og_duplication <- rbindlist(all_duplications, fill = TRUE)
    per_spp_og_loss <- rbindlist(all_losses, fill = TRUE)
    transf_donors <- rbindlist(all_transf_donors, fill = TRUE)
    transf_recips <- rbindlist(all_transf_recips, fill = TRUE)

    # Clean up chunk lists
    rm(all_og_event_res, all_speciations, all_duplications, all_losses,
       all_transf_donors, all_transf_recips)
    gc()

    # Write outputs
    message("Writing output files...")
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

message("DEBUG: About to call get_per_spp_og_counts"); flush(stderr())
per_spp_og_counts <- get_per_spp_og_counts(orthogroup_dir)
message(paste("DEBUG: Loaded per_spp_og_counts, dim:", nrow(per_spp_og_counts), "x", ncol(per_spp_og_counts))); flush(stderr())

message("DEBUG: About to call summarize_generax_per_species"); flush(stderr())
generax_res_per_species <-
  summarize_generax_per_species(
    per_og_events,
    per_spp_og_events,
    spp_tranf_rates_fpaths,
    ogs,
    per_spp_og_counts
  )
