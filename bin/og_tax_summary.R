#!/usr/bin/env Rscript

# TODO - BUILD IN OPTPARSE TO STREAMLINE THE FOLLOWING SECTION
# Read in commandline arguments
args = commandArgs(trailingOnly=TRUE)

# Get path to the orthogroup gene count summary file
ogcounts <- args[1]
# And then the path to the samplesheet
samples <- args[2]

num_seq_filt <- as.numeric(args[3])
num_spp_filt <- as.numeric(args[4])
prop_spp_spptree_filt <- as.numeric(args[5])
num_grp_filt <- as.numeric(args[6])
copy_num_filt1 <- as.numeric(args[7])
copy_num_filt2 <- as.numeric(args[8])

ogs <- read.delim(ogcounts, check.names = FALSE)
samples <- read.delim(samples, sep = ",")

# Convert the % species filter to the number of species required
num_spp_spptree_filt <- round(nrow(samples) * prop_spp_spptree_filt)

colnames(ogs) <- gsub("\\..*", "", colnames(ogs))

# Extract and remove the totals column
totals <- ogs['Total']
ogs <- ogs[,-which(colnames(ogs) == 'Total')]

# Convert zero-counts to NA
ogs[-1][ogs[-1] == 0] <- NA

# Calculate the mean copy number per species
copynum <- apply(ogs[-1], 1, mean, na.rm=T)

# Convert counts to a binary presence/absence
ogs[-1][ogs[-1] > 0] <- 1

# For each orthogroup, count the number of species included
numspp <- rowSums(ogs[-1], na.rm = T)

# Convert species names to taxon group
for(i in 2:ncol(ogs)){
    # Identify the species
    spp <- colnames(ogs)[i]
    
    # Strip the EukProt ID (if relevant)
    spp <- sub("EP0.*?_", "", spp)

    # And the taxonomic group for this species
    grp <- unique(as.character(samples$taxonomy[which(samples$species == spp)]))

    # Now replace the species name with group name
    colnames(ogs)[i] <- grp
}

# Count number of taxonomic groups included in each orthogroup
taxcount <- t(rowsum(t(ogs[-1]),
              group = colnames(ogs)[-1],
              na.rm = T))
# Convert counts to binary for easy counts of species in each og.
taxcount[taxcount > 1] <- 1
taxcount <- rowSums(taxcount)

res <-
    data.frame(
        orthogroup = ogs$Orthogroup,
        num_spp = numspp,
        total_copy_num = totals[,1],
        mean_copy_num = copynum,
        num_tax_grps = taxcount
    )

# Create the subsets
spptree_core <-
    res[which(res$total_copy_num >= num_seq_filt &
              res$mean_copy_num <= copy_num_filt1 &
              res$num_spp >= num_spp_filt &
              res$num_spp >= num_spp_spptree_filt &
              res$num_tax_grps >= num_grp_filt),]
genetree_core <-
    res[which(res$total_copy_num >= num_seq_filt &
              res$mean_copy_num <= copy_num_filt2 &
              res$num_spp >= num_spp_filt &
              res$num_tax_grps >= num_grp_filt),]

# And remove the species tree core ogs from the remnants we'll just infer gene
# family trees for (to eliminate redundant computational effort)
genetree_core <- genetree_core[-which(genetree_core$orthogroup %in% spptree_core$orthogroup),]

write.csv(res, paste0("all_ogs_counts.csv"), quote = F, row.names = F)
write.csv(spptree_core, paste0("spptree_core_ogs_counts.csv"), quote = F, row.names = F)
write.csv(genetree_core, paste0("genetree_core_ogs_counts.csv"), quote = F, row.names = F)
