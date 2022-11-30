#!/usr/bin/env Rscript

# Read in commandline arguments
args = commandArgs(trailingOnly=TRUE)

# Get path to the orthogroup gene count summary file
ogcounts <- args[1]
# And then the path to the samplesheet
samples <- args[2]

numsppFilt <- args[3]
numgrpFilt <- args[4]
copynumFilt1 <- args[5]
copynumFilt2 <- args[6]

ogs <- read.delim(ogcounts)
samples <- read.delim(samples, sep = ",")

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
    
    # And the taxonomic group for this species
    grp <- unique(as.character(samples$taxonomy[which(samples$species == spp)]))
    
    # Now replace the species name with group name
    colnames(ogs)[i] <- grp
}

# Count number of taxonomic groups included in each orthogroup
taxcount <- t(rowsum(t(ogs[-1]), 
              group = colnames(ogs)[-1], 
              na.rm = T))
taxcount[taxcount > 1] <- 1
taxcount <- rowSums(taxcount)

res <- 
    data.frame(
        orthogroup = ogs$Orthogroup,
        num_spp = numspp,
        total_copy_num = totals,
        mean_copy_num = copynum,
        num_tax_grps = taxcount
    )

# Create the subsets
extreme_core <- 
    res[which(res$mean_copy_num <= copynumFilt1 & 
              res$num_spp <= numsppFilt & 
              res$num_tax_grps >= numgrpFilt),]
remaining_core <- 
    res[which(res$mean_copy_num <= copynumFilt2 & 
              res$num_spp <= numsppFilt & 
              res$num_tax_grps >= numgrpFilt),]

# And remove the extreme core ogs from the remnants (to eliminate redundant
# computational effort)
remaining_core <- remaining_core[-which(remaining_core$orthogroup %in% extreme_core$orthogroup),]

write.csv(res, paste0("all_ogs_counts.csv"), quote = F, row.names = F)
write.csv(extreme_core, paste0("extreme_core_ogs_counts.csv"), quote = F, row.names = F)
write.csv(remaining_core, paste0("remaining_core_ogs_counts.csv"), quote = F, row.names = F)
