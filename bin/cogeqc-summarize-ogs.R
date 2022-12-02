#!/usr/bin/env Rscript
library(parallel)
library(cogeqc)

# A quick function to pull out annotations and reformat for cogeqc
get_annots <- 
    function(spp = NULL, ids = NULL, anns = NULL){
        # Name annotations according to their corresponding protein accessions
        names(anns) <- ids
        
        # Remove trailing semicolons
        anns <- gsub(";$", "", anns)
        
        # Convert to a list, with a named entry of multiple annotations per accession
        anns <- strsplit(anns, split = ";")
        
        # Reformat for cogeqc (data.frame with row for genes, row for annotation)
        anns <- data.frame(t(stack(anns)[,2:1]))
        
        # Rename rows
        rownames(anns) <- c("Gene", "Annotation")
        
        return(anns)
    }

# Slighly modify the cogeqc read_orthofinder_stats() function to ignore the
# duplications file, since we're comparing orthogroups prior to inferring
# the species tree
read_orthofinder_stats <-
    function (stats_dir = NULL, species = NULL)
        {
            og_overlap <- file.path(stats_dir, "Orthogroups_SpeciesOverlaps.tsv")
            og_overlap <- read.csv(og_overlap, sep = "\t", header = TRUE, 
                row.names = 1)
            rownames(og_overlap) <- species
            colnames(og_overlap) <- species
            stats <- file.path(stats_dir, "Statistics_PerSpecies.tsv")
            stats <- read.csv(stats, sep = "\t", nrows = 10, header = TRUE, 
                row.names = 1)
            stats <- as.data.frame(t(stats))[, -c(3, 5, 6, 7)]
            colnames(stats) <- c("n_genes", "n_genes_in_ogs", "perc_genes_in_ogs", 
                "n_ss_ogs", "n_genes_in_ss_ogs", "perc_genes_in_ss_ogs")
            rownames(stats) <- species
            stats <- cbind(data.frame(Species = rownames(stats)), stats)
            stats$Species <- as.factor(stats$Species)
            rownames(stats) <- NULL
            result_list <- list(stats = stats, og_overlap = og_overlap)
            return(result_list)
        }

# Additionally, edit the assess_orthogroups function to run in parallel,
# using mclappy from the parallel package
assess_orthogroups <- 
function (orthogroups = NULL, annotation = NULL, correct_overclustering = TRUE, mc.cores = 1) 
{
    og_list <- split(orthogroups, orthogroups$Species)
    og_list <- mclapply(seq_along(og_list), mc.cores = mc.cores, function(x) {
        species <- names(og_list)[x]
        idx <- which(names(annotation) == species)
        merged <- merge(og_list[[x]], annotation[[idx]])
        names(merged)[4] <- "Annotation"
        H <- calculate_H(merged, correct_overclustering = correct_overclustering)
        names(H) <- c("Orthogroups", paste0(species, "_score"))
        return(H)
    })
    merge_func <- function(x, y) {
        merge(x, y, by = "Orthogroups", all = TRUE)
    }
    final_df <- Reduce(merge_func, og_list)
    means <- apply(final_df[, -1], 1, mean, na.rm = TRUE)
    final_df$mean_score <- means
    return(final_df)
}

args = commandArgs(trailingOnly=TRUE)

# Get the base directory to where the orthofinder results are.
# Again, specified from the commandline
og_dir <- args[1]

# Pull out the inflation parameter from the filepath
inflation <- gsub(".*_", "", og_dir)

# Now, create a variable using these to specify paths to the orthofinder 
# orthogroups file
og_file <- paste0(og_dir, '/Orthogroups/Orthogroups.tsv')

# Specify the path to orthofinders orthogroup summary stats to be used by cogeqc.
og_stat_dir <- paste0(og_dir, '/Comparative_Genomics_Statistics/')

# Go ahead and read in the orthogroups file
orthogroups <- read_orthogroups(og_file)

# Strip trailing text from species name - may not need in full implementation. 
# Names are determined in orthofinder using the file name, so just include 
# the species here. 
orthogroups$Species <- gsub('[.].*', '', orthogroups$Species)

# Pull out the list of species - we are only running this script on a 
# on a subset of species, so we want to be sure that we're not reading in 
# annotations for species other than those we're testing inflation 
# parameters with. 
species <- unique(orthogroups$Species)

# and remove the Species name from the gene name - this will create issues when
# pairing with the annotations. 
orthogroups$Gene <- gsub('^(?:[^_]*_)*\\s*(.*)', '\\1', orthogroups$Gene)

# Read in the annotations. 
annots <- list.files('./', pattern = "annotations.tsv")
spps <- gsub("-protein-annotations.tsv", "", annots)
annots <- annots[which(spps %in% unique(species))]
interpro <- list()

for(i in 1:length(annots)){
    spp <- gsub("-protein-annotations.tsv", "", annots[i]) # So we can name the entry
    
    # read in their annotations
    annotations <- read.delim(paste0('./', annots[i]), sep = "\t", header = T)
    
    # Pull out the InterPro annotations 
    interpro[[i]] <- get_annots(spp, annotations$From, annotations$InterPro)
}

# Now name all the entries according to their OrthoFinder species ID (second column)
names(interpro) <- spps[which(spps %in% species)]

# Great, now we can pair these annotations with the orthogroups, assessing how 
# well each inflation parameter infers sensible orthogroups with respect to the
# homogeneity and dispersal of annotations
mc_cores <- detectCores()-2
makeForkCluster(nnodes = mc.cores)
interpro_assess <- assess_orthogroups(orthogroups, interpro, mc.cores = mc_cores)

# Read in the orthofinder orthogroup statistics
ortho_stats <- read_orthofinder_stats(og_stat_dir, spps[which(spps %in% species)])

# Let's focus on a subset of particularly informative summary statistics 
# that can characterize the "quality" of our inferred orthogroups. 
# Namely, we'll look closely at:
#   1) The number of orthogroups
#   2) The proportion of orthogroups with >= 4 spp
#   3) The proportion of OGs with all species
#   4) The mean per-species gene count per orthogroup for OGs with >= 4 spp
#   5) The mean OG composition score for InterPro protein domains
#   6) The mean per-species percentage genes in orthogroups
#   7) The mean per-species percentage of single-species orthogroups
#   8) The mean pairwise species overlap of orthogroups

# Determine how many species are in each orthogroup
og_freqs <- table(as.factor(unique(orthogroups[,1:2])$Orthogroup))

# Get the number of gene copies per species, per orthogroup
per_spp_og_counts <- table(orthogroups[,1:2])

# And from this, get the mean per-species count per orthogroup with at least 4 spp
per_spp_og_counts <- rowMeans(per_spp_og_counts[!rowSums(per_spp_og_counts == 0) > 4,])

# pull out proportional overlap between species
overlap <- ortho_stats$og_overlap/do.call(pmax, ortho_stats$og_overlap)
overlap <- overlap[lower.tri(overlap)]

og_quality <- 
    data.frame(
        inflation_param = inflation, 
        num_ogs = length(unique(orthogroups$Orthogroup)),
        num_ogs_gt_4spp = length(which(og_freqs >= 4)),
        num_ogs_all_spp = length(which(og_freqs == max(og_freqs))),
        per_spp_4spp_og_counts = mean(per_spp_og_counts),
        interpro_score = mean(interpro_assess$mean_score),
        perc_genes_in_ogs = mean(ortho_stats$stats$perc_genes_in_ogs),
        perc_genes_in_ss_ogs = mean(ortho_stats$stats$perc_genes_in_ss_ogs),
        pairwise_overlap = mean(overlap)
    )

# Write out to a tsv.
write.table(og_quality, file = paste0('MCL-Inflation-', inflation, '-cogeqc-summary.tsv'), col.names = T, row.names = F, sep = '\t', quote = F)

sink("version.txt")
packageVersion("cogeqc")
sink()
