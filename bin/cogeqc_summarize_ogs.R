#!/usr/bin/env Rscript
library(tidyverse)
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

        # Reformat for cogeqc (data.frame with column each for genes and annotation)
        anns <- stack(anns)[,2:1]
        colnames(anns) <- c("Gene", "Annotation")
        anns$Gene <- as.character(anns$Gene)

        return(anns)
    }

# A function to read in the orthogroup statistics produced by orthofinder
# that are needed for this analysis.
get_orthofinder_stats <-
    function(og_stats_dir = NULL, species = NULL){
        # Read in stats
        spp_overlap <-
            read_tsv(paste0(og_stats_dir, "Orthogroups_SpeciesOverlaps.tsv"),
                     skip = 1, col_names=c("V1", species),
                     show_col_types = FALSE)
        og_stats <-
            read_tsv(paste0(og_stats_dir, "Statistics_PerSpecies.tsv"),
                     skip = 1, n_max = 10, col_names=c("Statistic", species),
                     show_col_types = FALSE)

        og_stats <-
            data.frame(Species = as.factor(species),
                       n_genes = unlist(og_stats[1,-1]),
                       n_genes_in_ogs = unlist(og_stats[2,-1]),
                       perc_genes_in_ogs = unlist(og_stats[4,-1]),
                       n_ss_ogs = unlist(og_stats[8,-1]),
                       n_genes_in_ss_ogs = unlist(og_stats[9,-1]),
                       perc_genes_in_ss_ogs = unlist(og_stats[10,-1]),
                       row.names = NULL)

        result_list <- list(stats = og_stats, og_overlap = spp_overlap)
        return(result_list)
    }

args = commandArgs(trailingOnly=TRUE)

# Get the base directory to where the orthofinder results are.
# Again, specified from the commandline
og_dir <- args[1]

# And the minimum number of species for orthogroup phylogenetic inference
min_spp <- args[2]

# Read in the annotations.
annots <- list.files('./', pattern = "cogeqc_annotations.tsv")
spps <- gsub("_cogeqc_annotations.tsv", "", annots)

# Pull out the inflation parameter from the filepath
inflation <- gsub(".*_", "", og_dir)

# Now, create a variable using these to specify paths to the orthofinder
# orthogroups file
og_file <- paste0(og_dir, '/Orthogroups/Orthogroups.tsv')

# Specify the path to orthofinders orthogroup summary stats to be used by cogeqc.
og_stat_dir <- paste0(og_dir, '/Comparative_Genomics_Statistics/')

# Go ahead and read in the orthogroups file
orthogroups <- read_orthogroups(og_file)

# Convert periods to hyphens in species names to match annotation file naming convention.
# OrthoFinder uses periods (e.g., "Agaricus.bisporus") but annotations use hyphens (e.g., "Agaricus-bisporus").
orthogroups$Species <- gsub('[.]', '-', orthogroups$Species)

# Get the complete list of species included here
all_species <- unique(orthogroups$Species)

# Remove any orthogroup members for which we do not have annotations
orthogroups <- orthogroups[which(orthogroups$Species %in% spps),]

# Pull out the list of species - we are only running this script on a
# on a subset of species, so we want to be sure that we're not reading in
# annotations for species other than those we're testing inflation
# parameters with.
species <- unique(orthogroups$Species)

# and remove the Species name from the gene name - this will create issues when
# pairing with the annotations.
orthogroups$Gene <- gsub('^(?:[^_]*_)*\\s*(.*)', '\\1', orthogroups$Gene)
# Also handle pipe-delimited format (e.g., "tr|K8EBU5|K8EBU5-9CHLO" -> "K8EBU5")
orthogroups$Gene <- gsub('^[^|]*\\|([^|]+)\\|.*$', '\\1', orthogroups$Gene)

# Initialize
interpro <- list()

for(i in 1:length(spps)){
    spp <- spps[i] # So we can name the entry

    # read in their annotations
    annotations <- read.delim(paste0('./', spp, '_cogeqc_annotations.tsv'), sep = "\t", header = T)

    # Identify which we have annotations for this species.
    # Check if ALL annotations are either NA, empty strings, or just whitespace
    non_missing <-
        all(is.na(annotations$xref_interpro) | annotations$xref_interpro == "" | grepl("^\\s*$", annotations$xref_interpro))

    # Pull out the InterPro annotations
    interpro[[i]] <-
        if(non_missing == FALSE){
            interpro[[i]] <- get_annots(spp, annotations$accession, annotations$xref_interpro)
        }else{
            interpro[[i]] <- NA
        }
}

# Now name all the entries according to their OrthoFinder species ID (second column)
names(interpro) <- spps

# And drop species that are missing the annotations
interpro <- Filter(function(a) any(!is.na(a)), interpro)

# And lastly intersect these with the species used for MCL-testing
interpro <- interpro[which(names(interpro) %in% species)]

# Filter out species that don't have sufficient annotation coverage
# assess_orthogroups requires at least 1 orthogroup with >=2 unique annotated genes
# to calculate Sorensen-Dice scores (exact requirement from calculate_H function)
check_min_annotations <- function(spp_name, ann_df, og_df) {
    if(is.na(ann_df)[1]) return(FALSE)

    # Get orthogroups for this species
    spp_ogs <- og_df[which(og_df$Species == spp_name), ]

    # Merge with annotations
    merged <- merge(spp_ogs, ann_df, by = "Gene")
    merged <- merged[!is.na(merged$Annotation), ]

    # Count unique annotated genes per orthogroup (not rows, since genes can have multiple annotations)
    by_og <- split(merged, merged$Orthogroup)
    unique_genes_per_og <- sapply(by_og, function(x) length(unique(x$Gene)))

    # Check if at least one orthogroup has >=2 unique annotated genes
    return(sum(unique_genes_per_og >= 2) >= 1)
}

# Apply the check to filter species
interpro_valid <- sapply(names(interpro), function(s) {
    check_min_annotations(s, interpro[[s]], orthogroups)
})

interpro <- interpro[interpro_valid]

# Great, now we can pair these annotations with the orthogroups, assessing how
# well each inflation parameter infers sensible orthogroups with respect to the
# homogeneity and dispersal of InterPro domain annotations
if(length(interpro) > 1){
    og_sub <- orthogroups[which(orthogroups$Species %in% names(interpro)),]
    interpro_score <- tryCatch({
        assess <- assess_orthogroups(og_sub, interpro)
        mean(assess$Mean_score)
    }, error = function(e) return(NA))
}else{
    interpro_score <- NA
}

# Read in the orthofinder orthogroup statistics
ortho_stats <- get_orthofinder_stats(og_stats_dir = og_stat_dir, species = all_species)

# Let's focus on a subset of particularly informative summary statistics
# Determine how many species are in each orthogroup
og_freqs <- table(as.factor(unique(orthogroups[,1:2])$Orthogroup))

# Get the number of gene copies per species, per orthogroup
per_spp_og_counts <- table(orthogroups[,1:2])

# And from this, get the mean per-species count per orthogroup with at least
# the user-specified minimum number of species
per_spp_og_counts <- rowMeans(per_spp_og_counts[!rowSums(per_spp_og_counts == 0) >= min_spp,])

# pull out proportional overlap between species
overlap <- ortho_stats$og_overlap[,-1]/do.call(pmax, ortho_stats$og_overlap[,-1])
overlap <- overlap[lower.tri(overlap)]

num_ogs <- length(unique(orthogroups$Orthogroup))
og_quality <-
    data.frame(
        inflation_param = inflation,
        num_ogs = num_ogs,
        perc_ogs_gt_min_spp = length(which(og_freqs >= min_spp)) / num_ogs,
        per_spp_og_counts = mean(per_spp_og_counts),
        interpro_score = interpro_score,
        perc_genes_in_ss_ogs = mean(ortho_stats$stats$perc_genes_in_ss_ogs),
        mean_num_ss_ogs = mean(ortho_stats$stats$n_ss_ogs),
        pairwise_overlap = mean(overlap) / 100
    )

# Write out to a tsv.
write.table(og_quality, file = paste0('MCL_Inflation_', inflation, '_cogeqc_summary.tsv'), col.names = T, row.names = F, sep = '\t', quote = F)

sink("version.txt")
packageVersion("cogeqc")
sink()
