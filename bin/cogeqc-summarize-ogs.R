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
annots <- list.files('./', pattern = "cogeqc-annotations.tsv")
spps <- gsub("-cogeqc-annotations.tsv", "", annots)

# Reduce down to the species included in the MCL test dataset
spps <- spps[which(spps %in% species)]
    
# Create a function that will extract annotations for a given species. 
get_spp_annots <- 
    function(species_index){
        spp <- spps[species_index] # So we can name the entry

        # read in their annotations
        annotations <- read.delim(paste0('./', spp, '-cogeqc-annotations.tsv'), sep = "\t", header = T)
    
        # Identify which we have annotations for this species.
        non_missing <-
            c(sum(is.na(annotations$InterPro)) == length(annotations$InterPro),
            sum(is.na(annotations$SUPFAM)) == length(annotations$SUPFAM),
            sum(is.na(annotations$PROSITE)) == length(annotations$PROSITE),
            sum(is.na(annotations$HOGENOM)) == length(annotations$HOGENOM),
            sum(is.na(annotations$OMA)) == length(annotations$OMA),
            sum(is.na(annotations$OrthoDB)) == length(annotations$OrthoDB))
    
        # Pull out the InterPro annotations
        interpro <-
            if(non_missing[1] == FALSE){
                interpro <- get_annots(spp, annotations$From, annotations$InterPro)
            }else{
                interpro <- NA
            }
        supfam <-
            if(non_missing[2] == FALSE){
                supfam <- get_annots(spp, annotations$From, annotations$SUPFAM)
            }else{
                supfam <- NA
            }
        prosite <-
            if(non_missing[3] == FALSE){
                prosite <- get_annots(spp, annotations$From, annotations$PROSITE)
            }else{
                prosite <- NA
            }
        hogenom <-
            if(non_missing[4] == FALSE){
                hogenom <- get_annots(spp, annotations$From, annotations$HOGENOM)
            }else{
                hogenom <- NA
            }
        oma <-
            if(non_missing[5] == FALSE){
                oma <- get_annots(spp, annotations$From, annotations$OMA)
            }else{
                oma <- NA
            }
        orthodb <-
            if(non_missing[6] == FALSE){
                orthodb <- get_annots(spp, annotations$From, annotations$OrthoDB)
            }else{
                orthodb <- NA
            }
        annotation_res <- 
            list(
                interpro = interpro,
                supfam = supfam,
                prosite = prosite,
                hogenom = hogenom,
                oma = oma,
                orthodb = orthodb
            )
        annotation_res
    }

# Now extract annotations for each species, doing this in parallel cause we can. 
annotation_list <- mclapply(1:length(spps), get_spp_annots)

# And pull out each individual annotation, creating a list of species
interpro <- list()
supfam <- list()
prosite <- list()
hogenom <- list()
oma <- list()
orthodb <- list()

for(i in 1:length(spps)){
    interpro[[spps[i]]] <- annotation_list[[i]]$interpro
    supfam[[spps[i]]] <- annotation_list[[i]]$supfam
    prosite[[spps[i]]] <- annotation_list[[i]]$prosite
    hogenom[[spps[i]]] <- annotation_list[[i]]$hogenom
    oma[[spps[i]]] <- annotation_list[[i]]$oma
    orthodb[[spps[i]]] <- annotation_list[[i]]$orthodb
}

# And drop species for each that are missing the annotations
interpro <- Filter(function(a) any(!is.na(a)), interpro)
supfam <- Filter(function(a) any(!is.na(a)), supfam)
prosite <- Filter(function(a) any(!is.na(a)), prosite)
hogenom <- Filter(function(a) any(!is.na(a)), hogenom)
oma <- Filter(function(a) any(!is.na(a)), oma)
orthodb <- Filter(function(a) any(!is.na(a)), orthodb)

# Great, now we can pair these annotations with the orthogroups, assessing how
# well each inflation parameter infers sensible orthogroups with respect to the
# homogeneity and dispersal of annotations

# We will run all orthogroup assessments in parallel (six simulatenously) to 
# speed things up. 
# create lists of:
# 1) orthogroup subsets
# 2) annotations

og_sets <- 
    list(
        interpro = orthogroups[which(orthogroups$Species %in% names(interpro)),],
        supfam = orthogroups[which(orthogroups$Species %in% names(supfam)),],
        prosite = orthogroups[which(orthogroups$Species %in% names(prosite)),],
        hogenom = orthogroups[which(orthogroups$Species %in% names(hogenom)),],
        oma = orthogroups[which(orthogroups$Species %in% names(oma)),],
        orthodb = orthogroups[which(orthogroups$Species %in% names(orthodb)),]
    )
    
annotations <- 
    list(
        interpro = interpro,
        supfam = supfam, 
        prosite = prosite, 
        hogenom = hogenom,
        oma = oma, 
        orthodb = orthodb
    )

# Now run cogeqc on all 6 summary stats simultaneously 
assess_orthogroups_parallel <- 
    function(i){
        assess_orthogroups(og_sets[[i]], annotations[[i]])
    }

# And get the cogeqc scores
assessments <- 
    mclapply(1:6, assess_orthogroups_parallel)
names(assessments) <- 
    c('interpro', 'supfam', 'prosite',
    'hogenom', 'oma', 'orthodb')
   
interpro_assess <-
    assess_orthogroups(orthogroups[which(orthogroups$Species %in% names(interpro)),],
                       interpro)
supfam_assess <-
    assess_orthogroups(orthogroups[which(orthogroups$Species %in% names(supfam)),],
                       supfam)
prosite_assess <-
    assess_orthogroups(orthogroups[which(orthogroups$Species %in% names(prosite)),],
                       prosite)
hogenom_assess <-
    assess_orthogroups(orthogroups[which(orthogroups$Species %in% names(hogenom)),],
                       hogenom)
oma_assess <-
    assess_orthogroups(orthogroups[which(orthogroups$Species %in% names(oma)),],
                       oma)
orthodb_assess <-
    assess_orthogroups(orthogroups[which(orthogroups$Species %in% names(orthodb)),],
                       orthodb)

# Read in the orthofinder orthogroup statistics
ortho_stats <- get_orthofinder_stats(og_stats_dir = og_stat_dir, species = spps[which(spps %in% species)])

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
per_spp_og_counts <- rowMeans(per_spp_og_counts[!rowSums(per_spp_og_counts == 0) >= 4,])

# pull out proportional overlap between species
overlap <- ortho_stats$og_overlap[,-1]/do.call(pmax, ortho_stats$og_overlap[,-1])
overlap <- overlap[lower.tri(overlap)]

og_quality <-
    data.frame(
        inflation_param = inflation,
        num_ogs = length(unique(orthogroups$Orthogroup)),
        num_ogs_gt_4spp = length(which(og_freqs >= 4)),
        num_ogs_all_spp = length(which(og_freqs == max(og_freqs))),
        per_spp_4spp_og_counts = mean(per_spp_og_counts),
        interpro_score = mean(assessments$interpro$Mean_score),
        supfam_score = mean(assessments$supfam$Mean_score),
        prosite_score = mean(assessments$prosite$Mean_score),
        hogenom_score = mean(assessments$hogenom$Mean_score),
        oma_score = mean(assessments$oma$Mean_score),
        orthodb_score = mean(assessments$orthodb$Mean_score),
        perc_genes_in_ss_ogs = mean(ortho_stats$stats$perc_genes_in_ss_ogs),
        total_num_ss_ogs = sum(ortho_stats$stats$n_ss_ogs),
        mean_num_ss_ogs = mean(ortho_stats$stats$n_ss_ogs),
        pairwise_overlap = mean(overlap)
    )

# Write out to a tsv.
write.table(og_quality, file = paste0('MCL-Inflation-', inflation, '-cogeqc-summary.tsv'), col.names = T, row.names = F, sep = '\t', quote = F)

sink("version.txt")
packageVersion("cogeqc")
sink()
