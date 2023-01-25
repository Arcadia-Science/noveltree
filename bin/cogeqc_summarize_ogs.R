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

# Initialize
interpro <- list()
supfam <- list()
prosite <- list()
hogenom <- list()
oma <- list()
orthodb <- list()

for(i in 1:length(spps)){
    spp <- spps[i] # So we can name the entry

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
    interpro[[i]] <-
        if(non_missing[1] == FALSE){
            interpro[[i]] <- get_annots(spp, annotations$From, annotations$InterPro)
        }else{
            interpro[[i]] <- NA
        }
    supfam[[i]] <-
        if(non_missing[2] == FALSE){
            supfam[[i]] <- get_annots(spp, annotations$From, annotations$SUPFAM)
        }else{
            supfam[[i]] <- NA
        }
    prosite[[i]] <-
        if(non_missing[3] == FALSE){
            prosite[[i]] <- get_annots(spp, annotations$From, annotations$PROSITE)
        }else{
            prosite[[i]] <- NA
        }
    hogenom[[i]] <-
        if(non_missing[4] == FALSE){
            hogenom[[i]] <- get_annots(spp, annotations$From, annotations$HOGENOM)
        }else{
            hogenom[[i]] <- NA
        }
    oma[[i]] <-
        if(non_missing[5] == FALSE){
            oma[[i]] <- get_annots(spp, annotations$From, annotations$OMA)
        }else{
            oma[[i]] <- NA
        }
    orthodb[[i]] <-
        if(non_missing[6] == FALSE){
            orthodb[[i]] <- get_annots(spp, annotations$From, annotations$OrthoDB)
        }else{
            orthodb[[i]] <- NA
        }
}

# Now name all the entries according to their OrthoFinder species ID (second column)
names(interpro) <- spps
names(supfam) <- spps
names(prosite) <- spps
names(hogenom) <- spps
names(oma) <- spps
names(orthodb) <- spps

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
# Count the number of species for which we have each summary statistic - we can
# only calculate these if there is >= 2 species. 
og_assess_list <- list(
    list(
        og_set = orthogroups[which(orthogroups$Species %in% names(interpro)),], 
        ann_set = interpro, spp_count = length(names(interpro))),
    list(
        og_set = orthogroups[which(orthogroups$Species %in% names(supfam)),], 
        ann_set = supfam, spp_count = length(names(supfam))),
    list(
        og_set = orthogroups[which(orthogroups$Species %in% names(prosite)),], 
        ann_set = prosite, spp_count = length(names(prosite))),
    list(
        og_set = orthogroups[which(orthogroups$Species %in% names(hogenom)),], 
        ann_set = hogenom, spp_count = length(names(hogenom))),
    list(
        og_set = orthogroups[which(orthogroups$Species %in% names(oma)),], 
        ann_set = oma, spp_count = length(names(oma))),
    list(
        og_set = orthogroups[which(orthogroups$Species %in% names(orthodb)),], 
        ann_set = orthodb, spp_count = length(names(orthodb)))
)

# A quick function to run the assessment in parallel, checking that there are
# enough species
get_assessments <- 
    function(i){
        if(og_assess_list[[i]]$spp_count > 1){
            assess <- assess_orthogroups(og_assess_list[[i]]$og_set, og_assess_list[[i]]$ann_set)
            assess <- mean(assess$Mean_score)
        }else{
            assess <- NA
        }
        return(assess)
    }

# Now, run each assessment simultaneously to save time
assessment_res <- mclapply(1:6, get_assessments, mc.cores = 6)

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
        interpro_score = assessment_res[[1]],
        supfam_score = assessment_res[[2]],
        prosite_score = assessment_res[[3]],
        hogenom_score = assessment_res[[4]],
        oma_score = assessment_res[[5]],
        orthodb_score = assessment_res[[6]],
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
