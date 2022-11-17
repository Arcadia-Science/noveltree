#!/usr/bin/env Rscript
# Libraries to get annotations using the uniprot web service
# The docker image on quay.io is out of date, so pull a docker image of
# bioconductor, and then install from the binary that comes with it.

# Be sure to set the library path to a temporary directory first. 
dir.create('./Rlibs', showWarnings = FALSE)
.libPaths('./Rlibs')
BiocManager::install('cogeqc', lib = './Rlibs')
library(cogeqc, lib = './Rlibs')

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
            colnames(stats) <- c("N_genes", "N_genes_in_OGs", "Perc_genes_in_OGs", 
                "N_ssOGs", "N_genes_in_ssOGs", "Perc_genes_in_ssOGs")
            rownames(stats) <- species
            stats <- cbind(data.frame(Species = rownames(stats)), stats)
            stats$Species <- as.factor(stats$Species)
            rownames(stats) <- NULL
            result_list <- list(stats = stats, og_overlap = og_overlap)
            return(result_list)
        }

args = commandArgs(trailingOnly=TRUE)

# Identify the MCL clustering parameter that we're summarizing. 
inflation <- args[1]


# Get the base directory to where the orthofinder results are.
# Again, specified from the commandline
ofDir <- args[2]

# Specify the path to protein annotations to be used by cogeqc.
# Again, specified from the commandline
annotDir <- args[3]
# annotDir <- './results/protein-annotations/'

# Now, create a variable using these to specify paths to the orthofinder 
# orthogroups file
ogFile <- paste0(ofDir, 'Results_Inflation_', inflation, '/Orthogroups/Orthogroups.tsv')

# Specify the path to orthofinders orthogroup summary stats to be used by cogeqc.
ogStatDir <- paste0(ofDir, 'Results_Inflation_', inflation, '/Comparative_Genomics_Statistics/')

# Go ahead and read in the orthogroups file
orthogroups <- read_orthogroups(ogFile)

# Strip trailing text from species name - may not need in full implementation. 
# Names are determined in orthofinder using the file name, so just include 
# the species here. 
orthogroups$Species <- gsub('[.].*', '', orthogroups$Species)

# and remove the Species name from the gene name - this will create issues when
# pairing with the annotations. 
orthogroups$Gene <- gsub('^(?:[^_]*_)*\\s*(.*)', '\\1', orthogroups$Gene)

# Read in the annotations. 
annots <- list.files(annotDir, pattern = "annotations.tsv")
interpro <- list()

for(i in 1:length(annots)){
    spp <- gsub("-protein-annotations.tsv", "", annots[i]) # So we can name the entry
    annotations <- read.delim(paste0(annotDir, annots[i]), sep = "\t", header = T)
    
    # Pull out the specific annotations individually. 
    interpro[[i]] <- get_annots(spp, annotations$From, annotations$InterPro)
}

# Now name all the entries according to their OrthoFinder species ID (second column)
species <- gsub("-protein-annotations.tsv", "", annots)
names(interpro) <- species

# Great, now we can pair these annotations with the orthogroups, assessing how 
# well each inflation parameter infers sensible orthogroups with respect to the
# homogeneity and dispersal of annotations
interpro_assess <- assess_orthogroups(orthogroups, interpro)

# Read in the orthofinder orthogroup statistics
ortho_stats <- read_orthofinder_stats(ogStatDir, species)

# Let's focus on a subset of particularly informative summary statistics 
# that can characterize the "quality" of our inferred orthogroups. 
# Namely, we'll look closely at:
#   1) The number of orthogroups
#   2) The proportion of orthogroups with >= 4 spp
#   3) The proportion of OGs with 90% of species
#   4) The proportion of OGs with all species
#   5) The mean OG composition score for InterPro protein domains
#   6) The mean per-species percentage genes in orthogroups
#   7) The mean per-species percentage of single-species orthogroups
#   8) The mean pairwise species overlap of orthogroups

og.freqs <- summary(as.factor(unique(orthogroups[,1:2])$Orthogroup))

# pull out proportional overlap between species
overlap <- ortho_stats$og_overlap/do.call(pmax, ortho_stats$og_overlap)
overlap <- overlap[lower.tri(overlap)]

og_quality <- 
    data.frame(
        NumOGs = length(unique(orthogroups$Orthogroup)),
        NumOGs_GT_4spp = length(which(og.freqs >= 4)),
        NumOGs_90perc_spp = length(which(og.freqs >= quantile(og.freqs, 0.9)[[1]])),
        NumOGs_All_spp = length(which(og.freqs == max(og.freqs))),
        interpro_score = mean(interpro_assess$Mean_score),
        perc_genes_in_OGs = mean(ortho_stats$stats$Perc_genes_in_OGs),
        perc_genes_in_ssOGs = mean(ortho_stats$stats$Perc_genes_in_ssOGs),
        pairwise_overlap = mean(overlap)
    )
# Write out to a tsv.
write.table(og_quality, file = paste0(annotDir, spp, '-protein-annotations.tsv'), col.names = T, row.names = F, sep = '\t', quote = F)

sink("version.txt")
packageVersion("cogeqc")
sink()
