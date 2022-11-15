#!/usr/bin/env Rscript

# Libraries to get annotations using the uniprot web service
# The docker image on quay.io is out of date, so pull a docker image of
# bioconductor, and then install
# Be sure to set the library path to a temporary directory first. 
dir.create('./Rlibs', showWarnings = FALSE)
.libPaths('./Rlibs')
BiocManager::install('UniProt.ws', lib = './Rlibs')
library(UniProt.ws, lib = './Rlibs')

args = commandArgs(trailingOnly=TRUE)

# Pull out the species names - this is to be specified from the commandline
spp <- args[1]

# Get the filepath to the protein accessions we'll be mapping. 
# Also specified from the commandline. 
ids <- args[2]

# Because this is going to be a fair bit of information, let's make a new directory to house these outputs/gene ontologies. 
annotDir <- './annotations/'
dir.create(file.path(annotDir), showWarnings = FALSE) # Will only create the directory if it doesn't already exist. 

# Specify the fields that we would like to download. 
seqMeta <- 
    c('organism_name', 'organism_id', 'accession', 
    'protein_name', 'length', 'mass')

gos <- 
    c('go_p', 'go_c', 'go_f', 'go_id')
    
famDomains <- 
    c('xref_disprot','xref_ideal','xref_interpro',
    'xref_panther','xref_pirsf','xref_pfam',
    'xref_sfld','xref_supfam','xref_tigrfams')

ogDBs <- 
    c('xref_inparanoid','xref_orthodb','xref_phylomedb',
    'xref_treefam','xref_eggnog')

# get the list of accessions for this species. 
accessions <- read.table(ids, sep = '\t')$V1

# Get the species tax ID by pulling out into for the first accession
taxID <- UniProt.ws::queryUniProt(paste0("accession:", accessions[1]), fields = 'organism_id')[[1]]
# Pull down annotations for proteins in this species 
up <- UniProt.ws::UniProt.ws(taxID)

# And pull down annotations, removing rows for proteins without any annotations. 
annots <- UniProt.ws::select(up, accessions, c(seqMeta, gos, famDomains, ogDBs), 'UniProtKB')
toDrop <- which(rowSums(is.na(annots[,-c(1:2)])) == ncol(annots[,-c(1:2)]))
if(length(toDrop) > 0){
    annots <- annots[-toDrop,]
}
colnames(annots) <- 
    gsub("[.][.]", "_", colnames(annots)) |> 
        gsub(pattern = "[.]", replacement = "_") |> 
        sub(pattern = "_$", replacement = "")

# Write out to a tsv.
write.table(annots, file = paste0(annotDir, spp, '-protein-annotations.tsv'), col.names = T, row.names = F, sep = '\t', quote = F)

sink("version.txt")
packageVersion("UniProt.ws")
sink()