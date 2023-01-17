#!/usr/bin/env Rscript
# Library to get annotations using the uniprot web service
library(UniProt.ws)

args = commandArgs(trailingOnly=TRUE)

# Pull the string specifying annotations to download
annots_to_download <- args[1]

# Pull out the species names - this is to be specified from the commandline
spp <- args[2]

# Get the filepath to the protein accessions we'll be mapping.
# Also specified from the commandline.
ids <- args[3]

# Because this is going to be a fair bit of information,
# let's make a new directory to house these outputs/gene ontologies,
# one for each species.
dir.create(file.path(spp), showWarnings = FALSE)

# Specify the fields that we would like to download.
common_cols <- c('organism_name', 'organism_id')
seq_cogeqc <-
    c('xref_interpro', 'xref_supfam', 'xref_prosite', 'xref_hogenom',
    'xref_oma', 'xref_orthodb')

seq_meta <-
    c('protein_name', 'length', 'mass', 'cc_mass_spectrometry',
    'virus_hosts', 'organelle', 'cc_rna_editing')

seq_gos <-
    c('go_p', 'go_c', 'go_f', 'go_id')

seq_funcs <-
    c('absorption', 'ft_act_site', 'cc_activity_regulation',
    'ft_binding', 'cc_catalytic_activity', 'cc_cofactor',
    'ft_dna_bind', 'ec', 'cc_function', 'kinetics', 'cc_pathway',
    'ph_dependence', 'redox_potential', 'rhea', 'ft_site',
    'temp_dependence')

seq_inter <-
    c('cc_interaction', 'cc_subunit')

seq_interdbs <-
    c('xref_biogrid', 'xref_corum', 'xref_complexportal', 'xref_dip',
    'xref_elm', 'xref_intact', 'xref_mint', 'xref_string')

seq_biotech <-
    c('cc_allergen','cc_biotechnology', 'cc_disruption_phenotype',
    'cc_disease','ft_mutagen', 'cc_pharmaceutical',
    'cc_toxic_dose')

seq_localize <-
    c('ft_intramem', 'cc_subcellular_location', 'ft_topo_dom', 'ft_transmem')

seq_ptm <-
    c('ft_chain', 'ft_crosslnk', 'ft_disulfid', 'ft_carbohyd',
    'ft_init_met', 'ft_lipid', 'ft_mod_res', 'ft_peptide',
    'cc_ptm', 'ft_propep', 'ft_signal', 'ft_transit')

seq_xptmdb <-
    c('xref_carbonyldb', 'xref_depod', 'xref_glyconnect', 'xref_glygen',
    'xref_metosite', 'xref_phosphositeplus', 'xref_swisspalm', 'xref_unicarbkb',
    'xref_iptmnet')

seq_domains <-
    c('ft_coiled', 'ft_compbias', 'cc_domain', 'ft_domain',
    'ft_motif', 'protein_families', 'ft_region', 'ft_repeat',
    'ft_zn_fing')

seq_xfamdom <-
    c('xref_cdd', 'xref_disprot', 'xref_gene3d', 'xref_hamap',
    'xref_ideal', 'xref_interpro', 'xref_panther', 'xref_pirsf',
    'xref_prints', 'xref_prosite', 'xref_pfam', 'xref_prodom',
    'xref_sfld', 'xref_smart', 'xref_supfam', 'xref_tigrfams')

seq_xseqdbs <-
    c('xref_ccds', 'xref_embl', 'xref_pir', 'xref_refseq')

seq_x3ddbs <-
    c('xref_alphafolddb', 'xref_bmrb', 'xref_pcddb', 'xref_pdb',
    'xref_pdbsum', 'xref_sasbdb', 'xref_smr')

seq_xenzpath <-
    c('xref_brenda', 'xref_biocyc', 'xref_pathwaycommons', 'xref_plantreactome',
    'xref_reactome', 'xref_sabio-rk', 'xref_signor', 'xref_signalink',
    'xref_unipathway')

seq_xortho <-
    c('xref_genetree', 'xref_hogenom', 'xref_inparanoid', 'xref_ko',
    'xref_oma', 'xref_orthodb', 'xref_phylomedb', 'xref_treefam',
    'xref_eggnog')

annotations <- list(
    seq_cogeqc, seq_meta, seq_gos, seq_funcs, seq_inter,
    seq_interdbs, seq_biotech, seq_localize, seq_ptm,
    seq_xptmdb, seq_domains, seq_xfamdom, seq_xseqdbs,
    seq_x3ddbs, seq_xenzpath, seq_xortho
    )

annot_names <-
    c('-cogeqc-annotations.tsv', '-prot-metadat.tsv', '-prot-gene-ontologies.tsv', 
    '-prot-functions.tsv', '-prot-interactions.tsv', '-prot-interactions-xref.tsv',
    '-prot-biotech-annots.tsv', '-prot-localization.tsv', '-prot-post-trans-mods.tsv', 
    '-prot-post-trans-mods-xref.tsv', '-prot-fams-domains.tsv', '-prot-fams-domains-xref.tsv',
    '-prot-seq-dbs-xref.tsv', '-prot-3d-dbs-xref.tsv', '-prot-enzyme-paths-xref.tsv', 
    '-prot-orthology-dbs.tsv')

# get the list of accessions for this species.
accessions <- read.table(ids, sep = '\t')$V1

# Get the species tax ID by pulling out into for the first accession
tax_id <- UniProt.ws::queryUniProt(paste0("accession:", accessions[1]),
                                   fields = 'organism_id')[[1]]
# Pull down annotations for proteins in this species
up <- UniProt.ws::UniProt.ws(tax_id)

# And pull down annotations, removing rows for proteins without any annotations.
# We now use the annots_to_download user variable to determine which we are
# downloading. 
if (annots_to_download == "all") {
    anns <- 1:16
} else if (annots_to_download == "none") {
    anns <- 1
} else {
    anns <- annots_to_download
}
for(i in anns){
    annots <- UniProt.ws::select(up, accessions, c(common_cols, annotations[[i]]), 'UniProtKB')
    to_drop <- which(rowSums(is.na(annots[,-c(1:4)])) == ncol(annots[,-c(1:4)]))

    if(length(to_drop) < nrow(annots)){
        if(length(to_drop) > 0){
            annots <- annots[-to_drop,]
        }

        colnames(annots) <-
        gsub("[.][.]", "_", colnames(annots)) |>
            gsub(pattern = "[.]", replacement = "_") |>
            sub(pattern = "_$", replacement = "")

        # Write out to a tsv assuming we haven't removed every protein
        # Put the cogeqc annotations in its own directory
        if(i == 1){
            write.table(annots, file = paste0(spp, annot_names[i]),
                        col.names = T, row.names = F, sep = '\t', quote = F)
        }else{
            write.table(annots, file = paste0(spp, '/', spp, annot_names[i]),
                        col.names = T, row.names = F, sep = '\t', quote = F)
        }
    }
}

sink("version.txt")
packageVersion("UniProt.ws")
sink()
