# Load the necessary libraries
require(dplyr)
require(geiger)
require(phytools)
Rcpp::sourceCpp("calculate_dist_stats.cpp")
# The following function calculates pairwise multivariate distances between
# species, pulls out the distances between each species and the reference
# species (i.e. human) and conducts statistical tests of the hypothesis
# that each protein or each species is exceptionally similar
calc_prot_spp_dists <-
  function(gene_family, gf_stats_path,
           ref_spp, spp_tree,
           max_treepl_treesize,
           keep_stats, out_dir,
           clinvar) {
    # gf_stats_path: the file path to the file containing protein summary
    # statistics for a single gene family
    # gf_tree_path: the file path to the corresponding gene family tree
    # ref_spp: the species name of the "reference" species, the protein of
    # which we will compare all non-reference proteins to
    # spp_tree: the time-calibrated species tree, that if provided, will be
    # used to time-calibrate the gene family tree
    # max_treepl_treesize: the maximum number of proteins in a gene family
    # tree to use treePL to time calibrate. Larger trees will be
    # time-calibrated with pathD8
    # - treePL: more accurate, but slower for larger trees.
    #   - https://doi.org/10.1093/bioinformatics/bts492
    # - pathD8: less accurate, but significantly faster for large trees.
    #   - https://doi.org/10.1080/10635150701613783
    # keep_stats: the AA summary statistics we wish to retain in analyses
    # out_dir: the base directory to write the output files to
    # Prep output directories
    dir.create(paste0(out_dir, "/congruified-gfts/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/protein-dist-mats/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/protein-phylo-dist-mats/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/phylo-corrected-data/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/protein-dists-to-reference/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/species-dists-to-reference/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/protein-pvals/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/species-pvals/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/pairwise-protein-dist-perm-test/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/final_protein_pair_summary_tables/"),
               recursive = TRUE, showWarnings = FALSE)

    # Read in the gene family protein stats and tree
    gf_stats <- as.matrix(read.csv(gf_stats_path, row.names = 1))
    gf_tree <- ape::read.tree(gene_family["gft"])

    # Make sure that the sequences in each match
    gf_stats <-
      gf_stats[which(rownames(gf_stats) %in% gf_tree$tip.label), keep_stats]
    gf_stats <-
      gf_stats[match(gf_tree$tip.label, rownames(gf_stats)), ]

    # Scale the data
    gf_stats <- apply(gf_stats, 2, scale)
    rownames(gf_stats) <- gf_tree$tip.label

    # Check if any columns have been replaced with NaN's due to all
    # values in the original matrix being zero, and replace these with 0.
    gf_stats[, which(colSums(!is.finite(gf_stats)) == nrow(gf_stats))] <- 0
    gf_stats <- na.omit(gf_stats)
    gf_tree <- ape::keep.tip(gf_tree, rownames(gf_stats))

    if (!is.null(spp_tree) && ape::is.ultrametric(spp_tree)) {
      pathd8_path <-
        paste0(Sys.getenv("CONDA_PREFIX"), "/bin:", Sys.getenv("CONDA_PREFIX"),
               "/lib64:", Sys.getenv("LD_LIBRARY_PATH"))
      Sys.setenv(LD_LIBRARY_PATH = pathd8_path)
      dir.create(paste0(gene_family["family"]))
      setwd(paste0(gene_family["family"]))
      # Determine if we're using treePL or PATHd8
      if (length(gf_tree$tip.label) < max_treepl_treesize) {
        scale_method <- "treePL"
      } else {
        scale_method <- "PATHd8"
      }
      taxonomy <-
        matrix(gsub("_.*", "", gf_tree$tip.label),
               dimnames = list(gf_tree$tip.label, NULL),
               ncol = 1)
      gf_tree <-
        suppressWarnings(geiger::congruify.phylo(reference = spp_tree,
                                                 target = gf_tree,
                                                 scale = scale_method,
                                                 taxonomy = taxonomy)$phy)
      setwd("../")
      unlink(paste0(gene_family["family"]), recursive = TRUE)
    }
    # And add a small value to the tree edge lengths to ensure it plays
    # nicely in the case that there are any zero-branch lengths. This
    # number is based on the smallest branch lengths typically inferred
    # by phylogenetic software.
    if (min(gf_tree$edge.length) == 0) {
      gf_tree$edge.length <- gf_tree$edge.length + 1e-6
    }

    # Now calculate the pairwise mahalanobis distances between proteins,
    # using a phylogenetic GLS transformation of the protein features,
    # effectively giving us the residual trait variation after removing
    # the variance explained by phylogeny (evolutionary history) alone.
    transf_data <-
      phylo_gls_transform(gf_stats, gf_tree) # nolint

    # Note that when calculating the mahalanobis distances, the traits
    # should be stored in columns
    dist_mat <-
      pairwise_mahalanobis(transf_data) # nolint

    # Now, pull out the distances between each species and our reference species
    focal_dists_idx <- which(grepl(ref_spp, rownames(dist_mat)))
    focal_dists <-
      matrix(dist_mat[-focal_dists_idx, focal_dists_idx],
             nrow = nrow(dist_mat[-focal_dists_idx, ]),
             ncol = length(focal_dists_idx),
             dimnames = list(rownames(dist_mat)[-focal_dists_idx],
                             colnames(dist_mat)[focal_dists_idx]))
    focal_prots <- rownames(focal_dists)

    # Now, on a protein-by-protein basis, conduct permutation tests to assess
    # whether individual non-reference species / reference species protein pairs
    # are significantly less dissimilar than expected.
    per_prot_dist_res <-
      dist_permute_test(dist_mat, focal_dists, n_permutations = 10000) # nolint

    # Now, for each gene family, obtain the mean, median, and SD of the
    # distances to reference proteins
    stats_list <-
      calculate_dist_stats(focal_dists) # nolint
    focal_stats <- do.call(cbind, stats_list)
    dimnames(focal_stats) <-
      list(focal_prots, paste0(names(stats_list), "_dist_to_ref"))
    focal_dist_prot_res <-
      data.frame(protein = focal_prots, focal_stats, row.names = NULL)
    focal_dist_spp_res <-
      data.frame(species = gsub("_.*", "", focal_prots),
                 focal_stats, row.names = NULL)

    # Summarize the data getting the mean of each stat for each species
    obs_spp_dists <-
      data.frame(cbind(n_proteins =
                         summary(as.factor(focal_dist_spp_res$species)),
                       aggregate(focal_dist_spp_res[, -1],
                                 by = list(species =
                                             focal_dist_spp_res$species),
                                 mean))[c(2, 1, 3:6)],
                 row.names = NULL)

    # Test how significantly similar proteins are to reference versions:
    per_prot_signif <-
      get_prot_zscore_pvals(focal_dist_prot_res) # nolint

    # Test the extent to which species are unusually similar to reference:
    per_spp_signif <-
      per_spp_wilcox(focal_dist_spp_res) # nolint

    # Get the phylogenetic distance matrix
    prot_phylo_dists <- ape::cophenetic.phylo(gf_tree)

    # Now pull together into a "final" summary table containing key bits
    # of information and are useful for handing off to translation.

    # Identify the non-reference species for each protein pair
    # Species names always precede the first `_`:
    # remove everything that follows
    nonref_spp <-
      gsub("_.*", "", per_prot_dist_res$observation)

    # And get the UniProt protein IDs for each protein
    # Additional substitutions needed to handle exceptions for
    # non-human proteins (i.e. not from uniprot).
    # Some example proteins are shown below:

    # Examp. 1: "Nematostella-vectensis_tr|A7RG13|A7RG13-NEMVE"
    # - This is from UniProt - we want the middle identifier between pipes
    # - Human proteins will always match this format
    # - For these we simply need to remove everything before/including the
    #   first pipe, and then everything after/including the last
    # Examp 2: "Mnemiopsis-leidyi_gb|GFAT01090803.1|.p1"
    # - This is an oddball transcriptome - we want "GFAT01090803.1.p1" and so
    #   will start out by replacing instances of "|.p" with ".p" and then
    #   proceed as we did before.
    # - Some non-human proteins don't include pipes, we we start by removing
    #   the species names (everything before "_")

    # First the non-human proteins with the extra step
    nonref_prot <-
      gsub(".*\\_", "", per_prot_dist_res$observation) |>
      gsub(pattern = "\\|\\.p", replacement = "\\.p") |>
      gsub(pattern = "^[^|]*\\|", replacement = "") |>
      gsub(pattern = "\\|.*", replacement = "")
    # Then for humans
    ref_prot <-
      gsub(pattern = "^[^|]*\\|", replacement = "",
           per_prot_dist_res$reference) |>
      gsub(pattern = "\\|.*", replacement = "")

    # Get a vector containing the pairwise phylogenetic distances among
    # proteins
    phyl_dists <-
      apply(per_prot_dist_res[, 1:2], 1, function(x) {
        obs <- x[1]
        ref <- x[2]
        prot_phylo_dists[obs, ref]
      })

    # Assemble into a combined table, populating disease info with
    # empty entries to start
    final_summary_table <-
      data.frame(
        gene_family = gene_family["family"],
        nonref_species = nonref_spp,
        nonref_protein = nonref_prot,
        ref_protein = ref_prot,
        phylo_dist = phyl_dists,
        trait_dist = per_prot_dist_res$distance,
        rank_trait_dist = per_prot_dist_res$rank_distance,
        pvalue_rowwise = per_prot_dist_res$pvalue_across_nonref,
        pvalue_colwise = per_prot_dist_res$pvalue_within_nonref,
        associated_gene = NA,
        disease_mim = NA,
        disease_name = NA,
        concept_id = NA,
        source_name = NA,
        source_id = NA
      )

    # Identify whether there are any disease-associated human genes within
    # this gene family, and if so, retain this info.
    target_prots <- unique(ref_prot)

    # Only process clinvar data if it's provided (not NULL)
    if (!is.null(clinvar)) {
      disease_info <-
        clinvar[which(clinvar$uniprot_id %in% target_prots), ]
      # Make sure there are no lingering semicolons, since we'll be appending
      # multiple disease entries into a single row.
      disease_info$disease_name <-
        gsub(";", ":", disease_info$disease_name)
    } else {
      # If clinvar is NULL, create an empty data frame
      disease_info <- data.frame()
    }

    # If there are any disease-associated genes in this gene family,
    # incorporate this info into the table
    if (nrow(disease_info) > 0) {
      # Store the names of clinvar columns we'll be accessing
      # and in the appropriate order
      clinvar_columns <-
        c(
          "disease_mim",
          "disease_name",
          "concept_id",
          "source_name",
          "source_id"
        )

      for (d in unique(disease_info$uniprot_id)) {
        # Pull out the focal disease gene
        idx <-
          which(disease_info$uniprot_id == d)
        diseases <- disease_info[idx, ]
        associated_gene <- unique(disease_info$associated_genes[idx])
        # Collapse column entries for this gene, appending with a semicolon
        diseases <-
          sapply(diseases, function(x) {
            paste(x, collapse = ";")
          })
        # Add to the combined table
        idx <-
          which(final_summary_table$ref_protein == d)
        final_summary_table[idx, "associated_gene"] <- associated_gene
        final_summary_table[idx, 11:15] <-
          matrix(rep(diseases[clinvar_columns], length(idx)),
                 ncol = 5,
                 byrow = TRUE)
      }
    }
    # And return all outputs
    return(
      list(
        phylo_corrected_data = transf_data,
        protein_dist_mat = dist_mat,
        congruified_phylo = gf_tree,
        prot_phylo_dists = prot_phylo_dists,
        prot_dists_to_ref = focal_dist_prot_res,
        spp_dists_to_ref = obs_spp_dists,
        protein_pvals = per_prot_signif,
        species_pvals = per_spp_signif,
        per_protein_dist_res = per_prot_dist_res,
        final_summary_table = final_summary_table
      )
    )
  }

# The following function calls the above function for a single gene
# family, allowing calculations to be done in parallel.
genefam_aa_conservation <-
  function(gene_family, ref_spp = ref_spp, aa_stat_basedir = aa_stat_basedir,
           spp_tree = spp_tree, max_treepl_treesize = 300, clinvar = clinvar,
           keep_stats =
           c("molecular_weight", "aromaticity", "instability", "flexibility",
             "gravy_bm", "isoelectric_point", "charge_at_pH_7", "helix_fract",
             "sheet_fract", "molar_ext_coef_cysteines"),
           out_dir = "gf-aa-multivar-distances") {
    # gene_family: a named vector containing "family", the gene family ID/name,
    # and "gft", the file path to the corresponding gene family tree
    # ref_spp: the species name of the "reference" species, the protein
    # of which we will compare all non-reference proteins to
    # spp_tree: the time-calibrated species tree, that if provided, will
    # be used to time-calibrate the gene family tree
    # max_treepl_treesize: the maximum number of proteins in a gene family
    # tree to use treePL to time calibrate. Larger trees will be
    # time-calibrated with PATHd8
    # keep_stats: the AA summary statistics we wish to retain in analyses
    # out_dir: the base directory to write the output files to
    gf_dist_res <-
      calc_prot_spp_dists(
        gene_family = gene_family,
        gf_stats_path = paste0(aa_stat_basedir,
                               gene_family["family"],
                               "_summary_statistics.csv"),
        spp_tree = spp_tree,
        ref_spp = ref_spp,
        max_treepl_treesize = max_treepl_treesize,
        keep_stats = keep_stats,
        out_dir = out_dir,
        clinvar = clinvar
      )
    ape::write.tree(gf_dist_res$congruified_phylo,
                    file = paste0(out_dir, "/congruified-gfts/",
                                  gene_family["family"],
                                  "_congruified.newick"))
    write.table(gf_dist_res$phylo_corrected_data, sep = "\t",
                file = paste0(out_dir, "/phylo-corrected-data/",
                              gene_family["family"],
                              "_phylo_corr_dat.tsv"),
                quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(gf_dist_res$protein_dist_mat, sep = "\t",
                file = paste0(out_dir, "/protein-dist-mats/",
                              gene_family["family"],
                              "_protein_dists.tsv"),
                quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(gf_dist_res$prot_phylo_dists, sep = "\t",
                file = paste0(out_dir, "/protein-phylo-dist-mats/",
                              gene_family["family"],
                              "_phylo_dists.tsv"),
                quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(gf_dist_res$prot_dists_to_ref, sep = "\t",
                file = paste0(out_dir, "/protein-dists-to-reference/",
                              gene_family["family"],
                              "_protein_dists.tsv"),
                quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(gf_dist_res$spp_dists_to_ref, sep = "\t",
                file = paste0(out_dir, "/species-dists-to-reference/",
                              gene_family["family"],
                              "_species_dists.tsv"),
                quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(gf_dist_res$protein_pvals, sep = "\t",
                file = paste0(out_dir, "/protein-pvals/",
                              gene_family["family"],
                              "_protein_reference_dist_pvals.tsv"),
                quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(gf_dist_res$species_pvals, sep = "\t",
                file = paste0(out_dir, "/species-pvals/",
                              gene_family["family"],
                              "_species_reference_dist_pvals.tsv"),
                quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(gf_dist_res$per_protein_dist_res, sep = "\t",
                file = paste0(out_dir, "/pairwise-protein-dist-perm-test/",
                              gene_family["family"],
                              "_protein_protein_dist_permutation_test.tsv"),
                quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(gf_dist_res$final_summary_table, sep = "\t",
                file = paste0(out_dir, "/final_protein_pair_summary_tables/",
                              gene_family["family"],
                              "_final_summary_table.tsv"),
                quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
