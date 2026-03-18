# Load the necessary libraries
require(dplyr)
require(phytools)
Rcpp::sourceCpp("calculate_dist_stats.cpp")
source("centroid_distance_functions.R")

# calc_universal_dists:
# Compute phylo-GLS transformed data, pairwise Mahalanobis distances,
# inverse covariance, centroid distances, and cophenetic distances.
# Runs on ALL gene families regardless of reference species.
calc_universal_dists <-
  function(gene_family, gf_stats_path, keep_stats, out_dir) {
    # Prep output directories for universal outputs
    dir.create(paste0(out_dir, "/protein-dist-mats/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/protein-phylo-dist-mats/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(out_dir, "/phylo-corrected-data/"),
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

    # Compute inverse covariance for centroid distance calculation
    inv_cov <- compute_inverse_covariance(transf_data) # nolint

    # Compute centroid distances for all proteins
    centroid_results <- calc_centroid_dists(transf_data, inv_cov, gf_tree)

    # Get the phylogenetic distance matrix
    prot_phylo_dists <- ape::cophenetic.phylo(gf_tree)

    return(list(
      transf_data = transf_data,
      dist_mat = dist_mat,
      prot_phylo_dists = prot_phylo_dists,
      inv_cov = inv_cov,
      centroid_results = centroid_results,
      gf_tree = gf_tree
    ))
  }

# calc_ref_dists:
# Compute reference-specific distances, permutation tests, z-scores,
# Wilcoxon tests, and final summary table. Only called when the
# reference species has proteins in the gene family.
calc_ref_dists <-
  function(universal_results, ref_spp, gene_family, out_dir) {
    # Prep output directories for reference-specific outputs
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

    dist_mat <- universal_results$dist_mat
    prot_phylo_dists <- universal_results$prot_phylo_dists

    # Pull out the distances between each species and our reference species
    focal_dists_idx <- which(grepl(ref_spp, rownames(dist_mat)))
    focal_dists <-
      matrix(dist_mat[-focal_dists_idx, focal_dists_idx],
             nrow = nrow(dist_mat[-focal_dists_idx, ]),
             ncol = length(focal_dists_idx),
             dimnames = list(rownames(dist_mat)[-focal_dists_idx],
                             colnames(dist_mat)[focal_dists_idx]))
    focal_prots <- rownames(focal_dists)

    # Conduct permutation tests to assess whether individual non-reference
    # species / reference species protein pairs are significantly less
    # dissimilar than expected.
    per_prot_dist_res <-
      dist_permute_test(dist_mat, focal_dists, n_permutations = 10000) # nolint

    # For each gene family, obtain the mean, median, and SD of the
    # distances to reference proteins
    stats_list <-
      calculate_dist_stats(focal_dists) # nolint
    focal_stats <- do.call(cbind, stats_list)
    dimnames(focal_stats) <-
      list(focal_prots, paste0(names(stats_list), "_dist_to_ref"))
    focal_dist_prot_res <-
      data.frame(protein = focal_prots, focal_stats, row.names = NULL)
    focal_dist_spp_res <-
      data.frame(species = gsub("_[^_]+$", "", focal_prots),
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

    # Identify the non-reference species for each protein pair
    nonref_spp <-
      gsub("_[^_]+$", "", per_prot_dist_res$observation)

    # Extract protein IDs
    nonref_prot <-
      gsub(".*\\_", "", per_prot_dist_res$observation) |>
      gsub(pattern = "\\|\\.p", replacement = "\\.p") |>
      gsub(pattern = "^[^|]*\\|", replacement = "") |>
      gsub(pattern = "\\|.*", replacement = "")
    ref_prot <-
      gsub(pattern = "^[^|]*\\|", replacement = "",
           per_prot_dist_res$reference) |>
      gsub(pattern = "\\|.*", replacement = "")

    # Get pairwise phylogenetic distances among proteins
    phyl_dists <-
      apply(per_prot_dist_res[, 1:2], 1, function(x) {
        obs <- x[1]
        ref <- x[2]
        prot_phylo_dists[obs, ref]
      })

    # Assemble into a combined summary table
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
        pvalue_colwise = per_prot_dist_res$pvalue_within_nonref
      )

    return(list(
      prot_dists_to_ref = focal_dist_prot_res,
      spp_dists_to_ref = obs_spp_dists,
      protein_pvals = per_prot_signif,
      species_pvals = per_spp_signif,
      per_protein_dist_res = per_prot_dist_res,
      final_summary_table = final_summary_table
    ))
  }

# genefam_aa_conservation:
# Main wrapper function for a single gene family.
# Always runs universal + centroid analysis on all families.
# Conditionally runs reference-specific analysis when ref_spp is present.
genefam_aa_conservation <-
  function(gene_family, ref_spp = ref_spp, aa_stat_basedir = aa_stat_basedir,
           keep_stats =
           c("molecular_weight", "aromaticity", "instability", "flexibility",
             "gravy_bm", "isoelectric_point", "charge_at_pH_7", "helix_fract",
             "sheet_fract", "molar_ext_coef_cysteines"),
           out_dir = "gf-aa-multivar-distances") {

    # --- Tier 1: Universal analysis (all gene families) ---
    universal_res <-
      calc_universal_dists(
        gene_family = gene_family,
        gf_stats_path = paste0(aa_stat_basedir,
                               gene_family["family"],
                               "_summary_statistics.csv"),
        keep_stats = keep_stats,
        out_dir = out_dir
      )

    # Write universal outputs
    write.table(universal_res$transf_data, sep = "\t",
                file = paste0(out_dir, "/phylo-corrected-data/",
                              gene_family["family"],
                              "_phylo_corr_dat.tsv"),
                quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(universal_res$dist_mat, sep = "\t",
                file = paste0(out_dir, "/protein-dist-mats/",
                              gene_family["family"],
                              "_protein_dists.tsv"),
                quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(universal_res$prot_phylo_dists, sep = "\t",
                file = paste0(out_dir, "/protein-phylo-dist-mats/",
                              gene_family["family"],
                              "_phylo_dists.tsv"),
                quote = FALSE, row.names = TRUE, col.names = TRUE)

    # Write centroid outputs (always)
    write_centroid_outputs(universal_res$centroid_results,
                          gene_family, out_dir)

    # --- Tier 2: Reference-specific analysis (conditional) ---
    # Check that there are enough non-reference species with ≥2 proteins
    # for the per-species Wilcoxon test (needs ≥2 non-ref species to compare)
    all_prots <- rownames(universal_res$dist_mat)
    nonref_prots <- all_prots[!grepl(ref_spp, all_prots)]
    nonref_spp_counts <- table(gsub("_[^_]+$", "", nonref_prots))
    n_nonref_with_enough <- sum(nonref_spp_counts >= 2)

    if (ref_spp != "none" &&
        any(grepl(ref_spp, rownames(universal_res$dist_mat))) &&
        n_nonref_with_enough >= 2) {
      ref_res <-
        calc_ref_dists(
          universal_results = universal_res,
          ref_spp = ref_spp,
          gene_family = gene_family,
          out_dir = out_dir
        )

      write.table(ref_res$prot_dists_to_ref, sep = "\t",
                  file = paste0(out_dir, "/protein-dists-to-reference/",
                                gene_family["family"],
                                "_protein_dists.tsv"),
                  quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(ref_res$spp_dists_to_ref, sep = "\t",
                  file = paste0(out_dir, "/species-dists-to-reference/",
                                gene_family["family"],
                                "_species_dists.tsv"),
                  quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(ref_res$protein_pvals, sep = "\t",
                  file = paste0(out_dir, "/protein-pvals/",
                                gene_family["family"],
                                "_protein_reference_dist_pvals.tsv"),
                  quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(ref_res$species_pvals, sep = "\t",
                  file = paste0(out_dir, "/species-pvals/",
                                gene_family["family"],
                                "_species_reference_dist_pvals.tsv"),
                  quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(ref_res$per_protein_dist_res, sep = "\t",
                  file = paste0(out_dir, "/pairwise-protein-dist-perm-test/",
                                gene_family["family"],
                                "_protein_protein_dist_permutation_test.tsv"),
                  quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(ref_res$final_summary_table, sep = "\t",
                  file = paste0(out_dir, "/final_protein_pair_summary_tables/",
                                gene_family["family"],
                                "_final_summary_table.tsv"),
                  quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
  }
