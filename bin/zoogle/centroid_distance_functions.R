# Centroid-based distance analysis functions
# These functions compute family-wide centroid distances for all proteins
# in a gene family, independent of any reference species.

require(Rcpp)
Rcpp::sourceCpp("centroid_species_perm_test.cpp")

# calc_centroid_dists:
# Compute Mahalanobis distance from each protein to the family centroid,
# protein-level rank p-values, and species-level permutation p-values.
#
# Input:
#   transf_data: phylo-GLS transformed data matrix (proteins x traits)
#   inv_cov: inverse covariance matrix from compute_inverse_covariance()
#   gf_tree: gene family tree (ape phylo object, tip-pruned to match data)
#
# Output: list with:
#   centroid_dists_df: data.frame (protein, species, centroid_dist,
#                      rank_centroid_dist, pvalue_protein)
#   centroid_species_pvals: named vector of species-level p-values
calc_centroid_dists <- function(transf_data, inv_cov, gf_tree) {
  # Compute family centroid (column means of transformed data)
  centroid <- colMeans(transf_data)

  # Compute Mahalanobis distance from each protein to centroid
  diffs <- sweep(transf_data, 2, centroid)
  centroid_dists <- sqrt(rowSums((diffs %*% inv_cov) * diffs))

  # Extract species labels from protein names
  species_labels <- gsub("_[^_]+$", "", names(centroid_dists))

  # Protein-level p-value: rank-based empirical percentile
  # High = far from centroid ("divergent"), low = close ("typical")
  rank_dists <- rank(centroid_dists)
  pvalue_protein <- rank_dists / length(centroid_dists)

  # Species-level permutation test (two-tailed)
  centroid_species_pvals <-
    centroid_species_perm_test(centroid_dists, species_labels, # nolint
                              n_permutations = 1000L)

  # Build protein-level data frame
  centroid_dists_df <- data.frame(
    protein = names(centroid_dists),
    species = species_labels,
    centroid_dist = as.numeric(centroid_dists),
    rank_centroid_dist = as.numeric(rank_dists),
    pvalue_protein = as.numeric(pvalue_protein),
    row.names = NULL
  )

  return(list(
    centroid_dists_df = centroid_dists_df,
    centroid_species_pvals = centroid_species_pvals
  ))
}

# write_centroid_outputs:
# Write centroid distance results to TSV files.
#
# Input:
#   centroid_results: output from calc_centroid_dists()
#   gene_family: named vector with "family" element (OG name)
#   out_dir: base output directory
#
# Output: writes two TSV files (centroid-dists/ and centroid-summary-tables/)
write_centroid_outputs <- function(centroid_results, gene_family, out_dir) {
  dir.create(paste0(out_dir, "/centroid-dists/"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(out_dir, "/centroid-summary-tables/"),
             recursive = TRUE, showWarnings = FALSE)

  # Write per-protein centroid distances
  write.table(centroid_results$centroid_dists_df, sep = "\t",
              file = paste0(out_dir, "/centroid-dists/",
                            gene_family["family"],
                            "_centroid_dists.tsv"),
              quote = FALSE, row.names = FALSE, col.names = TRUE)

  # Build summary table with species-level p-values merged in
  spp_pvals <- centroid_results$centroid_species_pvals
  summary_table <- data.frame(
    gene_family = gene_family["family"],
    centroid_results$centroid_dists_df,
    pvalue_species = spp_pvals[centroid_results$centroid_dists_df$species],
    row.names = NULL
  )

  write.table(summary_table, sep = "\t",
              file = paste0(out_dir, "/centroid-summary-tables/",
                            gene_family["family"],
                            "_centroid_summary_table.tsv"),
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}
