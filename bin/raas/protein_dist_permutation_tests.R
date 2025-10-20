# Load the necessary libraries
require(dplyr)
require(tidyr)
require(Rcpp)
require(RcppArmadillo)

Rcpp::sourceCpp("./perm_test_across_non_ref.cpp")
Rcpp::sourceCpp("./perm_test_within_non_ref.cpp")
# `dist_permute_test`:
# Function to conduct permutation tests of distances within gene family.
# Input:
# `dist_mat`: the complete distance matrix, used in permutation test to
#             establish expected distributions.
# `focal_dists`: subsetted distance matrix between focal proteins, with
#                reference proteins in columns, and all other proteins
#                in rows.
# `n_permutations`: Default 10000, number of permutations for the test.
#
# Process: Performs permutation tests on each column of `focal_dists` using
#          `perm_test_across_non_ref` and `perm_test_within_non_ref`,
#          calculates ranks, p-values for each distance.
#
# Output: Returns a reshaped and sorted DataFrame with focal proteins,
#         reference proteins, distances, ranks, p-values.
dist_permute_test <- function(dist_mat, focal_dists, n_permutations = 10000) {
  # Extract protein names from row names
  obs_names <- rownames(focal_dists)

  # Perform permutation test. Here, p-values are calculated per-reference
  # protein, and correspond to the probability of observing a distance as
  # small or smaller than observed for a given non-reference protein. In
  # other words, permutations are conducted within column (per-reference
  # protein), across non-reference proteins.

  # Create a list of data frames, each containing the ranks and p-values for
  # each reference protein
  ref_dists_p_values_list <-
    lapply(seq_len(ncol(focal_dists)),
           function(reference) {
             p_values <-
                perm_test_across_non_ref(focal_dists[,reference], n_permutations) # nolint
             data.frame(
               observation = obs_names,
               reference = colnames(focal_dists)[reference],
               distance = focal_dists[, reference],
               rank_distance = rank(focal_dists[, reference]),
               pvalue_across_nonref = p_values,
               row.names = NULL
             )
           })

  # Combine the individual data frames into one
  p_values_df <- do.call(rbind, ref_dists_p_values_list)

  # Now, for each species/protein test whether they are more significantly more
  # similar to each human protein than they are to non-human proteins
  # In other words, permutations are conducted within column (per-reference
  # protein), across non-reference proteins
  obs_dists_p_values <-
    perm_test_within_non_ref(dist_mat, focal_dists, n_permutations) # nolint
  obs_dists_p_values <- setNames(
    lapply(seq_len(ncol(obs_dists_p_values)),
           function(x) obs_dists_p_values[, x]), colnames(obs_dists_p_values)
  )
  obs_dists_p_values <-
    do.call(rbind,
            lapply(names(obs_dists_p_values), function(name) {
              df <- as.data.frame(stack(obs_dists_p_values[[name]]))
              df$reference <- name
              df
            }))

  # Prepare to combine:
  obs_dists_p_values$key <-
    paste0(obs_dists_p_values$ind, obs_dists_p_values$reference, sep = "_")
  key <- paste0(p_values_df$observation, p_values_df$reference, sep = "_")
  # And order to match the rest of the table:
  obs_dists_p_values <-
    obs_dists_p_values[match(key, obs_dists_p_values$key), ]

  # And combine:
  p_values_df$pvalue_within_nonref <- obs_dists_p_values$values

  return(p_values_df)
}
