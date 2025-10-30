require(Rcpp)
require(RcppArmadillo)

# A function to calculate Z-scores for normalized distances between proteins
# within gene family. The expectation is that these are a vector of distances
# of a set of proteins to some reference (e.g. human protein). However, this
# is a general function to calculate Z-scores that can be used in many other
# contexts.
calc_zscores <- function(values) {
  # values: a vector of measurements for each observation/sample
  overall_mean <- mean(values)
  overall_sd <- sd(values)

  # Calculate Z-scores
  z_scores <- (values - overall_mean) / overall_sd
  return(z_scores)
}

# A function to obtain per-protein p-values corresponding to z-scores obtained
# from normalized distances between proteins
get_prot_zscore_pvals <- function(focal_dist_prot_res) {
  # focal_dist_prot_res: a data.frame that reports, for each protein, the
  # minimum, mean, median and StdDev of distances to the reference species'
  # (e.g. human) proteins.

  # Log transform distances, as the  distributions of distances are heavily
  # right-skewed. We would like to normalize so that the distribution is
  # statistically consistent with the definition of z-scores. Add 1e-6 (a
  # sufficiently, but not excessively small number) to the minimums to avoid
  # the presence of zeroes.
  prot_mins <- log10(focal_dist_prot_res$min_dist_to_ref + 1e-6)
  prot_means <- log10(focal_dist_prot_res$mean_dist_to_ref)
  prot_medians <- log10(focal_dist_prot_res$median_dist_to_ref)
  prot_sds <- log10(focal_dist_prot_res$sd_dist_to_ref)

  zscore_res <-
    data.frame("protein" = focal_dist_prot_res$protein,
               "min_dist_zscores" = calc_zscores(prot_mins),
               "mean_dist_zscores" = calc_zscores(prot_means),
               "median_dist_zscores" = calc_zscores(prot_medians),
               "sd_dist_zscores" = calc_zscores(prot_sds))
  pvals <- lapply(zscore_res[, -1], pnorm)
  qvals <- lapply(pvals, function(x) p.adjust(x, method = "fdr"))
  zscore_res <- as.data.frame(do.call(cbind, c(zscore_res, pvals, qvals)))
  zscore_res[, -1] <- apply(zscore_res[, -1], 2, as.numeric)
  return(zscore_res)
}

# A function to test, for each species, the significance of
# protein similarity to humans
per_spp_wilcox <- function(focal_dist_spp_res) {
  # focal_dist_spp_res: a data.frame that reports, for each species
  # (summarized across all of their respective proteins), the  minimum, mean,
  # median and StdDev of distances to the reference species' (e.g. human)
  # proteins.

  # List to store p-values
  p_values_list <- list()

  # Perform the Wilcoxon Rank Sum test for each species
  for (unique_species in unique(focal_dist_spp_res$species)) {

    # Split the data into two groups: the current species and all other species
    current_species_data <-
      focal_dist_spp_res[focal_dist_spp_res$species == unique_species, ]
    other_species_data <-
      focal_dist_spp_res[focal_dist_spp_res$species != unique_species, ]

    # Perform the Wilcoxon Rank Sum test on minimum distances
    test_min <-
      wilcox.test(current_species_data$min_dist_to_ref,
                  other_species_data$min_dist_to_ref,
                  alternative = "less", exact = FALSE)

    # Perform the Wilcoxon Rank Sum test on mean distances
    test_mean <-
      wilcox.test(current_species_data$mean_dist_to_ref,
                  other_species_data$mean_dist_to_ref,
                  alternative = "less", exact = FALSE)

    # Perform the Wilcoxon Rank Sum test on median distances
    test_median <-
      wilcox.test(current_species_data$median_dist_to_ref,
                  other_species_data$median_dist_to_ref,
                  alternative = "less", exact = FALSE)

    # Check for NA values in standard deviation of distances
    # before performing the test
    if (any(is.na(current_species_data$sd_dist_to_ref)) ||
          any(is.na(other_species_data$sd_dist_to_ref))) {
      test_sd_p_value <- NA
    } else {
      # Perform the Wilcoxon Rank Sum test on standard deviation of distances
      test_sd <-
        wilcox.test(current_species_data$sd_dist_to_ref,
                    other_species_data$sd_dist_to_ref,
                    alternative = "less", exact = FALSE)
      test_sd_p_value <- test_sd$p.value
    }

    # Store the p-values
    p_values_list[[unique_species]] <-
      c(test_min$p.value, test_mean$p.value,
        test_median$p.value, test_sd_p_value)
  }

  # Combine p-values into a data frame
  p_values_df <- do.call(rbind, p_values_list)
  colnames(p_values_df) <-
    c("p_value_min", "p_value_mean", "p_value_median", "p_value_sd")
  q_values_df <- apply(p_values_df, 2, function(x) p.adjust(x, method = "fdr"))
  colnames(q_values_df) <-
    c("q_value_min", "q_value_mean", "q_value_median", "q_value_sd")

  spp_stat_signif <-
    data.frame(species = rownames(p_values_df),
               cbind(p_values_df, q_values_df), row.names = NULL)
  return(spp_stat_signif)
}
