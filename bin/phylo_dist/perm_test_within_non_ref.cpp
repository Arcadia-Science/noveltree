#include <Rcpp.h>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix perm_test_within_non_ref(NumericMatrix dist_mat, NumericMatrix focal_dists, int n_permutations) {
  // Input: Takes a full distance matrix (dist_mat) between all proteins and a
  //        matrix (focal_dists) of distances from non-human to human proteins,
  //        along with the number of permutations (n_permutations).
  // Process: Conducts a permutation test by shuffling non-diagonal elements in
  //          dist_mat and counts occurrences where the permuted distances are
  //          less than or equal to the distances in focal_dists.
  // Output: Returns a matrix of p-values (pvals), representing the likelihood of
  //         each non-human to human protein distance being as extreme as observed
  //         under the null hypothesis.

  int n_rows = focal_dists.nrow();
  int n_cols = focal_dists.ncol();
  NumericMatrix count_mat(n_rows, n_cols);

  // Find indices of reference (e.g. human) proteins in dist_mat
  IntegerVector ref_indices(n_cols);
  CharacterVector dist_mat_colnames = colnames(dist_mat);
  CharacterVector focal_dists_colnames = colnames(focal_dists);

  for (int i = 0; i < n_cols; i++) {
    std::string ref_protein = Rcpp::as<std::string>(focal_dists_colnames[i]);
    for (int j = 0; j < dist_mat_colnames.size(); j++) {
      if (ref_protein == Rcpp::as<std::string>(dist_mat_colnames[j])) {
        ref_indices[i] = j;
        break;
      }
    }
  }

  // Permutation test
  for (int perm = 0; perm < n_permutations; perm++) {
    for (int i = 0; i < n_rows; i++) {
      NumericVector rowVec(dist_mat.nrow());
      // Copy the non-diagonal elements to a separate vector
      for (int j = 0; j < rowVec.size(); j++) {
        if (i != j) {
          rowVec[j] = dist_mat(i, j);
        } else {
          rowVec[j] = NA_REAL;
        }
      }
      std::shuffle(rowVec.begin(), rowVec.end(), std::default_random_engine(std::rand()));

      // Counting the number of times permuted distances are smaller than observed distances
      for (int j = 0; j < n_cols; j++) {
        int ref_col = ref_indices[j];
        if (!NumericVector::is_na(rowVec[ref_col]) && rowVec[ref_col] < focal_dists(i, j)) {
          count_mat(i, j)++;
        }
      }
    }
  }

  // Calculate p-values
  NumericMatrix pvals = count_mat / n_permutations;
  rownames(pvals) = rownames(focal_dists);
  colnames(pvals) = colnames(focal_dists);

  return pvals;
}
