#include <Rcpp.h>
#include <random>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector perm_test_across_non_ref(NumericVector distances, int n_permutations) {
  // Input: `distances` (NumericVector of observed distances), `n_permutations`
  //        (integer number of random shuffles).
  // Process: Conducts a permutation test on each distance by shuffling an index
  //          vector (shared across all elements per permutation) and comparing
  //          permuted distances to observed values.
  // Output: Returns a NumericVector `p_values`, each representing the probability
  //         of observing a distance as extreme as the actual value under the null
  //         hypothesis.
  int n = distances.size();
  NumericVector p_values(n);
  IntegerVector count(n);

  // Build index vector once, reuse across all permutations
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::mt19937 rng(std::random_device{}());

  // Shared shuffle: one permutation evaluated across all elements
  for (int perm = 0; perm < n_permutations; perm++) {
    std::shuffle(indices.begin(), indices.end(), rng);
    for (int i = 0; i < n; i++) {
      if (distances[indices[i]] < distances[i]) {
        count[i]++;
      }
    }
  }

  for (int i = 0; i < n; i++) {
    p_values[i] = (double)count[i] / n_permutations;
  }

  return p_values;
}
