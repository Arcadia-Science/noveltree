#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector perm_test_across_non_ref(NumericVector distances, int n_permutations) {
  // Input: `distances` (NumericVector of observed distances), `n_permutations` 
  //        (integer number of random shuffles).
  // Process: Conducts a permutation test on each distance, comparing it to 
  //          shuffled distributions to calculate p-values.
  // Output: Returns a NumericVector `p_values`, each representing the probability 
  //         of observing a distance as extreme as the actual value under the null 
  //         hypothesis.
  int n = distances.size();
  NumericVector p_values(n);
  
  // Count the number of times permuted distances are smaller than observed distances
  for (int i = 0; i < n; i++) {
    int count = 0;
    for (int j = 0; j < n_permutations; j++) {
      NumericVector permuted_distances = Rcpp::sample(distances, n, false);
      if (permuted_distances[i] < distances[i]) {
        count++;
      }
    }
    p_values[i] = (double)count / n_permutations;
  }
  
  return p_values;
}
