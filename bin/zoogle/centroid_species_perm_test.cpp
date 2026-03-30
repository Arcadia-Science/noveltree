#include <Rcpp.h>
#include <random>
#include <numeric>
#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector centroid_species_perm_test(
    NumericVector centroid_dists,
    CharacterVector species_labels,
    int n_permutations) {
  // Input: centroid_dists (distance per protein to family centroid),
  //        species_labels (species name per protein),
  //        n_permutations (number of random shuffles).
  // Process: For each unique species, computes observed mean centroid distance.
  //          Permutes species labels n_permutations times. Two-tailed p-value:
  //          proportion of permutations where |permuted mean - global mean|
  //          >= |observed mean - global mean|.
  // Output: Named NumericVector of p-values, one per unique species.

  int n = centroid_dists.size();

  // Compute global mean centroid distance
  double global_mean = 0.0;
  for (int i = 0; i < n; i++) {
    global_mean += centroid_dists[i];
  }
  global_mean /= n;

  // Map species to their indices and compute observed means
  std::vector<std::string> spp_vec(n);
  for (int i = 0; i < n; i++) {
    spp_vec[i] = Rcpp::as<std::string>(species_labels[i]);
  }

  // Get unique species in order of first appearance
  std::vector<std::string> unique_spp;
  std::unordered_map<std::string, int> spp_index;
  for (int i = 0; i < n; i++) {
    if (spp_index.find(spp_vec[i]) == spp_index.end()) {
      spp_index[spp_vec[i]] = unique_spp.size();
      unique_spp.push_back(spp_vec[i]);
    }
  }
  int n_spp = unique_spp.size();

  // Compute observed mean centroid distance per species
  std::vector<double> obs_means(n_spp, 0.0);
  std::vector<int> spp_counts(n_spp, 0);
  for (int i = 0; i < n; i++) {
    int idx = spp_index[spp_vec[i]];
    obs_means[idx] += centroid_dists[i];
    spp_counts[idx]++;
  }
  for (int s = 0; s < n_spp; s++) {
    obs_means[s] /= spp_counts[s];
  }

  // Compute observed |mean - global_mean| per species
  std::vector<double> obs_dev(n_spp);
  for (int s = 0; s < n_spp; s++) {
    obs_dev[s] = std::abs(obs_means[s] - global_mean);
  }

  // Permutation test
  std::vector<int> perm_count(n_spp, 0);
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::mt19937 rng(std::random_device{}());

  for (int perm = 0; perm < n_permutations; perm++) {
    std::shuffle(indices.begin(), indices.end(), rng);

    // Compute permuted mean centroid distance per species
    // (species labels are shuffled, so protein i gets species of indices[i])
    std::vector<double> perm_means(n_spp, 0.0);
    std::vector<int> perm_counts(n_spp, 0);
    for (int i = 0; i < n; i++) {
      int perm_spp_idx = spp_index[spp_vec[indices[i]]];
      perm_means[perm_spp_idx] += centroid_dists[i];
      perm_counts[perm_spp_idx]++;
    }

    for (int s = 0; s < n_spp; s++) {
      perm_means[s] /= perm_counts[s];
      double perm_dev = std::abs(perm_means[s] - global_mean);
      if (perm_dev >= obs_dev[s]) {
        perm_count[s]++;
      }
    }
  }

  // Compute p-values
  NumericVector p_values(n_spp);
  CharacterVector names(n_spp);
  for (int s = 0; s < n_spp; s++) {
    p_values[s] = (double)perm_count[s] / n_permutations;
    names[s] = unique_spp[s];
  }
  p_values.attr("names") = names;

  return p_values;
}
