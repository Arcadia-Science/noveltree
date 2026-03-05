#include <RcppArmadillo.h>
using namespace Rcpp;

// A function to, given a matrix of trait observations and a corresponding 
// phylogenetic variance-covariance matrix, transform the data using a 
// phylogenetic GLS transformation using the method described in 
// Butler et al., 2000 - https://doi.org/10.1111/j.0014-3820.2000.tb00026.x
// Note, by default, this assumes a model of brownian motion wherein: 
// "the elements of the G matrix (phylogenetic vcv) are simply the
// amount of time from the root of the phylogeny to the most recent 
// common ancestor of the pair of taxa. The diagonal entries (species variances) 
// are the amounts of time from the root of the tree to each species (the depth 
// of the tree if a molecular clock is assumed)"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List phylo_correction(NumericMatrix trait_matrix, NumericMatrix phylo_vcv) {
  // trait_matrix: matrix of traits (stored in columns) measured for each node (stored in rows) in phylo_vcv
  // phylo_vcv: the phylogenetic covariance matrix (phylo_vcv)
  // Get the phylogenetic variance-covariance matrix
  arma::mat g_matrix = as<arma::mat>(phylo_vcv);
  
  // Compute the Cholesky decomposition of the inverse of the Phylogenetic VCV (g_matrix)
  arma::mat correction_factor;
  try {
    correction_factor = arma::chol(arma::inv(g_matrix));
  } catch (std::runtime_error& e) {
    stop("Error in Cholesky decomposition: ", e.what());
  }
  
  // Get the number of traits and samples
  int num_traits = trait_matrix.ncol();
  int num_samples = trait_matrix.nrow();
  
  // Initialize a list to store the corrected traits
  List U(num_traits);
  
  // Loop through each trait and perform the correction
  for (int t = 0; t < num_traits; t++) {
    
    // Get the current trait
    arma::vec current_trait(num_samples);
    for (int i = 0; i < num_samples; i++) {
      current_trait(i) = trait_matrix(i, t);
    }
    
    // Create a model matrix for the current trait
    arma::mat mod_matrix = arma::ones(num_samples, 1);
    mod_matrix.col(0) = current_trait;
    
    // Perform the correction by post-multiplying the model matrix
    // by the correction factor matrix
    arma::mat corrected_traits = correction_factor * mod_matrix;
    
    // Store the corrected trait in the list
    U[t] = wrap(corrected_traits);
  }
  
  return U;
}
