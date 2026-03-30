#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat compute_chunk_distances(const arma::mat& data, const arma::mat& inv_cov, const arma::umat& pair_indices) {
  // This function calculates the actual mahalanobis distances between the
  // set of observating within a matrix of paired indices (i.e. a chunk) given 
  // the phenotypic data, the phenotypic covariance matrix, and the matrix 
  // of indices among which we'd like to make the calculations. 
  
  // data: trait data as class matrix with observations in rows, traits stored in columns
  // inv_cov: inverse of the phenotypic covariance matrix. # of rows and columns equal to number of traits
  // pair_indices: list containing the indices of each pairwise comparison of observations we will calculate distance between
  
  int num_pairs = pair_indices.n_cols;
  arma::mat distances(num_pairs, 1);
  
  for(int p = 0; p < num_pairs; p++) {
    int i = pair_indices(0, p) - 1;  // Adjust for 0-based indexing in C++
    int j = pair_indices(1, p) - 1;  // Adjust for 0-based indexing in C++
    
    // Compute the difference between the two observations
    arma::rowvec diff = data.row(i) - data.row(j);
    
    // Ensure the matrix operation is correct
    arma::mat temp = diff * inv_cov * arma::trans(diff);  // This should result in a 1x1 matrix
    distances(p) = std::sqrt(temp(0,0));  // Extract the scalar value from the 1x1 matrix
  }
  
  return distances;
}
