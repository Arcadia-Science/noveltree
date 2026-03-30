#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat compute_inverse_covariance(const arma::mat& transformed_data) {
  // transformed_data: matrix with traits/phenotypes stored in columns, and species/nodes in rows
  // Regularize to ensure the matrix is positive definite and compute the corrected VCV
  arma::mat corrected_vcv = arma::cov(transformed_data) + arma::eye(transformed_data.n_cols, transformed_data.n_cols) * 1e-6;
  
  // Compute the inverse of the corrected VCV
  arma::mat inv_cov = arma::inv(corrected_vcv);
  return inv_cov;
}
