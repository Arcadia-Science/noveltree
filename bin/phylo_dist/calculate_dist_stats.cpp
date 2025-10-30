#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List calculate_dist_stats(arma::mat focal_dists)
{
  // The following Rcpp function efficiently calculates summary stats
  // (i.e. minimum, max, mean, etc) of each protein to the corresponding
  // proteins in the "reference" species (i.e. humans).
  int n = focal_dists.n_rows;
  arma::vec focal_dist_mins(n);
  arma::vec focal_dist_means(n);
  arma::vec focal_dist_medians(n);
  arma::vec focal_dist_sds(n);

  for (int i = 0; i < n; ++i)
  {
    arma::vec row = focal_dists.row(i).t();          // transpose the row to get a column vector
    row.replace(arma::datum::nan, arma::datum::inf); // Replace NaN with Inf for arma::find_finite

    arma::uvec finite_indices = arma::find_finite(row); // Get indices of finite values
    arma::vec finite_values = row.elem(finite_indices); // Get finite values

    if (finite_values.n_elem > 0)
    {
      focal_dist_mins(i) = arma::min(finite_values);
      focal_dist_means(i) = arma::mean(finite_values);
      focal_dist_medians(i) = arma::median(finite_values);
      if (finite_values.n_elem > 1)
      {
        focal_dist_sds(i) = arma::stddev(finite_values);
      }
      else
      {
        focal_dist_sds(i) = arma::datum::nan;
      }
    }
    else
    {
      focal_dist_mins(i) = arma::datum::nan;
      focal_dist_means(i) = arma::datum::nan;
      focal_dist_medians(i) = arma::datum::nan;
      focal_dist_sds(i) = arma::datum::nan;
    }
  }

  return List::create(
      Named("min") = focal_dist_mins,
      Named("mean") = focal_dist_means,
      Named("median") = focal_dist_medians,
      Named("sd") = focal_dist_sds);
}
