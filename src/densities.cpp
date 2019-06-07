// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
const double log2pi = std::log(2.0 * M_PI);

//' Multivariate normal density function
//'
//' This function returns the log-density for a multivariate Gaussian distribution.
//' The data must be imputed as a matrix, using e.g., \code{as.matrix}, with each row
//' representing an observation.
//'
//' @param x matrix of observations
//' @param mu mean vector
//' @param sigma positive definite covariance matrix
//' @param logd logical; whether log-density should be returned (default to \code{FALSE})
//' @return density or log-density of the \code{nrow(x)} sample
//' @keywords internal
//' @export
// [[Rcpp::export(.dmvnorm_arma)]]
arma::vec dmvnorm_arma(arma::mat x,
                      arma::rowvec mu,
                      arma::mat sigma,
                      bool logd = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mu) ;
    out(i)  = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

//' Multivariate Student density function
//'
//' This function returns the log-density for a multivariate Student distribution.
//' The data must be imputed as a matrix, using e.g., \code{as.matrix}, with each row
//' representing an observation.
//'
//' @param x matrix of observations
//' @param mu location vector
//' @param sigma positive definite scale matrix
//' @param df degrees of freedom
//' @param logd logical; whether log-density should be returned (default to \code{FALSE})
//' @return density or log-density of the \code{nrow(x)} sample
//' @keywords internal
//' @export
// [[Rcpp::export(.dmvt_arma)]]
arma::vec dmvt_arma(arma::mat x,
                    arma::rowvec mu,
                    arma::mat sigma,
                    Rcpp::NumericVector df,
                    bool logd = false) {
  int n = x.n_rows;
  double xdim = static_cast<double>(x.n_cols);
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag())); // log-determinant
  double constants = -0.5 * xdim * (log(df)[0] + std::log(M_PI)) + Rcpp::lgamma(0.5 * (xdim + df))[0] - Rcpp::lgamma(0.5 * df)[0];

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mu) ;
    out(i)      = constants - 0.5 * (df[0] + xdim) * log(1.0 + arma::sum(z%z) / df[0]) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}


