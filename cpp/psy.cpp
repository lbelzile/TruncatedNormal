// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double psy(NumericVector x, NumericMatrix L, NumericVector l, NumericVector u, NumericVector mu) {
  
  int d = u.size();
  
  // last element = 0
  x[d-1]  = 0.0;
  mu[d-1] = 0.0;
  
  // convert to Armadillo views
  arma::vec x_ = as<arma::vec>(x);
  arma::mat L_ = as<arma::mat>(L);
  arma::vec l_ = as<arma::vec>(l);
  arma::vec u_ = as<arma::vec>(u);
  arma::vec mu_ = as<arma::vec>(mu);
  
  // cv = L * x
  arma::vec cv = L_ * x_;
  
  // l - mu - cv, u - mu - cv
  NumericVector lower = wrap(l_ - mu_ - cv);
  NumericVector upper = wrap(u_ - mu_ - cv);
  
  NumericVector logprobs = lnNpr(lower, upper, false);
  
  arma::vec lp = as<arma::vec>(logprobs);
  
  arma::vec terms = lp + 0.5 * square(mu_) - (x_ % mu_);
  
  //sum
  return accu(terms);
}
