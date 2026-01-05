// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericVector trandn_cpp(const arma::vec& l, const arma::vec& u) {
  
  int n = l.n_elem;
  arma::vec x(n, arma::fill::zeros);
  const double a = 0.4;
  
  arma::uvec I = arma::find(l > a);
  arma::uvec J = arma::find(u < -a);
  
  // call arma::vec ntail(const arma::vec& l, const arma::vec& u)
  if (!I.is_empty()) x.elem(I) = ntail(l.elem(I), u.elem(I));
  if (!J.is_empty()) x.elem(J) = -ntail(-u.elem(J), -l.elem(J));
  
  arma::uvec mask = arma::zeros<arma::uvec>(n);
  if (!I.is_empty()) mask.elem(I).ones();
  if (!J.is_empty()) mask.elem(J).ones();
  arma::uvec L = arma::find(mask == 0);
  
  // call arma::vec tn(const arma::vec& l, const arma::vec& u)
  if (!L.is_empty()) x.elem(L) = tn(l.elem(L), u.elem(L));
  
  // wrap into NumericVector
  Rcpp::NumericVector out = Rcpp::wrap(x);
  out.attr("dim") = R_NilValue;  // make it into a vector
  return out;
}
