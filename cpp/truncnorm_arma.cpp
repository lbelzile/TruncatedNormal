// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

arma::vec truncnorm_arma(const arma::vec& l, const arma::vec& u) {
  int n = l.n_elem;
  arma::vec x(n, arma::fill::zeros);
  const double a = 0.4; // threshold for switching between methods
  
  arma::uvec I = find(l > a);
  arma::uvec J = find(u < -a);
  if (!I.is_empty()) x.elem(I) = ntail(l.elem(I), u.elem(I));
  if (!J.is_empty()) x.elem(J) = -ntail(-u.elem(J), -l.elem(J));
  
  arma::uvec mask = arma::zeros<arma::uvec>(n);
  if (!I.is_empty()) mask.elem(I).ones();
  if (!J.is_empty()) mask.elem(J).ones();
  arma::uvec Lset = find(mask == 0);
  if (!Lset.is_empty()) x.elem(Lset) = tn(l.elem(Lset), u.elem(Lset));
  
  return x;
}
