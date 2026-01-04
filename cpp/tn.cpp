// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

arma::vec tn(const arma::vec& l, const arma::vec& u) {
  int n = l.n_elem;
  arma::vec x(n, fill::zeros);
  const double tol = 2.05; // controls switch between methods
  
  arma::uvec wide = find(abs(u - l) > tol);
  
  for (size_t k = 0; k < wide.n_elem; ++k) {
    size_t idx = wide[k];
    double val;
    do {
      val = R::rnorm(0.0, 1.0);
    } while (val < l[idx] || val > u[idx]);
    x[idx] = val;
  }
  
  arma::uvec narrow = find(abs(u - l) <= tol);
  
  for (size_t k = 0; k < narrow.n_elem; ++k) {
    size_t idx = narrow[k];
    double pl = R::pnorm(l[idx], 0.0, 1.0, 1, 0);
    double pu = R::pnorm(u[idx], 0.0, 1.0, 1, 0);
    double u01 = pl + (pu - pl) * R::runif(0.0, 1.0);
    x[idx] = R::qnorm(u01, 0.0, 1.0, 1, 0);
  }
  
  return x;
  
}
