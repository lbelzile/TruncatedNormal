// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

arma::vec ntail(const arma::vec& l, const arma::vec& u) {
  int n = l.n_elem;
  
  if (n == 0) return arma::vec();
  
  arma::vec c = arma::square(l) / 2.0;
  arma::vec f = arma::exp(c - arma::square(u) / 2.0) - 1.0;
  arma::vec x = c - arma::log(1.0 + arma::randu<arma::vec>(n) % f);
  arma::uvec I = find(arma::square(arma::randu<arma::vec>(n)) % x > c);
  
  while (!I.is_empty()) {
    int d = I.n_elem;
    arma::vec cy = c.elem(I);
    arma::vec fy = f.elem(I);
    arma::vec y = cy - arma::log(1.0 + arma::randu<arma::vec>(d) % fy);
    
    arma::vec u2 = arma::randu<arma::vec>(d);
    arma::uvec accepted = find(arma::square(u2) % y < cy);
    
    if (!accepted.is_empty()) {
      x.elem(I.elem(accepted)) = y.elem(accepted);
    }
    
    arma::uvec rejected = find(arma::square(u2) % y >= cy);
    
    if (rejected.is_empty()) break;
    
    I = I.elem(rejected);
  }
  
  return arma::sqrt(2.0 * x);
  
  
}
