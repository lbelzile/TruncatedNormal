// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// lnNpr using arma
inline arma::vec lnNpr_arma(const arma::vec& a, const arma::vec& b, bool check = true) {
  
  if (check) {
    if (a.n_elem != b.n_elem) {
      Rcpp::stop("In function `lnNpr_arma`, vectors `a` and `b` do not have the same size.");
    }
    if (arma::any(a >= b)) {
      Rcpp::stop("In function `lnNpr_arma`, inequality `a < b` not fulfilled for some component.");
    }
  }
  
  int n = a.n_elem;
  arma::vec out(n);
  for (int i = 0; i < n; i++) {
    double pa, pb;
    if (a[i] > 0) {
      pa = R::pnorm(a[i], 0.0, 1.0, 0, 1);
      pb = R::pnorm(b[i], 0.0, 1.0, 0, 1);
      out[i] = pa + log1p(-exp(pb - pa));
    } else if (b[i] < 0) {
      pa = R::pnorm(a[i], 0.0, 1.0, 1, 1);
      pb = R::pnorm(b[i], 0.0, 1.0, 1, 1);
      out[i] = pb + log1p(-exp(pa - pb));
    } else {
      pa = R::pnorm(a[i], 0.0, 1.0, 1, 0);
      pb = R::pnorm(b[i], 0.0, 1.0, 0, 0);
      out[i] = log1p(-pa - pb);
    }
  }
  return out;
}
