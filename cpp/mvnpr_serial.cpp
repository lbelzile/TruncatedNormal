// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// Serial version of MVNWorker loop
void mvnpr_serial(int n, const arma::mat& L,
                  const arma::vec& l, const arma::vec& u,
                  const arma::vec& mu, Rcpp::NumericVector& output) {
  int d = l.n_elem;
  
  for (int j = 0; j < n; j++) {
    arma::vec Z(d, fill::zeros);
    double p = 0.0;
    
    for (int k = 0; k < d - 1; k++) {
      double col = dot(L.row(k).subvec(0, k), Z.subvec(0, k));
      double tl = l[k] - mu[k] - col;
      double tu = u[k] - mu[k] - col;
      
      // arma::vec truncnorm_arma(const arma::vec& l, const arma::vec& u)
      arma::vec sample = truncnorm_arma(arma::vec(1).fill(tl),
                                        arma::vec(1).fill(tu));
      Z[k] = mu[k] + sample[0];
      
      // inline arma::vec lnNpr_arma(const arma::vec& a, const arma::vec& b, bool check = true)
      arma::vec lnpr = lnNpr_arma(arma::vec(1).fill(tl), arma::vec(1).fill(tu));
      p += lnpr[0] + 0.5 * mu[k] * mu[k] - mu[k] * Z[k];
    }
    
    double col_last = dot(L.row(d - 1), Z);
    double tl_last = l[d - 1] - col_last;
    double tu_last = u[d - 1] - col_last;
    
    // inline arma::vec lnNpr_arma(const arma::vec& a, const arma::vec& b, bool check = true)
    arma::vec lnpr_last = lnNpr_arma(arma::vec(1).fill(tl_last),
                                     arma::vec(1).fill(tu_last));
    p += lnpr_last[0];
    
    output[j] = std::exp(p);
  }
}
