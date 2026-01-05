// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List mvnrnd(int n,
                  const arma::mat& L,
                  arma::vec l,
                  arma::vec u,
                  arma::vec mu) {
  
  int d = l.n_elem;
  
  // last mu element is set to 0
  mu(d - 1) = 0.0;
  
  arma::mat Z(d, n, fill::zeros);
  arma::vec p(n, fill::zeros);    
  
  for (int k = 0; k < d; ++k) {
    
    arma::rowvec col(n, fill::zeros);
    for (int j = 0; j < n; ++j) {          
      double sum = 0.0;
      for (int i = 0; i <= k; ++i) {       
        sum += L(k, i) * Z(i, j);
      }
      col[j] = sum;
    }
    
    // truncation bounds
    arma::vec tl = l(k) - mu(k) - col.t();
    arma::vec tu = u(k) - mu(k) - col.t();
    
    // sample
    arma::vec z(n);
    // call Rcpp::NumericVector trandn_cpp(const arma::vec& l, const arma::vec& u)
    Rcpp::NumericVector sample = trandn_cpp(tl, tu);
    for (int j = 0; j < n; ++j) {
      z[j] = sample[j];
    }
    
    Z.row(k) = (mu(k) + z).t();
    
    // inline arma::vec lnNpr_arma(const arma::vec& a, const arma::vec& b, bool check = true)
    arma::vec lnpr = lnNpr_arma(tl, tu);
    p += lnpr + 0.5 * std::pow(mu(k), 2) - mu(k) * Z.row(k).t();
  }
  
  Rcpp::NumericVector logpr = Rcpp::wrap(p);
  logpr.attr("dim") = R_NilValue;
  
  return Rcpp::List::create(
    Rcpp::Named("logpr") = logpr,
    Rcpp::Named("Z") = Rcpp::wrap(Z)
  );
}
