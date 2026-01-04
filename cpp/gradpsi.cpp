// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gradpsi(NumericVector y_nv, const arma::mat &L, NumericVector l_nv, NumericVector u_nv) {
  int d = u_nv.size();
  
  arma::vec y = as<arma::vec>(y_nv);
  arma::vec l = as<arma::vec>(l_nv);
  arma::vec u = as<arma::vec>(u_nv);
  
  arma::vec cv = arma::zeros(d);
  arma::vec x  = arma::zeros(d);
  arma::vec mu = arma::zeros(d);
  
  // x[1:(d-1)] <- y[1:(d-1)]
  x.subvec(0, d - 2) = y.subvec(0, d - 2);
  
  // mu[1:(d-1)] <- y[d:(2*d-2)]
  mu.subvec(0, d - 2) = y.subvec(d - 1, 2 * d - 3);
  
  // cv[-1] <- L[-1,] %*% x
  cv.subvec(1, d - 1) = L.rows(1, d - 1) * x;
  
  arma::vec lt = l - mu - cv;
  arma::vec ut = u - mu - cv;
  
  // call NumericVector lnNpr(NumericVector a, NumericVector b, bool check = true)
  NumericVector lt_nv = wrap(lt);
  NumericVector ut_nv = wrap(ut);
  NumericVector w_nv  = lnNpr(lt_nv, ut_nv, true);
  arma::vec w = as<arma::vec>(w_nv);
  
  
  // compute gradients
  arma::vec pl = arma::exp(-0.5 * arma::square(lt) - w) / std::sqrt(2.0 * M_PI);
  arma::vec pu = arma::exp(-0.5 * arma::square(ut) - w) / std::sqrt(2.0 * M_PI);
  arma::vec P  = pl - pu;
  
  
  arma::vec dfdx = -mu.subvec(0, d - 2) + (P.t() * L.cols(0, d - 2)).t();
  arma::vec dfdm = mu - x + P;
  
  
  arma::vec grad(2 * d - 2);
  grad.subvec(0, d - 2) = dfdx;
  grad.subvec(d - 1, 2 * d - 3) = dfdm.subvec(0, d - 2);
  
  
  // make it a vector
  NumericVector out = wrap(grad);
  out.attr("dim") = R_NilValue;
  return out;
}
