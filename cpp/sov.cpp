// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace arma;
using namespace Rcpp;

inline double Phi(double x) {
  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

inline double PhiInv(double p) {
  return R::qnorm5(p, 0.0, 1.0, 1, 0);
}

void variable_reorder(vec& a0, vec& b0, mat& L) {
  int n = a0.n_elem;
  uvec perm = linspace<uvec>(0, n - 1, n);
  
  for (int i = 0; i < n - 1; ++i) {
    double min_width = std::numeric_limits<double>::infinity();
    int min_j = i;
    
    for (int j = i; j < n; ++j) {
      double Lj = L(j, j);
      double width = Phi(b0(j) / Lj) - Phi(a0(j) / Lj);
      if (width < min_width) {
        min_width = width;
        min_j = j;
      }
    }
    
    if (min_j != i) {
      a0.swap_rows(i, min_j);
      b0.swap_rows(i, min_j);
      L.swap_rows(i, min_j);
      L.swap_cols(i, min_j);
      std::swap(perm(i), perm(min_j));
    }
  }
}

Rcpp::List sov_mvn_prec_lower(const arma::vec& a,
                              const arma::vec& b,
                              const arma::vec& mean,
                              const arma::sp_mat& Q,
                              int N = 10000,
                              bool reorder = false) {
  int n = a.n_elem;
  vec a0 = a - mean;
  vec b0 = b - mean;
  
  mat T;
  if (!chol(T, mat(Q), "lower")) {
    Rcpp::stop("Precision matrix Q is not positive definite");
  }
  
  // no reordering
  
  vec x(n);
  vec weights(N);
  
  for (int k = 0; k < N; ++k) {
    double wprod = 1.0;
    
    for (int i = n - 1; i >= 0; --i) {
      double s = 0.0;
      if (i < n - 1) {
        for (int j = i + 1; j < n; ++j) {
          if (T(j, i) != 0.0) s += T(j, i) * (x(j) - mean(j));
        }
      }
      
      double lb_y = T(i, i) * a0(i) + s;
      double ub_y = T(i, i) * b0(i) + s;
      
      double Phi_lb = Phi(lb_y);
      double Phi_ub = Phi(ub_y);
      double width = Phi_ub - Phi_lb;
      if (width <= 0.0) { wprod = 0.0; break; }
      
      double u = R::runif(0.0, 1.0);
      double y = PhiInv(Phi_lb + u * width);
      x(i) = mean(i) + (y - s) / T(i, i);
      
      wprod *= width;
    }
    weights(k) = wprod;
  }
  
  double est = arma::mean(weights);
  double se  = arma::stddev(weights) / std::sqrt((double)N);
  
  return Rcpp::List::create(
    _["estimate"]  = est,
    _["stderr"]    = se,
    _["reordered"] = false,
    _["sparse"]    = true
  );
}

Rcpp::List sov_mvn_rcpp(const arma::vec& a,
                        const arma::vec& b,
                        const arma::vec& mean,
                        const arma::mat& Sigma,
                        int N = 10000,
                        bool reorder = true) {
  int n = a.n_elem;
  
  vec a0 = a - mean;
  vec b0 = b - mean;
  
  mat L;
  if (!chol(L, Sigma, "lower")) {
    Rcpp::stop("Sigma is not positive definite");
  }
  
  if (reorder) {
    variable_reorder(a0, b0, L);
  }
  
  vec y(n);
  vec weights(N);
  
  for (int k = 0; k < N; ++k) {
    double wprod = 1.0;
    for (int i = 0; i < n; ++i) {
      double sumLy = (i > 0) ? dot(L(i, span(0, i - 1)), y(span(0, i - 1))) : 0.0;
      double ai = (a0(i) - sumLy) / L(i, i);
      double bi = (b0(i) - sumLy) / L(i, i);
      
      double Phi_ai = Phi(ai);
      double Phi_bi = Phi(bi);
      double width = Phi_bi - Phi_ai;
      
      if (width <= 0.0) {
        wprod = 0.0;
        break;
      }
      
      double u = R::runif(0.0, 1.0);
      y(i) = PhiInv(Phi_ai + u * width);
      wprod *= width;
    }
    weights(k) = wprod;
  }
  
  double est = arma::mean(weights);
  double se = arma::stddev(weights) / std::sqrt((double)N);
  
  return Rcpp::List::create(
    _["estimate"]  = est,
    _["stderr"]    = se,
    _["reordered"] = reorder,
    _["sparse"]    = false
  );
}

// [[Rcpp::export]]
Rcpp::List sov_mvn(const arma::vec& a,
                   const arma::vec& b,
                   const arma::vec& mean,
                   const SEXP& mat_or_spmat,
                   int N = 10000,
                   bool reorder = true) {
  
  if (Rf_inherits(mat_or_spmat, "dgCMatrix")) {
    arma::sp_mat Q = Rcpp::as<arma::sp_mat>(mat_or_spmat);
    return sov_mvn_prec_lower(a, b, mean, Q, N, reorder);
  } else {
    arma::mat Sigma = Rcpp::as<arma::mat>(mat_or_spmat);
    return sov_mvn_rcpp(a, b, mean, Sigma, N, reorder);
  }
}
