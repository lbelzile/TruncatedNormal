// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// RcppParallel
struct MVNWorker : public Worker {
  const arma::mat& L;
  const arma::vec& l;
  const arma::vec& u;
  const arma::vec& mu;
  int d;
  
  RVector<double> output;
  
  MVNWorker(const arma::mat& L_, const arma::vec& l_, const arma::vec& u_,
            const arma::vec& mu_, Rcpp::NumericVector output_)
    : L(L_), l(l_), u(u_), mu(mu_), d(l_.n_elem), output(output_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (size_t j = begin; j < end; j++) {
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
};

// [[Rcpp::export]]
Rcpp::List mvnpr_faster(int n, const arma::mat& L,
                        const arma::vec& l, arma::vec mu,
                        const arma::vec& u,
                        bool useParallel = true,
                        int nThreads = 1) {
  int d = l.n_elem;
  mu[d - 1] = 0.0;
  
  Rcpp::NumericVector weights(n);
  
  if (useParallel) {
    if (nThreads > 1) {
      Rcpp::Function setThreadOptions("setThreadOptions", 
                                      Rcpp::Environment::namespace_env("RcppParallel"));
      setThreadOptions(Rcpp::Named("numThreads") = nThreads);
    }
    MVNWorker worker(L, l, u, mu, weights);
    parallelFor(0, n, worker);
  } else {
    // void mvnpr_serial(int n, const arma::mat& L, const arma::vec& l, const arma::vec& u,
      // const arma::vec& mu, Rcpp::NumericVector& output)
    mvnpr_serial(n, L, l, u, mu, weights);
  }
  
  arma::vec w = Rcpp::as<arma::vec>(weights);
  double prob = arma::mean(w);
  double sd_w = (n > 1) ? std::sqrt(arma::var(w)) : 0.0;
  double relErr = (prob != 0.0) ? (sd_w / std::sqrt((double)n)) / prob : R_PosInf;
  
  return Rcpp::List::create(
    Rcpp::Named("prob") = prob,
    Rcpp::Named("relErr") = relErr
  );
}
