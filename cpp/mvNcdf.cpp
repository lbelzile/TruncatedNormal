// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::export]]
Rcpp::List mvNcdf_cpp(Rcpp::NumericVector l, 
                      Rcpp::NumericVector u, 
                      Rcpp::NumericMatrix Sig, 
                      int n = 100000,
                      int numThreads = 1) {
  int d = l.size();
  if (u.size() != d || Sig.nrow() != d || Sig.ncol() != d)
    Rcpp::stop("l, u, and Sig must have matching dimensions.");
  
  // 1D case
  if (d == 1) {
    Rcpp::Function lnNpr("lnNpr");
    double s = std::sqrt(Sig(0,0));
    double logp = Rcpp::as<double>(lnNpr(l/s, u/s, true));
    return Rcpp::List::create(
      Rcpp::Named("prob") = std::exp(logp),
      Rcpp::Named("err") = NA_REAL,
      Rcpp::Named("relErr") = NA_REAL,
      Rcpp::Named("upbnd") = NA_REAL
    );
  }
  
  // Cholesky + permutation
  // // [[Rcpp::export('.cholpermGGE')]]
  // List cholperm(arma::mat Sigma, NumericVector l, NumericVector u)
  Rcpp::List out = cholperm(arma::mat(Sig.begin(), d, d, false), l, u);
  arma::mat L = Rcpp::as<arma::mat>(out["L"]);
  arma::vec lvec = Rcpp::as<arma::vec>(out["l"]);
  arma::vec uvec = Rcpp::as<arma::vec>(out["u"]);
  
  // scale L, l, u
  arma::vec D = L.diag();
  for (int i=0; i<d; i++) {
    L.row(i) /= D(i);
    lvec(i) /= D(i);
    uvec(i) /= D(i);
  }
  L -= arma::eye(d,d);
  
  // initial guess
  Rcpp::NumericVector x0(2*d - 2, 0.0);
  
  // first try nleqslv
  Rcpp::Environment nleqslv = Rcpp::Environment::namespace_env("nleqslv");
  Rcpp::Function nleqslv_fn = nleqslv["nleqslv"];
  
  Rcpp::List control = Rcpp::List::create(Rcpp::Named("maxit") = 500L);
  Rcpp::List solvneq = nleqslv_fn(x0,
                                  Rcpp::Named("fn") = Rcpp::InternalFunction(&gradpsi),
                                  Rcpp::Named("jac") = Rcpp::InternalFunction(&jacpsi),
                                  Rcpp::Named("L") = L,
                                  Rcpp::Named("l") = lvec,
                                  Rcpp::Named("u") = uvec,
                                  Rcpp::Named("global") = "pwldog",
                                  Rcpp::Named("method") = "Broyden",
                                  Rcpp::Named("control") = control);
  
  Rcpp::NumericVector xmu_nv = solvneq["x"];
  int exitflag = solvneq["termcd"];
  
  bool ok = (exitflag == 1 || exitflag == 2);
  
  // if nleqslv failed, try auglag
  if (!ok) {
    Rcpp::Environment alabama = Rcpp::Environment::namespace_env("alabama");
    Rcpp::Function auglag = alabama["auglag"];
    Rcpp::List solvneqc = auglag(x0,
                                 Rcpp::Named("fn") = Rcpp::InternalFunction(&psy),
                                 Rcpp::Named("gr") = Rcpp::InternalFunction(&gradpsi),
                                 Rcpp::Named("heq") = Rcpp::InternalFunction(&gradpsi),
                                 Rcpp::Named("hin") = R_NilValue,
                                 Rcpp::Named("control.outer") = Rcpp::List::create(Rcpp::Named("trace")=false));
    xmu_nv = solvneqc["par"];
  }
  
  // parse x and mu
  arma::vec xmu = Rcpp::as<arma::vec>(xmu_nv);
  arma::vec x(d, arma::fill::zeros);
  arma::vec mu(d, arma::fill::zeros);
  for (int i=0; i<d-1; i++) {
    x(i) = xmu[i];
    mu(i) = xmu[d-1+i];
  }
  
  // call mvnpr_faster
  // Rcpp::List mvnpr_faster(int n, const arma::mat& L,
  // const arma::vec& l, arma::vec mu,
  // const arma::vec& u,
  // bool useParallel = true,
  // int nThreads = 1)
  arma::vec mu_full = mu;
  mu_full(d-1) = 0.0;
  Rcpp::List est = mvnpr_faster(n, L, lvec, mu_full, uvec, true, numThreads);
  
  // upper bound using psy
  arma::vec x_full = x;
  x_full(d-1) = 0.0;
  // double psy(NumericVector x, NumericMatrix L, NumericVector l, NumericVector u, NumericVector mu)
  double upbnd = std::exp(psy(Rcpp::wrap(x_full),
                              Rcpp::wrap(L),
                              Rcpp::wrap(lvec),
                              Rcpp::wrap(uvec),
                              Rcpp::wrap(mu_full)));
  
  return Rcpp::List::create(
    Rcpp::Named("prob") = est["prob"],
    Rcpp::Named("relErr") = est["relErr"],
    Rcpp::Named("upbnd") = upbnd
  );
}
