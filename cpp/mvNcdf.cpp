// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::export]]
double psy_neg_safe(NumericVector par,
                    NumericMatrix L,
                    NumericVector l,
                    NumericVector u) {
  int d = l.size();
  
  NumericVector x(d), mu(d);
  for (int i = 0; i < d - 1; i++) {
    x[i]  = par[i];
    mu[i] = par[d - 1 + i];
  }
  x[d - 1]  = 0.0;
  mu[d - 1] = 0.0;
  
  try {
    return -psy(x, L, l, u, mu);
  } catch (...) {
    return -1e10;
  }
}

// [[Rcpp::export]]
NumericVector heq_mu(NumericVector par,
                     NumericMatrix L_nv,
                     NumericVector l_nv,
                     NumericVector u_nv) {
  
  arma::mat L = Rcpp::as<arma::mat>(L_nv);
  arma::vec l = Rcpp::as<arma::vec>(l_nv);
  arma::vec u = Rcpp::as<arma::vec>(u_nv);
  
  NumericVector g = gradpsi(par, L, l_nv, u_nv);
  
  int d = l.n_elem;
  return g[ Rcpp::Range(d - 1, 2*d - 3) ];
}

// [[Rcpp::export]]
NumericVector hin_feasible(NumericVector par,
                           NumericMatrix L,
                           NumericVector l,
                           NumericVector u,
                           NumericMatrix Lorig,
                           NumericVector lorig,
                           NumericVector uorig) {
  int d = l.size();
  
  NumericVector x(d);
  for (int i = 0; i < d - 1; i++) x[i] = par[i];
  x[d - 1] = 0.0;
  
  NumericVector Lx = Lorig * x;
  
  NumericVector out(2*(d-1));
  for (int i = 0; i < d - 1; i++) {
    out[i]         = uorig[i] - Lx[i];
    out[d - 1 + i] = Lx[i] - lorig[i];
  }
  return out;
}

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
  NumericVector fvec_nv = solvneq["fvec"];
  bool flag = true;
  for (int i = 0; i < fvec_nv.size(); ++i) {
    if (std::abs(fvec_nv[i]) > 1e-6) { flag = false; break; }
  }
  if (!(exitflag == 1 || exitflag == 2)) flag = false;
  
  // if nleqslv failed, try auglag
  if (!flag) {
    Rcpp::Environment alabama = Rcpp::Environment::namespace_env("alabama");
    Rcpp::Function auglag = alabama["auglag"];
    
    Rcpp::List solvneqc = auglag(
      Rcpp::Named("par") = xmu_nv,
      Rcpp::Named("fn")  = Rcpp::InternalFunction(&psy_neg_safe),
      Rcpp::Named("gr")  = Rcpp::InternalFunction(&gradpsi),
      Rcpp::Named("heq") = Rcpp::InternalFunction(&heq_mu),
      Rcpp::Named("hin") = Rcpp::InternalFunction(&hin_feasible),
      Rcpp::Named("L")   = L,
      Rcpp::Named("l")   = lvec,
      Rcpp::Named("u")   = uvec,
      Rcpp::Named("Lorig") = out["L"],    // uses the original L
                                Rcpp::Named("lorig") = l,           // uses the original l
                                Rcpp::Named("uorig") = u,           // uses the original u
                                Rcpp::Named("control.outer") =
                                  Rcpp::List::create(
                                    Rcpp::Named("trace") = false,
                                    Rcpp::Named("method") = "nlminb"
                                  )
    );
    
    if (Rcpp::as<int>(solvneqc["convergence"]) != 0)
      Rcpp::stop("Did not find a solution to the nonlinear system in `mvNcdf`!");
    
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
  arma::vec mu_full = mu;
  mu_full(d-1) = 0.0;
  Rcpp::List est = mvnpr_faster(n, L, lvec, mu_full, uvec, true, numThreads);
  
  // upper bound via psy, since no difference in Cpp and R
  arma::vec x_full = x;
  x_full(d-1) = 0.0;
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
