// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List call_nleqslv(NumericVector x0,
                        Function fn,
                        Function jac,
                        SEXP L,
                        SEXP l,
                        SEXP u) {
  // Access nleqslv
  Environment pkg = Environment::namespace_env("nleqslv");
  Function nleqslv = pkg["nleqslv"];
  
  Function list("list");
  List control = list(_["maxit"] = 500L);
  
  List res = nleqslv(
    x0,
    Named("fn") = fn,
    Named("jac") = jac,
    Named("L") = L,
    Named("l") = l,
    Named("u") = u,
    Named("global") = "pwldog",
    Named("method") = "Broyden",
    Named("control") = control
  );
  return res;
}


// [[Rcpp::export]]
SEXP mvrandn_cpp(NumericVector l_in,
                 NumericVector u_in,
                 NumericMatrix Sig_in,
                 int n,
                 Nullable<NumericVector> mu_in = R_NilValue) {
  
  int d = l_in.size();
  
  if (u_in.size() != d) stop("l, u, and Sig have to match in dimension with u>l");
  
  if ((int)Sig_in.nrow() != d || (int)Sig_in.ncol() != d) stop("l, u, and Sig must match in dimension");
  
  arma::mat Sig = as<arma::mat>(Sig_in);
  
  NumericVector l = clone(l_in);
  NumericVector u = clone(u_in);
  
  NumericVector mu;
  bool has_mu = mu_in.isNotNull();
  if (has_mu) {
    mu = clone(mu_in);
    
    for (int i=0;i<d;i++){
      l[i] = l[i] - mu[i];
      u[i] = u[i] - mu[i];
    }
  }
  
  // d == 1
  if (d == 1) {
    double stddev = std::sqrt(Sig(0,0));
    
    //  expand bounds to length n
    arma::vec lstd(n, arma::fill::value(l[0] / stddev));
    arma::vec ustd(n, arma::fill::value(u[0] / stddev));
    
    // draw n samples
    // Rcpp::NumericVector trandn_cpp(const arma::vec& l, const arma::vec& u)
    Rcpp::NumericVector samp = trandn_cpp(lstd, ustd);
    
    // scale back and add mu
    NumericMatrix out(1, n);
    for (int i = 0; i < n; ++i) {
      double val = stddev * samp[i];
      if (has_mu) val += mu[0];
      out(0, i) = val;
    }
    
    Rcpp::NumericVector out_res = Rcpp::wrap(out);
    out_res.attr("dim") = R_NilValue;  // make it into a vector
    return out_res;
  }
  
  // using cholperm
  // List cholperm(arma::mat Sigma, NumericVector l, NumericVector u)
  List out = cholperm(Sig, l, u);
  
  
  arma::mat Lfull = as<arma::mat>(out["L"]);
  arma::vec lvec = as<arma::vec>(out["l"]);
  arma::vec uvec = as<arma::vec>(out["u"]);
  IntegerVector perm1based = out["perm"];
  
  int dd = d;
  arma::uvec perm_zero(dd);
  for (int i = 0; i < dd; ++i)
    perm_zero(i) = (unsigned int)(perm1based[i] - 1);
  
  // Diagonal D
  arma::vec D = Lfull.diag();
  
  // Warn if nearly singular
  if (arma::any(D < 1e-10))
    warning("Method may fail as covariance matrix is singular!");
  
  // L = Lfull / D (divide each row)
  arma::mat L = Lfull;
  for (int i = 0; i < dd; ++i)
    L.row(i) /= D(i);
  
  // L <- L - diag(d)
  L = L - arma::eye<arma::mat>(L.n_rows, L.n_cols);
  
  // Scale l and u
  lvec /= D;
  uvec /= D;
  
  
  int m = 2 * dd - 2;
  NumericVector x0(m);
  for (int i=0;i<m;i++) x0[i] = 0.0;
  
  
  Function gradpsi_R = Rcpp::Environment::global_env()["gradpsi"];
  Function jacpsi_R  = Rcpp::Environment::global_env()["jacpsi"];
  
  
  List solvneq = call_nleqslv(x0, gradpsi_R, jacpsi_R, wrap(L), wrap(lvec), wrap(uvec));
  NumericVector xmu = solvneq["x"];
  int exitflag = as<int>(solvneq["termcd"]);
  
  if (!(exitflag == 1 || exitflag == 2)) {
    if (solvneq.containsElementNamed("fvec")) {
      NumericVector fvec = solvneq["fvec"];
      bool allzero = true;
      for (int i=0;i<fvec.size();i++){
        if (std::abs(fvec[i]) > 1e-6) { allzero = false; break; }
      }
      if (!allzero) warning("Did not find a solution to the nonlinear system in `mvrandn`!");
    } else {
      warning("Did not find a solution to the nonlinear system in `mvrandn`!");
    }
  }
  
  
  arma::vec x = arma::zeros(dd);
  arma::vec muV = arma::zeros(dd);
  
  //sets last element to 0
  for (int i=0;i<dd-1;i++){
    x(i) = xmu[i];
    muV(i) = xmu[dd - 1 + i];
  }
  
  
  // call psy(x, L, l, u, muV)
  // double psy(NumericVector x, NumericMatrix L, NumericVector l, NumericVector u, NumericVector mu)
  double psistar = psy(wrap(x), wrap(L), wrap(lvec), wrap(uvec), wrap(muV));
  // there is an extra zero at the end of x and muV in Cpp version when passed to psy()
  // R version does not
  // but in psy(), the x[d], x[muV] version are set to zero again, so no harm done
  
  // output matrix Z
  arma::mat Z = arma::zeros<arma::mat>(dd, n);
  
  int accept = 0;
  int iter = 0;
  int nsim = n;
  long ntotsim = 0L;
  
  while (accept < n) {
    
    // Rcpp::List mvnrnd(int n, const arma::mat& L, arma::vec l, arma::vec u, arma::vec mu)
    List call = mvnrnd(nsim, L, lvec, uvec, muV);
    ntotsim += nsim;
    arma::vec re = as<arma::vec>(Rcpp::rexp(nsim));
    
    arma::vec logpr = as<arma::vec>(call["logpr"]);
    arma::mat Zcall = as<arma::mat>(call["Z"]);
    
    arma::vec crit = psistar - logpr;
    arma::uvec idx = arma::find(re > crit);
    
    int m_accept = idx.n_elem;
    
    if (m_accept > n - accept) {
      m_accept = n - accept;
      idx = idx.subvec(0, m_accept - 1);
    }
    if (m_accept > 0) {
      Z.cols(accept, accept + m_accept - 1) = Zcall.cols(idx);
    }
    
    accept += m_accept;
    iter += 1;
    
    if (m_accept > 0) {
      nsim = std::min((int)1e6,
                      std::min(n, (int)std::ceil((double)nsim / (double)m_accept)));
    } else {
      nsim = std::min(1000000, n);
    }
    
    
    if ((ntotsim > 10000) && ((double)accept / (double)ntotsim < 0.001)) {
      warning("Acceptance probability smaller than 0.001");
    } else if (iter > 100000) {
      if (accept == 0) {
        stop("Could not sample from truncated Normal - check input");
      } else if (accept > 1) {
        Z = Z.cols(0, accept - 1);
        warning("Sample of size smaller than n returned.");
        break;
      }
    }
  } 
  
  
  arma::uvec order(dd);
  for (int i = 0; i < dd; i++) {
    order(perm_zero(i)) = i;
  }
  arma::mat finalZ = Lfull * Z;
  
  
  arma::mat finalZ_reordered(dd, finalZ.n_cols);
  for (unsigned i = 0; i < (unsigned)dd; ++i) {
    finalZ_reordered.row(i) = finalZ.row(order(i));
  }
  
  
  // add mu back
  NumericMatrix ret(dd, finalZ_reordered.n_cols);
  for (unsigned j = 0; j < finalZ_reordered.n_cols; ++j) {
    for (unsigned i = 0; i < (unsigned)dd; ++i) {
      double val = finalZ_reordered(i, j);
      if (has_mu) val += mu[i];
      ret(i, j) = val;
    }
  }
  return ret;
}
