// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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

// [[Rcpp::export]]
arma::mat jacpsi_cpp(const arma::vec& y, const arma::mat& L, const arma::vec& l, const arma::vec& u) {
  int d = u.n_elem;
  
  // Checks
  // if ((int)l.n_elem != d) stop("l and u must have same length");
  // if ((int)L.n_rows != d || (int)L.n_cols != d) stop("L must be d x d");
  
  // Expect y length 2*d - 2
  // if ((int)y.n_elem != 2*d - 2) stop("y must have length 2*d - 2");
  
  const double SQRT2PI = std::sqrt(2.0 * M_PI);
  arma::vec cv = arma::zeros<arma::vec>(d);
  arma::vec x  = arma::zeros<arma::vec>(d);
  arma::vec mu = arma::zeros<arma::vec>(d);
  
  // fill x[0..d-2] and mu[0..d-2] from y
  if (d > 1) {
    x.subvec(0, d-2) = y.subvec(0, d-2);
    mu.subvec(0, d-2) = y.subvec(d-1, 2*d - 3);
  }
  
  // cv[-1] <- L[-1,] %*% x
  if (d > 1) {
    // rows 1 ... d-1 of cv
    cv.rows(1, d-1) = L.rows(1, d-1) * x;
  }
  //truncated limits
  arma::vec lt = l - mu - cv;
  arma::vec ut = u - mu - cv;
  
  // inline arma::vec lnNpr_arma(const arma::vec& a, const arma::vec& b, bool check = true)
  arma::vec w = lnNpr_arma(lt, ut);
  
  // pl and pu
  arma::vec lt2 = arma::square(lt);
  arma::vec ut2 = arma::square(ut);
  
  arma::vec pl = arma::exp(-0.5 * lt2 - w) / SQRT2PI;
  arma::vec pu = arma::exp(-0.5 * ut2 - w) / SQRT2PI;
  arma::vec P  = pl - pu;
  
  // replace infinite lt/ut with 0
  arma::uvec idx_nonfinite_lt = arma::find_nonfinite(lt);
  if (!idx_nonfinite_lt.is_empty()) lt.elem(idx_nonfinite_lt).zeros();
  
  arma::uvec idx_nonfinite_ut = arma::find_nonfinite(ut);
  if (!idx_nonfinite_ut.is_empty()) ut.elem(idx_nonfinite_ut).zeros();
  
  // dP <- (-P^2) + lt * pl - ut * pu
  arma::vec dP = - arma::square(P) + lt % pl - ut % pu;
  
  // DL = diag(dP) * L
  arma::mat DL = arma::diagmat(dP) * L; // d x d
  
  // mx <- -diag(d) + DL
  arma::mat mx = -arma::eye<arma::mat>(d, d) + DL;
  
  // xx <- crossprod(L, DL) = L' * DL
  arma::mat xx = L.t() * DL;
  
  // (d-1) (d-1)
  arma::uword m = (d > 1) ? (d - 1) : 0;
  arma::mat xx_block, mx_block;
  if (m > 0) {
    xx_block = xx.submat(0, 0, m-1, m-1);
    mx_block = mx.submat(0, 0, m-1, m-1);
  } else {
    // For d==1, not going to happen
    return arma::mat(0,0);
  }
  
  // Build Jacobian:
  // top-left: xx_block 
  // top-right: t(mx_block) 
  // bottom-left: mx_block
  // bottom-right: diag(1 + dP[0:(m-1)])
  arma::mat top = arma::join_horiz(xx_block, mx_block.t());
  arma::mat bottom_right = arma::diagmat( arma::ones<arma::vec>(m) + dP.subvec(0, m-1) );
  arma::mat bottom = arma::join_horiz(mx_block, bottom_right);
  
  arma::mat Jac = arma::join_vert(top, bottom);
  
  return Jac;
}
