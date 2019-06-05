#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Calculate log of Gaussian distribution function accurately
//' 
//' This function returns the probability of a standard Gaussian variate between 
//' the interval \code{a} and \code{b}, avoiding numerical overflow. The function is vectorized
//' and is meant to be used only internally by the package \code{TruncatedNormal}.
//' 
//' @param a vector of lower bound
//' @param b vector of upper bound
//' @param check logical; should checks be performed? 
//' @return a vector of log probability.
//' @keywords internal
//' @useDynLib TruncatedNormal, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
NumericVector lnNpr(NumericVector a, NumericVector b, bool check = false){
  if(check){
    //Sanity checks
    if(a.size() != b.size()){
      Rcpp::stop("In function `lnNpr`, vectors `a` and `b` do not have the same size."); 
    }
    if(is_true(any(a >= b))){
      Rcpp::stop("In function `lnNpr`, inequality `a < b` not fulfilled for some component.");  
    }
  }
  //Define containers
  double pa(1);
  double pb(1);
  NumericVector p(a.size());
  for(int i = 0; i < a.size(); i++){
    // case b > a > 0
    // LogicalVector Apos = a > 0;
    if(a[i] > 0){
      // NumericVector pa = Rcpp::pnorm(a[Apos], 0, 1, 0, 1)[0];
      // NumericVector pb = Rcpp::pnorm(b[Apos], 0, 1, 0, 1)[0];
      // p[Apos] = pa + R::log1p(-exp(pb-pa));
      
      pa = R::pnorm(a[i], 0.0, 1.0, 0, 1);
      pb = R::pnorm(b[i], 0.0, 1.0, 0, 1);
      p[i] = pa + log1p(-exp(pb - pa));
    } else if(b[i] < 0){
      // // case a<b<0
      //LogicalVector Bneg = b < 0;
      pa = R::pnorm(a[i], 0.0, 1.0, 1, 1);
      pb = R::pnorm(b[i], 0.0, 1.0, 1, 1);
      p[i] = pb + log1p(-exp(pa - pb));
    } else{
      // case a<0<b
      pa = R::pnorm(a[i], 0.0, 1.0, 1, 0);
      pb = R::pnorm(b[i], 0.0, 1.0, 0, 0);
      p[i] = log1p(-pa - pb);
    }
  }
  return p;
}

// Variance of truncated Gaussian
// 
// The function computes the variance of the standard truncated Gaussian on the
// interval \code{a}, \code{b} with \code{a}< \code{b}.
// 
// @param a vector of lower bound
// @param b vector of upper bound
// @param check logical; should checks be performed? 
// @return vector of variances
NumericVector varTN(NumericVector a, NumericVector b, bool check = false){
  //Sanity checks
  if(a.size() != b.size()){
    Rcpp::stop("In function `varTN`, vectors `a` and `b` do not have the same size."); 
  }
  NumericVector phia = Rcpp::dnorm(a);
  NumericVector phib = Rcpp::dnorm(b);
  NumericVector Phibma = lnNpr(a, b);
  NumericVector varia = 1 + (a*phia - b*phib)/exp(Phibma) - exp(2.0 * (log(abs(phia - phib)) - Phibma));
  return varia;
}

//' Cholesky matrix decomposition with GGE ordering
//' 
//' This function computes the Cholesky decomposition of a covariance matrix
//' \code{Sigma} and returns a list containing the permuted bounds for integration. 
//' The prioritization of the variables follow the rule proposed in Gibson, Glasbey and Elston (1994)
//' and reorder variables to have outermost variables with smallest expected values.
//' 
//' The list contains an integer vector \code{perm} with the indices of the permutation, which is such that
//' \code{Sigma(perm, perm) == L \%*\% t(L)}.
//' The permutation scheme is described in Genz and Bretz (2009) in Section 4.1.3, p.37.
//' @param Sigma \code{d} by \code{d} covariance matrix
//' @param l \code{d} vector of lower bounds
//' @param u \code{d} vector of upper bounds
//' @export
//' @keywords internal
//' @return a list with components
//' \itemize{
//' \item{\code{L}: }{Cholesky root}
//' \item{\code{l}: }{permuted vector of lower bounds}
//' \item{\code{u}: }{permuted vector of upper bounds}
//' \item{\code{perm}: }{vector of integers with ordering of permutation}
//' }
//' @references Genz, A. and Bretz, F. (2009). Computations of Multivariate Normal and t Probabilities, volume 105. Springer, Dordrecht.
//' @references Gibson G.J., Glasbey C.A. and D.A. Elton (1994).  Monte Carlo evaluation of multivariate normal integrals and sensitivity to variate ordering. In: Dimon et al., Advances in Numerical Methods and Applications, WSP, pp. 120-126.
// [[Rcpp::export('.cholpermGB')]]
List cholpermGB(arma::mat Sigma, NumericVector l, NumericVector u){
  if(Sigma.n_cols != l.size() || Sigma.n_cols != u.size()){
    Rcpp::stop("Non conformal size for `l`, `u` and `Sigma`. Check input arguments");
  }
  int d = Sigma.n_cols;
  arma::mat Lc(d, d); //Cholesky matrix
  Lc.zeros(); // Initialize to zero matrix
  if(Sigma.n_rows <= 1){
    Lc(0,0) = pow(Sigma(0,0), 0.5);
    return Rcpp::List::create(Named("L") = Lc, Named("l") = l, Named("u") = u, Named("perm") = IntegerVector::create(1));
    //Rcpp::stop("Matrix must be larger than 1x1");
  }
  double tol = 1.0e-10;
  //Declare containers
  NumericVector a(d);
  NumericVector b(d);
  NumericVector a0(1);
  NumericVector b0(1);
  NumericVector pr0(1);
  arma::vec mu(d);
  IntegerVector perm(d);
  for(int i = 0; i < d; i++){
    a[i] = l[i] / sqrt(Sigma(i,i));
    b[i] = u[i] / sqrt(Sigma(i,i));
    perm[i] = i; //set initial vector with entries
  }
  NumericVector pr = varTN(a, b, false);
  int indmin = which_min(pr); // which_min uses the cpp increment (so returns zero for ordered vectors)
  perm[0] = indmin; perm[indmin] = 0; // swap indices
  // Set Cholesky entries
  double cii = sqrt(Sigma(indmin, indmin));
  Lc.col(0) = Sigma.col(indmin) / cii;
  Lc.swap_rows(0, indmin);
  Lc(0, 0) = cii;
  a0[0] = a[indmin];
  b0[0] = b[indmin];
  pr0 = lnNpr(a0, b0);
  //Expectation of Truncated Normal on (a, b)
  mu(0) = (exp(-0.5 * pow(a[indmin], 2) - pr0[0]) - exp(-0.5 * pow(b[indmin], 2) - pr0[0]))/pow(2.0 * M_PI, 0.5);
  // END OF LOOP FOR FIRST ITERATION
  double mui; double denomi;
  if(d > 1){
  for(int j = 1; j < d; j++){
    NumericVector a = NumericVector(d - j);
    NumericVector b = NumericVector(d - j);
    for(int i0 = j; i0 < d; i0++){
      int i = perm[i0];
      mui = dot(Lc.submat(i0, 0, i0, j - 1), mu.subvec(0, j - 1));
      denomi = Sigma(i, i) - dot(Lc.submat(i0, 0, i0, j - 1), Lc.submat(i0, 0, i0, j - 1));
      if(denomi < -0.001){
        Rcpp::stop("`Sigma` is not positive definite");
      } else if(denomi < 0){
        denomi = tol;
      } else{
        denomi = pow(denomi, 0.5); 
      }
      a[i0 - j] = (l[i] - mui) / denomi;
      b[i0 - j] = (u[i] - mui) / denomi;
    }
    NumericVector pr = varTN(a, b);
    int min0 = which_min(pr);
    int indmin = perm[min0 + j]; // index of minimum amongst remaining entries
    if(min0 > 0){
      Lc.swap_rows(min0 + j, j);
      perm[min0 + j] = perm[j]; // swap indices in permutation
      perm[j] = indmin;
    }
    //Update Cholesky, jth column
    Lc(j, j) = sqrt(Sigma(indmin, indmin) - dot(Lc.submat(j, 0, j, j - 1), Lc.submat(j, 0, j, j - 1)));
    if(j < (d - 1)){ // non-empty loop
      for(int i0 = j + 1; i0 < d; i0++){
        int i = perm[i0];
        Lc(i0,j) = (Sigma(i, indmin)  - dot(Lc.submat(j, 0, j, j - 1), Lc.submat(i0, 0, i0, j - 1)))/ Lc(j, j);
      }
    }
    a0[0] = a[min0];
    b0[0] = b[min0];
    pr0 = lnNpr(a0, b0);
    // Compute E(a,b)
    mu(j) = (exp(-0.5 * pow(a0[0], 2) - pr0[0]) - exp(-0.5 * pow(b0[0], 2) - pr0[0]))/pow(2.0 * M_PI, 0.5);
  }
  }
  return Rcpp::List::create(Named("L") = Lc, Named("l") = l[perm], Named("u") = u[perm], Named("perm") = perm + IntegerVector(d, 1));
}


//' Cholesky matrix decomposition with GGE ordering
//' 
//' This function computes the Cholesky decomposition of a covariance matrix
//' \code{Sigma} and returns a list containing the permuted bounds for integration. 
//' The prioritization of the variables follow the rule proposed in Gibson, Glasbey and Elston (1994)
//' and reorder variables to have outermost variables with smallest expected values.
//' 
//' The list contains an integer vector \code{perm} with the indices of the permutation, which is such that
//' \code{Sigma(perm, perm) == L \%*\% t(L)}.
//' The permutation scheme is described in Genz and Bretz (2009) in Section 4.1.3, p.37.
//' @param Sigma \code{d} by \code{d} covariance matrix
//' @param l \code{d} vector of lower bounds
//' @param u \code{d} vector of upper bounds
//' @export
//' @keywords internal
//' @return a list with components
//' \itemize{
//' \item{\code{L}: }{Cholesky root}
//' \item{\code{l}: }{permuted vector of lower bounds}
//' \item{\code{u}: }{permuted vector of upper bounds}
//' \item{\code{perm}: }{vector of integers with ordering of permutation}
//' }
//' @references Genz, A. and Bretz, F. (2009). Computations of Multivariate Normal and t Probabilities, volume 105. Springer, Dordrecht.
//' @references Gibson G.J., Glasbey C.A. and D.A. Elton (1994).  Monte Carlo evaluation of multivariate normal integrals and sensitivity to variate ordering. In: Dimon et al., Advances in Numerical Methods and Applications, WSP, pp. 120-126.
// [[Rcpp::export('.cholpermGGE')]]
List cholperm(arma::mat Sigma, NumericVector l, NumericVector u){
  if(Sigma.n_cols != l.size() || Sigma.n_cols != u.size()){
    Rcpp::stop("Non conformal size for `l`, `u` and `Sigma`. Check input arguments");
  }
  int d = Sigma.n_cols;
  arma::mat Lc(d, d); //Cholesky matrix
  Lc.zeros(); // Initialize to zero matrix
  if(Sigma.n_rows <= 1){
    Lc(0,0) = pow(Sigma(0,0), 0.5);
    return Rcpp::List::create(Named("L") = Lc, Named("l") = l, Named("u") = u, Named("perm") = IntegerVector::create(1));
    //Rcpp::stop("Matrix must be larger than 1x1");
  }
  double tol = 1.0e-10;
  //Declare containers
 
  NumericVector a(d);
  NumericVector b(d);
  arma::vec mu(d);
  IntegerVector perm(d);

  for(int i = 0; i < d; i++){
    a[i] = l[i] / sqrt(Sigma(i,i));
    b[i] = u[i] / sqrt(Sigma(i,i));
    perm[i] = i; //set initial vector with entries
  }
  NumericVector pr = lnNpr(a, b);
  int indmin = which_min(pr); // which_min uses the cpp increment (so returns zero for ordered vectors)
  perm[0] = indmin; perm[indmin] = 0; // swap indices
  // Set Cholesky entries
  double cii = sqrt(Sigma(indmin, indmin));
  Lc.col(0) = Sigma.col(indmin) / cii;
  Lc.swap_rows(0, indmin);
  Lc(0, 0) = cii;
  //Expectation of Truncated Normal on (a, b)
  mu(0) = (exp(-0.5 * pow(a[indmin], 2) - pr[indmin]) - exp(-0.5 * pow(b[indmin], 2) - pr[indmin]))/pow(2.0 * M_PI, 0.5);
  // END OF LOOP FOR FIRST ITERATION
  double mui; double denomi;
  for(int j = 1; j < d; j++){
    NumericVector a = NumericVector(d - j);
    NumericVector b = NumericVector(d - j);
    for(int i0 = j; i0 < d; i0++){
      int i = perm[i0];
      mui = dot(Lc.submat(i0, 0, i0, j - 1), mu.subvec(0, j - 1));
      denomi = Sigma(i, i) - dot(Lc.submat(i0, 0, i0, j - 1), Lc.submat(i0, 0, i0, j - 1));
      if(denomi < -0.001){
        Rcpp::stop("`Sigma` is not positive definite");
      } else if(denomi < 0){
        denomi = tol;
      } else{
       denomi = pow(denomi, 0.5); 
      }
      a[i0 - j] = (l[i] - mui) / denomi;
      b[i0 - j] = (u[i] - mui) / denomi;
    }
    NumericVector pr = lnNpr(a, b);
    int min0 = which_min(pr);
    int indmin = perm[min0 + j]; // index of minimum amongst remaining entries
    if(min0 > 0){
      Lc.swap_rows(min0 + j, j);
      perm[min0 + j] = perm[j]; // swap indices in permutation
      perm[j] = indmin;
    }
    //Update Cholesky, jth column
    Lc(j, j) = sqrt(Sigma(indmin, indmin) - dot(Lc.submat(j, 0, j, j - 1), Lc.submat(j, 0, j, j - 1)));
    if(j < (d - 1)){ // non-empty loop
      for(int i0 = j + 1; i0 < d; i0++){
        int i = perm[i0];
        Lc(i0,j) = (Sigma(i, indmin)  - dot(Lc.submat(j, 0, j, j - 1), Lc.submat(i0, 0, i0, j - 1)))/ Lc(j, j);
      }
    }
    // Compute E(a,b)
    mu(j) = (exp(-0.5 * pow(a[min0], 2) - pr[min0]) - exp(-0.5 * pow(b[min0], 2) - pr[min0]))/pow(2.0 * M_PI, 0.5);
  }
  return Rcpp::List::create(Named("L") = Lc, Named("l") = l[perm], Named("u") = u[perm], Named("perm") = perm + IntegerVector(d, 1));
}

//' Quantile of truncated Gaussian
//' 
//' The function compute the \code{p}th quantile associated to the truncated standard Gaussian 
//' variate on the interval (\code{l},\code{u}).
//' 
//' @param p vector of probabilities
//' @param l \code{d} vector of lower bounds
//' @param u \code{d} vector of upper bounds
//' @return vector of quantiles
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector Phinv(NumericVector p, NumericVector l, NumericVector u){
  if(p.size() != l.size() || p.size() != u.size()){
    Rcpp::stop("Non-conformal sizes in `p`, `l` or `u` in function `Phinv`");
  }
  LogicalVector flip  = LogicalVector(u < 0);
     for(int i = 0; i < u.size(); i++){
       if(flip[i] == true){
         l[i] = -l[i];
         u[i] = -u[i]; 
       }
     }
  //use symmetry of normal
    NumericVector pl = Rcpp::pnorm(l, 0.0, 1.0, 0, 0);
    NumericVector x = Rcpp::qnorm(pl + (Rcpp::pnorm(u, 0.0, 1.0, 0, 0) - pl) * p, 0.0, 1.0, 0, 0);
     for(int i = 0; i < u.size(); i++){
       if(flip[i] == true){
        x[i] = -x[i];
     }
    }
    return x;
  }
