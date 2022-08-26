#include <Rcpp.h>
using namespace Rcpp;

// 

// [[Rcpp::export]]
NumericVector tn2(NumericVector l, NumericVector u) {
  //boundary and vector dimension checks
  if(l.size() != u.size()){
    Rcpp::stop("Non-conformal size in vectors l and u");
  }
  
  int n = l.size();
  
  double tol = 2.05;
  NumericVector lowlarge;
  NumericVector lowsmall;
  NumericVector uplarge;
  NumericVector upsmall;
  
  NumericVector ldiffidx;
  NumericVector sdiffidx;
  
  int j = 0;
  int k = 0;
  
  for(int i = 0; i < n; i++){
    if(std::abs(u[i]-l[i]) > tol){
      lowlarge[j] = l[i]; 
      uplarge[j] = u[i];
      ldiffidx[j] = i;
      j++;
    } else {
      lowsmall[k] = l[i]; 
      upsmall[k] = u[i];
      sdiffidx[k] = i;
      k++;
    }
  }
  
  NumericVector xl(ldiffidx.size());
  NumericVector xs(sdiffidx.size());
  
  Function trnd("trnd2");
  xl = trnd(lowlarge, uplarge);
  
  NumericVector pl = R::pnorm(lowsmall, 0.0, 1.0, 1, 0);
  NumericVector pu = R::pnorm(upsmall, 0.0, 1.0, 1, 0);
  for(int i = 0; i < xs.size(); i++){
    xs[i] = R::qnorm(pl[i] + (pu[i]-pl[i])*R::runif(0.0, 1.0));
  }
  
  //need to merge xs and xl back together
  
  NumericVector x(n);
  
  for(int i = 0; i < ldiffidx; i++){
    x[ldiffidx[i]] = xl[i];
  }
  
  for(int i = 0; i < sdiffidx; i++){
    x[sdiffidx[i]] = xs[i];
  }
  return x;
}

