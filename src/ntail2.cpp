#include <Rcpp.h>

using namespace Rcpp;

// samples a column vector from standard multivariate normal
// distribution truncated truncated over the region  [l,u]
// where l > 0
// uses acceptance-rejection from Marsaglia's Rayleigh
// distribution method (1964)

//' @export
// [[Rcpp::export]]
NumericVector ntail2(NumericVector l, NumericVector u) {
  //boundary and vector dimension checks
  if(l.size() != u.size()){
    Rcpp::stop("Non-conformal size in vectors l and u");
  }
  
  int n = l.size();
  
  if(is_true(any(l <= 0))){
    Rcpp::stop("The lower bound l has non-positive entries");
  }
  
  if(is_true(any(l >= u))){
    Rcpp::stop("The lower bound l is larger than vector u for some entries");
  }
  
  //proceed with Marsaglia's Rayleigh distribution method
  NumericVector c(n); 
  NumericVector f(n);
  NumericVector x(n);
  bool rejected = false;
  
  for(int i = 0; i < n; i++){//propose values
    c[i] = std::pow(l[i],2)/2;
    f[i] = std::expm1(c[i] - std::pow(u[i],2)/2);
    x[i] = c[i] - std::log(1.0 + R::runif(0.0,1.0)*f[i]);
    
    if(std::pow(R::runif(0.0,1.0),2)*x[i] > c[i]){//checks if entry "i" is rejected
      rejected = true;//if so we say it is true
      while(rejected == true){//generate new values for this entry until accepted
        x[i] = c[i] - std::log(1.0 + R::runif(0.0,1.0)*f[i]);
        if(std::pow(R::runif(0.0,1.0),2)*x[i] <= c[i]){//checks if accepted
          rejected = false;
        }
      }
    }
    x[i] = std::sqrt(2*x[i]);
  }
  return x;
}



