#include <Rcpp.h>
using namespace Rcpp;


//samples from standard multivariate normal distribution
//truncated at [l,u]
//uses acceptance-rejection from standard multivariate normal
//distribution (naive method)

//intuition:
//loop through dimensions
//  sample x from standard normal
//    if l < x < u then accept value and go to next entry(dimension)
//    else sample until find acceptable value
//  return x

// benchmark:
//        test replications elapsed relative user.self sys.self
//1  trnd(l, u)          100    0.45     22.5      0.42     0.03
//2 trnd2(l, u)          100    0.02      1.0      0.02     0.00


//'@param l vector truncation lower bound
//'@param u vector truncation upper bound
//'@return sample vector x over truncated region [l,u]
//'@export
// [[Rcpp::export]]
NumericVector trnd2(NumericVector l, NumericVector u) {
  //boundary and vector dimension checks
  if(l.size() != u.size()){
    Rcpp::stop("Non-conformal size in vectors l and u");
  }
  
  int n = l.size();
  
  if(is_true(any(l >= u))){
    Rcpp::stop("The lower bound l is larger than vector u for some entries");
  }
  
  NumericVector x(n);
  
  for(int i = 0; i < n; i++){
    x[i] = R::rnorm(0.0, 1.0);
    //check if value is rejected
    //if so, samples until a new accepted value is found
    while(x[i] < l[i] || u[i] < x[i]){
      x[i] = R::rnorm(0.0,1.0);
    }
  }
  return x;
}
