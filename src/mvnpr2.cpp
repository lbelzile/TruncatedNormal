#include <Rcpp.h>
using namespace Rcpp;


//

// [[Rcpp::export]]
NumericVector mvnpr2(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


