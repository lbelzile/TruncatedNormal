#include <Rcpp.h>
using namespace Rcpp;

//' Calculate q-function of x (see Botev, L'Écuyer (2016))
//' Attention: not the same as the tail distribution function
//' 
//' @export
// [[Rcpp::export]]
NumericVector qfun2(NumericVector x) {
  NumericVector y(x.size());
  for(int i = 0; i < x.size(); i++){
    if(std::isinf(x[i])){y[i] = 0; std::cout<<"yes"<<std::endl;}//output 0 if entry is infinite
    y[i] = std::exp(.5*std::pow(x[i],2)+R::pnorm(x[i],0.0, 1.0, 0,1));
  }
  return y;// [[Rcpp::export]]
  ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
