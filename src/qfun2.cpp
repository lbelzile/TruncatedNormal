#include <Rcpp.h>
using namespace Rcpp;

//' Calculate q-function of x
//' Attention: not the same as the tail distribution function
//' @param x Argument vector 
//' @return y Returning vector
//' @references Botev, Z. and L'Écuyer, P. (2016). Simulation from the Normal distribution truncated to an interval in the tail

//' @export
// [[Rcpp::export]]
NumericVector qfun2(NumericVector x) {
  NumericVector y(x.size());
  for(int i = 0; i < x.size(); i++){
    if(std::isinf(x[i]) && x[i] > 0){
      y[i] = 0; //output 0 if entry is infinite
    } else {
      y[i] = std::exp(.5*std::pow(x[i],2)+R::pnorm(x[i],0.0, 1.0, 0,1));
    }
  }
  return y;// [[Rcpp::export]]
}

