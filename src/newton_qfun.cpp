#include <Rcpp.h>
using namespace Rcpp;


//' Calculate q-function of x
//' Attention: not the same as the tail distribution function
//' @param x Argument vector 
//' @return y Returning vector
//' @references Botev, Z. and L'Écuyer, P. (2016). Simulation from the Normal distribution truncated to an interval in the tail

//' @export
// [[Rcpp::export]]
NumericVector qfun3(NumericVector x) {
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

//' newton iteration method for finding p-quantile of standard
//' multivariate truncated normal distribution
//' newton2 assumes 0<l<u and 0<p<1

//' @export
//' @param p The probability vector for which we are trying to find a value
//' @param l The lower bound vector of the truncated standard normal distribution
//' @param u The upper bound vector of the truncated standard normal distribution
//' @return A bounded vector x for which P(X <= x) = p
//' @references Botev, Z. and L'Écuyer, P. (2016). Simulation from the Normal distribution truncated to an interval in the tail

// [[Rcpp::export]]
NumericVector newton3(NumericVector p, NumericVector l, NumericVector u) {
  //check boundaries and vector dimension
  if(l.size() != u.size()){
    Rcpp::stop("In function 'newton2', boundaries l and u are not the same size");
  }
  
  if(p.size() != l.size()){
    Rcpp::stop("In function 'newton2', vectors p and l or u are not the same size");
  }
  
  int n = l.size();
  NumericVector x(n);
  
  //Function qfun("qfun2");
  NumericVector ql = qfun3(l);
  NumericVector qu = qfun3(u);
  
  NumericVector ind;
  NumericVector val;
  double err = std::numeric_limits<double>::infinity();
  NumericVector del(n);
  
  for(int i = 0; i < n; i++){
    l[i] = std::pow(l[i],2);
    u[i] = std::pow(u[i],2);
    x[i] = std::sqrt(l[i] - 2*std::log(1+p[i]*std::expm1(((l[i]-u[i])/2))));
  } 
  
  NumericVector qx(n);
  double max_del;
  
  while(err > 1.0e-10){
    qx = qfun3(x);
    max_del = -1.0*std::numeric_limits<double>::infinity();
    
    for(int i = 0; i < n; i++){
      del[i] = (-1.0*qx[i])+((1-p[i])*std::exp(0.5*(std::pow(x[i],2)-l[i]))*ql[i])+(p[i]*std::exp(0.5*(std::pow(x[i],2)-u[i]))*qu[i]);
      x[i] = x[i]-del[i]; //Newton's step
      if(max_del < std::abs(del[i])){max_del = std::abs(del[i]);}
    }
    err = max_del;
  }
  
  return x;
}