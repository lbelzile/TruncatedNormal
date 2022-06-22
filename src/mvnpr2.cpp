#include <RcppArmadillo.h>
using namespace Rcpp;

//'@param n number of Monte Carlo trials
//'@param L Cov(X) = L*L
//'@param l truncation lower bound
//'@param u truncation upper bound
//'@param mu optimized exponential tilting parameter
//'@return P(l < x < u)
//'computes P(l<X<u), where X is normal with
// 'Cov(X)=L*L' and zero mean vector;
// exponential tilting uses parameter 'mu';
// Monte Carlo uses 'n' samples;
//' @export
// [[Rcpp::export]]
List mvnpr2(int n, arma::mat L, NumericVector l, NumericVector u, NumericVector mu) {
  int d = l.size();
  mu[d] = 0;
  int b = 100; // block size (set arbitrarily to 100)
  arma::mat Z(d,b);
  Z.zeros();
  arma::rowvec p(n);
  
  //helper vectors
  arma::colvec lvec;
  arma::colvec uvec;
  arma::colvec muvec;
  arma::colvec shiftvec;
  
  
  arma::colvec col;
  arma::colvec tl;
  arma::colvec tu;
  
  arma::colvec armatrandn;
  arma::colvec armalnNpr;
  Function lnNpr("lnNpr");
  Function trandn("trandn");
  
  NumericVector blocks(std::ceil(n/b));
  for(int i = 0; i < blocks.size(); i++){
    if(n-b >= 0){
      blocks[i] = b;
    } else {
      blocks[i] = n % b;
    }
    n -= b;
  }
  
  
  lvec.ones(blocks[0]);
  uvec.ones(blocks[0]);
  muvec.ones(blocks[0]);
  shiftvec.ones(blocks[0]);
  
  
  for(int j = 0; j < blocks.size();j++){
    
    for(int k = 0; k < d-1; k++){
      
      
      if(k == 0){
        col = trans(L(0,0) * Z.row(0).subvec(0,blocks[j]-1));
      } else {
        col = trans(L.row(k).subvec(0,k) * Z.rows(0,k).cols(0,blocks[j]-1));
      }
      
        
      lvec = l[k]*lvec;
      uvec = u[k]*uvec;
      muvec = mu[k]*muvec;
      
      tl = lvec - muvec - col;
      tu = uvec - muvec - col;
      
      armatrandn = as<arma::colvec>(wrap(trandn(tl,tu)));
      armalnNpr = as<arma::colvec>(wrap(lnNpr(tl,tu)));
        
      Z.row(k) = trans(muvec + armatrandn);
      
      shiftvec = 0.5*std::pow(mu[k],2)*shiftvec;
      p.subvec(b*j,b*j + blocks[j] - 1) = p.subvec(b*j,b*j + blocks[j] - 1) + trans(armalnNpr) + trans(shiftvec) - mu[k]*Z.row(k);
      
      lvec.ones(blocks[j+1]);
      uvec.ones(blocks[j+1]);
      muvec.ones(blocks[j+1]);
      shiftvec.ones(blocks[j+1]);
    }
  
    col = trans(L.row(d-1)*Z.cols(0,blocks[blocks.size()-1]-1));
    
    lvec.ones(100);
    uvec.ones(100);
    
    lvec = l[d-1]*lvec;
    uvec = u[d-1]*uvec;
    
    tl = lvec - col;
    tu = uvec - col;
    
    armalnNpr = as<arma::colvec>(wrap(lnNpr(tl,tu)));
    p.subvec(b*j,b*j + blocks[j] - 1) = p.subvec(b*j,b*j + blocks[j] - 1) + trans(armalnNpr);
  }
  
  for(int i = 0; i < n; i++){
    p(i) = std::exp(p(i));
  }
  
  p = exp(p);
  
  n = b*(blocks.size()-1) + blocks[blocks.size()-1];
  
  double prob = sum(p)/n;
  double sd;
  for(int i = 0; i < n; i++){
    sd += std::pow((p(i)-prob),2);
  }
  sd = std::sqrt(sd/n);
  double relErr = (sd/std::sqrt(n))/prob;
  return Rcpp::List::create(Named("prob") = prob, Named("relErr") = relErr);
}
