mvtrnd <- function(n, L, l, u, nu, mu)
{
  # generates the proposals from the exponentially tilted
  # sequential importance sampling pdf;
  # output:    'p', log-likelihood of sample
  #             Z, Gaussian sample
  #             R, random scale parameter so that sqrt(nu)*Z/R is student

  d <- length(l); # Initialization
  eta <- mu[d]; mu[d] <- 0;
  Z <- matrix(0, nrow = d, ncol = n); # create array for variables
  # precompute constants
  const <- log(2*pi) / 2 - lgamma(nu/2) - (nu / 2 - 1) * log(2) + lnNpr(-eta, Inf) + 0.5 * eta^2;
  R <- eta + trandn(rep(-eta, n), rep(Inf, n)); # simulate R~N(eta,1) with R>0
  p <- (nu - 1) * log(R) - eta * R + const; # compute Likelihood Ratio for R
  for(k in 1:d){
    # compute matrix multiplication L*Z
    col <- L[k,1:k] %*% Z[1:k,];
    # compute limits of truncation
    tl <- R * l[k] / sqrt(nu) - mu[k] - col;
    tu <- R * u[k] / sqrt(nu) - mu[k] - col;
    #simulate N(mu,1) conditional on [tl,tu]
    Z[k,] <- mu[k] + trandn(tl, tu);
    # update likelihood ratio
    p <- p + lnNpr(tl,tu)  + 0.5*mu[k]^2 - mu[k] * Z[k,];
  }
  list(p = p, Z = Z, R = R)
}
