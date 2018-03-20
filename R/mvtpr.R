mvtpr <- function(n, L, l, u, nu, mu){
  # computes P(l<X<u), where X is student with
  # 'Sig=L*L' and zero mean vector;
  # exponential tilting uses parameter 'mu';
  # Monte Carlo uses 'n' samples;
  d <- length(l); # Initialization
  eta <- mu[d];
  Z <- matrix(0, nrow = d, ncol = n); # create array for variables
  # precompute constants
  const <- log(2*pi) / 2 - lgamma(nu / 2) - (nu / 2 - 1) * log(2) +
    lnNpr(-eta, Inf) + 0.5 * eta^2;
  R <- eta + trandn(rep(-eta, n), rep(Inf, n)); # simulate R~N(eta,1) with R>0
  p <- (nu - 1) * log(R) - eta * R; # compute Likelihood Ratio for R
  R <- R / sqrt(nu); # scale parameter divided by nu
  for(k in 1:(d-1)){
    # compute matrix multiplication L*Z
    col <- L[k,1:k] %*% Z[1:k,];
    # compute limits of truncation
    tl <- R * l[k] - mu[k] - col;
    tu <- R * u[k] - mu[k] - col;
    #simulate N(mu,1) conditional on [tl,tu]
    Z[k,] <- mu[k] + trandn(tl, tu);
    # update likelihood ratio
    p <- p + lnNpr(tl, tu) + 0.5 * mu[k]^2 - mu[k] * Z[k,];
  }
  # deal with final Z(d) which need not be simulated
  col <- c(L[d,] %*% Z);
  tl <- R * l[d] - col
  tu <- R * u[d] - col;
  p <- p + lnNpr(tl, tu); # update LR corresponding to Z(d)
  # now switch back from logarithmic scale
  p <- exp(p)
  est.prob <- exp(const) * mean(p);
  est.relErr <- exp(const) * sd(p)/(sqrt(n) * est.prob); # relative error
  return(list(prob = est.prob, err = exp(const) * sd(p), relErr = est.relErr))
}
