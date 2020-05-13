mvtprqmc <- function(n, L, l, u, nu, mu){
  # computes P(l<X<u), where X is student with
  # Sig=L*L', zero mean vector, and degrees of freedom 'nu';
  # exponential tilting uses parameter 'mu';
  # Quasi Monte Carlo uses 'n' samples;
  d <- length(l); # Initialization
  eta <- mu[d];
  Z <- matrix(0, nrow = d, ncol = n); # create array for variables
  # QMC pointset
  if(n*(d-1) > 2e7){
    warning("High memory requirements for storage of QMC sequence\nConsider reducing n")
  }
  x <- as.matrix(randtoolbox::sobol(n, dim = d - 1, init = TRUE, scrambling = 1, seed = ceiling(1e6 * runif(1))))
  # x <- as.matrix(qrng::sobol(n = n, d = d - 1, randomize = TRUE))
  #Fixed 21.03.2018 to ensure that if d=2, no error returned
  # Monte Carlo uses 'n' samples;
  # precompute constants
  const <- log(2*pi) / 2 - lgamma(nu / 2) - (nu / 2 - 1) * log(2) +
    lnNpr(-eta, Inf) + 0.5 * eta^2;
  R <- eta + trandn(rep(-eta, n), rep(Inf, n));
  # simulate R~N(eta,1) with R>0
  p <- (nu - 1) * log(R) - eta * R; # compute Likelihood Ratio for R
  R <- R / sqrt(nu); # scale parameter divided by nu
  for(k in 1:(d-1)){
    # compute matrix multiplication L*Z
    col <- c(L[k,1:k] %*% Z[1:k,])
    #bottleneck, but hard to reduce
    # compute limits of truncation
    tl <- R * l[k] - mu[k] - col;
    tu <- R * u[k] - mu[k] - col;
    #simulate N(mu,1) conditional on [tl,tu]
    Z[k,] <- mu[k] + norminvp(x[1:n,k], tl, tu);
    # update likelihood ratio
    p <- p + lnNpr(tl,tu)  + .5*mu[k]^2-mu[k]*Z[k,];
  }
  # deal with final Z(d) which need not be simulated
  col <- c(L[d,] %*% Z);
  tl <- R * l[d] - col
  tu <- R * u[d] - col;
  p <- p + lnNpr(tl, tu); # update LR corresponding to Z(d)
  return(exp(const)*mean(exp(p)))
}
