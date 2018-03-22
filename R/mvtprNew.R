mvtpr <- function(n, L, l, u, nu, mu){
  # computes P(l<X<u), where X is student with
  # 'Sig=L*L' and zero mean vector;
  # exponential tilting uses parameter 'mu';
  # Monte Carlo uses 'n' samples;
  d <- length(l); # Initialization
  eta <- mu[d];
  # precompute constants
  const <- log(2*pi) / 2 - lgamma(nu / 2) - (nu / 2 - 1) * log(2) +
    lnNpr(-eta, Inf) + 0.5 * eta^2;
  # create array for variables
  nmax <- ceiling(n * d / 1e7)
  n0 <- ifelse(nmax > 1, nmax, 1) #number of blocks
  nrem <- n #initialize
  sumX <- 0; sumXsq <- 0
  for(iter in 1:n0){
    nst <- min(nrem, ceiling(n/nmax))
    Z <- matrix(0, nrow = d, ncol = nst); 
    R <- eta + trandn(rep(-eta, nst), rep(Inf, nst)); # simulate R~N(eta,1) with R>0
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
    if(n0 > 1){
      nrem <- nrem - nst
      sumX <- sumX + sum(p); sumXsq <- sumXsq + sum(exp(2*log(p)))
    }
  }
  if(n0 == 1){
    meanp <- mean(p)
    sdp <- sdp
  } else{
    meanp <- sumX / n
    sdp <- sqrt((sumXsq - sumX^2/n)/(n-1)) #numerically unstable, but there should not be too many blocks...
  }
  est.prob <- exp(const) * meanp;
  est.relErr <- exp(const) * sdp/(sqrt(n) * est.prob); # relative error
  return(list(prob = est.prob, err = exp(const) * sdp, relErr = est.relErr))
}
