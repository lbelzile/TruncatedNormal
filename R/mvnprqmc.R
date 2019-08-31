mvnprqmc <- function(n, L, l, u, mu){
    # computes P(l<X<u), where X is normal with
    # Cov(X)=L*L' and zero mean vector;
    # exponential tilting uses parameter 'mu';
    # Quasi Monte Carlo uses 'n' samples;
    d <- length(l) # Initialization
    Z <- matrix(0, d, n) # create array for variables
    # QMC pointset - 1 is Owen's scrambling
    if(n*(d-1) > 2e7){
      warning("High memory requirements for storage of QMC sequence\nConsider reducing n")
    }
    x <- as.matrix(randtoolbox::sobol(n, dim = d-1, init =TRUE, scrambling = 1, seed=ceiling(1e6*runif(1))))
    ## Similar option in fOptions package. Problem: sobol sequence can overflow (values above 1).
    # x <- qrng::sobol(n = n, d = d - 1, randomize = TRUE)
    p <- 0
    for (k in 1:(d-1)){
      # compute matrix multiplication L*Z
      if(k > 1){
        col <- crossprod(L[k, 1:(k-1)], Z[1:(k-1),]) 
      } else{
        col <- rep(0, n)
      }
      # compute limits of truncation
      tl <- l[k] - mu[k] - col
      tu <- u[k] - mu[k] - col
      # simulate N(mu,1) conditional on [tl, tu] via QMC
      if (d > 2){
        Z[k,] <- mu[k] + norminvp(x[, k], tl, tu)
      } else {
        Z[k,] <- mu[k] + norminvp(x, tl, tu)
      }
      # update likelihood ratio
      p <- p + lnNpr(tl, tu) + 0.5 * mu[k]^2 - mu[k] * Z[k,]
    }
    # deal with final Z(d) which need not be simulated
    col <- L[d,] %*% Z
    tl <- l[d] - col
    tu <- u[d] - col
    p <- p + lnNpr(tl, tu) # update LR corresponding to Z(d)
    mean(exp(p)) # now switch back from logarithmic scale
  }
