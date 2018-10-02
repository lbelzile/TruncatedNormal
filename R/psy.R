psy <-  function(x, L, l, u, mu)
  {# implements psi(x, mu) assume scaled 'L' without diagonal
    d <- length(u) 
    x[d] <- 0 
    mu[d] <- 0
    # compute now ~l and ~u
    cv <- L %*% x 
    sum(lnNpr(l - mu - cv, u - mu - cv) + 0.5 * mu^2 - x * mu)
  }

