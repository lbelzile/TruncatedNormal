tn <-
  function(l,u)
  { # samples a column vector of length=length(l)=length(u)
    # from the standard multivariate normal distribution,
    # truncated over the region [l,u], where -a<l<u<a for some
    # 'a' and l and u are column vectors;
    # uses acceptance rejection and inverse-transform method;
    tol <- 2.05 # controls switch between methods
    # threshold can be tuned for maximum speed for each platform
    # case: abs(u-l)>tol, uses accept-reject from randn
    x <- rep(0, length(l));
    I <- (abs(u-l)>tol)
    if (any(I)){
      x[I] <- trnd(l[I], u[I])
    }
    # case: abs(u-l)< tol, uses inverse-transform
    I <- !I
    if (any(I)){
      pl <- pnorm(l[I])
      pu <- pnorm(u[I]);
      x[I] <- qnorm(pl+(pu-pl)*runif(length(pl)))
    }
    return(x)
  }
