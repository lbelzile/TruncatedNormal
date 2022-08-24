mvNcdfP <-  function(l, u, SigInv, n = 1e5){
  d=length(l); # basic input check
  if  (length(u) !=d | d != sqrt(length(SigInv)) | any(l > u)){
    stop('l, u, and Sig have to match in dimension with u>l')
  }
  if(d == 1L){
    #warning("Univariate problem not handled; using `pnorm`")
    return(list(prob = exp(lnNpr(a = l * sqrt(SigInv[1]), b = u * sqrt(SigInv[1]))), err = NA, relErr = NA, upbnd = NA))
  }
  # Cholesky decomposition of matrix
  M <- chol(SigInv)
  D <- diag(M)
  if (any(D < 1e-10)){
    warning('Method may fail as covariance matrix is singular!')
  }
  
  u <- u*D;
  l <- l*D; # rescale
  # find optimal tilting parameter via non-linear equation solver
  x0 <- rep(0, 2 * length(l) - 2)
  solvneq <- nleqslv::nleqslv(x0, fn = gradpsiP, jac = jacpsiP,
                              M = M, l = l, u = u, global = "pwldog", method = "Broyden",
                              control = list(maxit = 500L))
  xmu <- solvneq$x
  exitflag <- solvneq$termcd
  flag <- TRUE
  if(!(exitflag %in% 1:2) || !isTRUE(all.equal(solvneq$fvec, rep(0, length(x0)), tolerance = 1e-6))){
    flag <- FALSE
  }
  x = rep(0,d);
  mu = rep(0,d);
  x[2:d] <- xmu[1:(d-1)]
  mu[2:d] <- xmu[d:(2*d-2)] # assign saddlepoint x* and mu*
  
  if(any((t(solve(M)) %*% x - u)[-1] > 0, (-t(solve(M)) %*% x + l)[-1] > 0)){
    warning("Solution to exponential tilting problem using Powell's dogleg method \n  does not lie in convex set l < (T^-T)x < u.")
    flag <- FALSE
  }
  # If Powell dogleg method fails, try constrained convex solver
  if(!flag){
    solvneqc <- alabama::auglag(par = xmu,
                                fn = function(par, l=l, M=M, u=u){
                                  ps <- try(-psyP(x = c(0,par[1:(d-1)]), mu = c(0,par[d:(2*d-2)]),
                                                 l = l, M = M, u = u))
                                  return(ifelse(is.character(ps), -1e10, ps))},
                                gr = function(x, l=l, M=M, u=u){gradpsiP(y=x, M=M,l=l, u=u)},
                                M=M, l=l, u=u,
                                # equality constraints d psi/d mu = 0
                                heq = function(x, l, M, u){gradpsiP(y = x, l=l, M=M, u = u)[d:(2*d-2)]},
                                hin = function(par,...){c((out$u - t(solve(out$L)) %*% c(0,par[1:(d-1)]))[-1] > 0, (t(solve(out$L)) %*% c(0,par[1:(d-1)]) - out$l)[-1])},
                                control.outer = list(trace = FALSE,method="nlminb"))
    if(solvneqc$convergence == 0){
      x <- solvneqc$par[1:(d-1)]
      mu <- solvneqc$par[d:(2*d-2)] 
    } else{
      stop('Did not find a solution to the nonlinear system in `mvNqmc`!') 
    }
  }
  est <- mvnprP(n, M, l, u, mu)
  # compute psi star
  est$upbnd <- exp(psyP(x, M, l, u, mu))
  return(est)
}
