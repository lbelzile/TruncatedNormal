#' Truncated multivariate normal cumulative distribution  (quasi-Monte Carlo)
#'
#'
#' Computes an estimate and a deterministic upper bound of the probability Pr\eqn{(l<X<u)},
#' where \eqn{X} is a zero-mean multivariate normal vector
#' with covariance matrix \eqn{\Sigma}, that is, \eqn{X} is drawn from \eqn{N(0,\Sigma)}.
#' Infinite values for vectors \eqn{u} and \eqn{l} are accepted.
#' The Monte Carlo method uses sample size \eqn{n}: 
#' the larger \eqn{n}, the smaller the relative error of the estimator.
#'
#' @inheritParams mvNcdf
#' @details Suppose you wish to estimate Pr\eqn{(l<AX<u)},
#'  where \eqn{A} is a full rank matrix
#'  and \eqn{X} is drawn from \eqn{N(\mu,\Sigma)}, then you simply compute
#'  Pr\eqn{(l-A\mu<AY<u-A\mu)},
#'  where \eqn{Y} is drawn from \eqn{N(0, A\Sigma A^\top)}.
#' @return  a list with components
#' \itemize{
#' \item{\code{prob}: }{estimated value of probability Pr\eqn{(l<X<u)}}
#' \item{\code{relErr}: }{estimated relative error of estimator}
#' \item{\code{upbnd}: }{ theoretical upper bound on true Pr\eqn{(l<X<u)}}
#' }
#' @author Zdravko I. Botev
#' @references Z. I. Botev (2017), \emph{The Normal Law Under Linear Restrictions:
#' Simulation and Estimation via Minimax Tilting}, Journal of the Royal
#' Statistical Society, Series B, \bold{79} (1), pp. 1--24.
#'
#' @note This version uses a Quasi Monte Carlo (QMC) pointset
#' of size \code{ceiling(n/12)} and estimates the relative error
#' using 12 independent randomized QMC estimators. QMC
#' is slower than ordinary Monte Carlo,
#' but is also likely to be more accurate when \eqn{d<50}.
#' For high dimensions, say \eqn{d>50}, you may obtain the same accuracy using
#' the (typically faster) \code{\link{mvNcdf}}.
#'
#' @seealso \code{\link{mvNcdf}}, \code{\link{mvrandn}}
#' @export
#' @keywords internal
#' @examples
#' d <- 15 
#' l <- 1:d
#' u <- rep(Inf, d)
#' Sig <- matrix(rnorm(d^2), d, d)*2 
#' Sig <- Sig %*% t(Sig)
#' mvNqmc(l, u, Sig, 1e4) # compute the probability
mvNqmc <- function(l, u, Sig, n = 1e5){
  d <- length(l) # basic input check
   if  (length(u) != d | d != sqrt(length(Sig)) | any(l > u)){
      stop('l, u, and Sig have to match in dimension with u>l')
    }
    if(d == 1L){
      #warning("Univariate problem not handled; using `pnorm`")
      return(list(prob = exp(lnNpr(a = l / sqrt(Sig[1]), b = u / sqrt(Sig[1]))), err = NA, relErr = NA, upbnd = NA))
    }
    # Cholesky decomposition of matrix
    out <- cholperm(Sig, l, u) 
    L <- out$L 
    l <- out$l 
    u <- out$u 
    D <- diag(L)
    if (any(D < 1e-10)){
      warning('Method may fail as covariance matrix is singular!')
    }
    L <- L / D
    u <- u / D
    l <- l / D # rescale
    diag(L) <- rep(0, d) # remove diagonal
  # find optimal tilting parameter via non-linear equation solver
  x0 <- rep(0, 2 * length(l) - 2)
  solvneq <- nleqslv::nleqslv(x0, 
                              fn = gradpsi, 
                              jac = jacpsi,
                              L = L, l = l, u = u, 
                              global = "pwldog", 
                              method = "Newton",
                              control = list(maxit = 500L))
  xmu <- solvneq$x
  exitflag <- solvneq$termcd
  flag <- TRUE
  if(!(exitflag %in% c(1,2)) || !isTRUE(all.equal(solvneq$fvec, rep(0, length(x0)), tolerance = 1e-6))){
    flag <- FALSE
  }
  # assign saddlepoint x* and mu*
  x <- xmu[1:(d-1)]
  mu <- xmu[d:(2*d-2)] 
  # check the constraints
  if(any((out$L %*% c(x,0) - out$u)[-d] > 0, (-out$L %*% c(x,0) + out$l)[-d] > 0)){
   warning("Solution to exponential tilting problem using Powell's dogleg method \n  does not lie in convex set l < Lx < u.")
   flag <- FALSE
  }
  # If Powell dogleg method fails, try constrained convex solver
  if(!flag){
    solvneqc <- alabama::auglag(par = xmu,
                    fn = function(par, l=l, L=L, u=u){
                      ps <- try(-psy(x = c(par[1:(d-1)],0), mu = c(par[d:(2*d-2)],0),
                                     l = l, L = L, u = u))
                      return(ifelse(is.character(ps), -1e10, ps))},
                    gr = function(x, l=l, L=L, u=u){gradpsi(y=x, L=L,l=l, u=u)},
                    L=L, l=l, u=u,
                    # equality constraints d psi/d mu = 0
                  heq = function(x, l, L, u){gradpsi(y = x, l=l, L=L, u = u)[d:(2*d-2)]},
                  hin= function(par,...){c((out$u - out$L %*% c(par[1:(d-1)],0))[-d] > 0, (out$L %*% c(par[1:(d-1)],0) - out$l)[-d])},
                  control.outer = list(trace = FALSE,method="nlminb"))
    if(solvneqc$convergence == 0){
      x <- solvneqc$par[1:(d-1)]
      mu <- solvneqc$par[d:(2*d-2)] 
    } else{
      stop('Did not find a solution to the nonlinear system in `mvNqmc`!') 
    }
  }
  p <- rep(0, 12)
  for (i in 1:12){ # repeat randomized QMC
    p[i] <- mvnprqmc(ceiling(n/12), L = L, l = l, u = u, mu = mu)
  }
  prob <- mean(p) # average of QMC estimates
  relErr <- sd(p)/(sqrt(12) * prob) # relative error
  upbnd <- exp(psy(x, L, l, u, mu)) # compute psi star
  est <- list(prob = prob, relErr = relErr, upbnd = upbnd)
  return(est)
}
