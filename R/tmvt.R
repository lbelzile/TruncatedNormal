#' Multivariate truncated Student distribution
#' 
#' Density, distribution function and random generation for the multivariate truncated Student distribution
#' with location vector \code{mu}, scale matrix \code{sigma}, lower truncation limit
#' \code{lb}, upper truncation limit \code{ub} and degrees of freedom \code{df}.
#'
#' The truncation limits can include infinite values. The Monte Carlo (\code{type = "mc"}) uses a sample of size \code{B}, while the
#' qausi Monte Carlo (\code{type = "qmc"}) uses a pointset of size \code{ceiling(n/12)} and estimates the relative error using 12 independent randomized QMC estimators. 
#' 
#' \code{pmvt} computes an estimate and a deterministic upper bound of the distribution function of the multivariate normal distribution.
#' Infinite values for vectors \eqn{u} and \eqn{l} are accepted. The Monte Carlo method uses sample size \eqn{n}: the larger \eqn{n}, the smaller the relative error of the estimator. 
#' 
#' @author Leo Belzile, R port from Matlab code by Z. I. Botev
#' @references  Z. I. Botev and P. L'Ecuyer (2015), Efficient probability estimation
#' and simulation of the truncated multivariate Student-t distribution,
#' Proceedings of the 2015 Winter Simulation Conference, pp. 380-391
#' 
#' @section Usage: \preformatted{
#' dtmvt(x, mu, sigma, df, lb, ub, type = c("mc", "qmc"), log = FALSE, B = 1e4)
#' ptmvt(q, mu, sigma, df, lb, ub, type = c("mc", "qmc"), log = FALSE, B = 1e4)
#' rtmvt(n, mu, sigma, df, lb, ub)
#' pmvt(mu, sigma, df, lb = -Inf, ub = Inf, type = c("mc", "qmc"), log = FALSE, B = 1e4)}
#'    
#' @name tmvt
#' @param n number of observations
#' @param x,q vector or matrix of quantiles
#' @param B number of replications for the (quasi)-Monte Carlo scheme
#' @param log logical; if \code{TRUE}, probabilities and density are given on the log scale.
#' @param mu vector of location parameters
#' @param sigma scale matrix
#' @param df degrees of freedom
#' @param lb vector of lower truncation limits
#' @param ub vector of upper truncation limits
#' @param type string, either of \code{mc} or \code{qmc} for Monte Carlo and quasi Monte Carlo, respectively
#' @return \code{dtmvt} gives the density, \code{ptmvt} gives the distribution function, \code{rtmvt} generate random deviates. 
#' @examples
#' d <- 4; lb <- rep(0, d)
#' mu <- runif(d)
#' sigma <- matrix(0.5, d, d) + diag(0.5, d)
#' samp <- rtmvt(n = 10, mu = mu, sigma = sigma, df = 2, lb = lb)
#' loglik <- dtmvt(samp, mu = mu, sigma = sigma, df = 2, lb = lb, log = TRUE)
#' cdf <- ptmvt(samp, mu = mu, sigma = sigma, df = 2, lb = lb, log = TRUE, type = "q")
NULL

#' Density function for the truncated multivariate Student distribution
#' 
#' This function returns the (log)-density of a matrix \code{x} of observations lying in the interval [\code{lb}, \code{ub}].
#' 
#' @seealso \code{\link{tmvt}}
#' @export
#' @inheritParams tmvt
#' @keywords internal
dtmvt <- function(x, mu, sigma, df, lb, ub, type = c("mc", "qmc"), log = FALSE, B = 1e4){
  if (any(missing(x), missing(sigma), missing(df))) {
    stop("Arguments missing in function call to `dtmvnorm`")
  }
  sigma <- as.matrix(sigma)
  if(missing(mu)){
    mu <- rep(0, length.out = ncol(sigma))
  }
  d <- length(mu)
  type <- match.arg(type)
  df <- as.vector(df)[1]
  if(isTRUE(all.equal(df, 0)) || isTRUE(all.equal(df, Inf))){
    return(dtmvnorm(x = x, mu = mu, sigma = sigma, lb = lb, ub = ub, B = B, log = log, type = type)) 
  }
  stopifnot(df > 0, length(mu) == ncol(sigma), nrow(sigma) == ncol(sigma), is.logical(log))
  if (is.vector(x)) {
    stopifnot(length(x) == length(mu))
    x <- matrix(x, nrow = 1, ncol = length(x))
  } else {
    stopifnot(ncol(x) == length(mu))
  }
  if(missing(lb)){
    lb <- rep(-Inf, d) 
  }
  if(missing(ub)){
    ub <- rep(Inf, d) 
  }
  ldens <- as.vector(TruncatedNormal::.dmvt_arma(x, mu = as.vector(mu), df = df, sigma = as.matrix(sigma), logd = TRUE))
  if(!isTRUE(all.equal(as.vector(c(lb, ub)), c(rep(-Inf, d), rep(Inf, d))))){
    kst <- switch(type,
                  mc = mvTcdf(l = lb - mu, u = ub - mu, df = df, Sig = sigma, n = B)$prob,
                  qmc = mvTqmc(l = lb - mu, u = ub - mu, df = df, Sig = sigma, n = B)$prob)
    ldens <- ldens - log(kst)
  }
  for(i in 1:nrow(x)){
    if(any(x[i,] > ub) || any(x[i,] < lb)){
      ldens[i] <- -Inf 
    }
  }
  if(log){
    return(ldens) 
  } else{
    return(exp(ldens))
  }
}

#' Cumulative distribution function of the truncated multivariate Student distribution.
#' 
#' This function returns the (log)-distribution function of a matrix \code{q} of observations lying in the interval [\code{lb}, \code{ub}].
#' 
#' @seealso \code{\link{tmvt}}
#' @export
#' @inheritParams tmvt
#' @keywords internal
ptmvt <- function(q, mu, sigma, df, lb, ub, type = c("mc", "qmc"), log = FALSE, B = 1e4){
  if (any(missing(q), missing(sigma), missing(df))) {
    stop("Arguments missing in function call to `ptmvt`")
  }
  sigma <- as.matrix(sigma)
  if(missing(mu)){
    mu <- rep(0, length.out = ncol(sigma))
  }
  d <- length(mu)
  if(missing(lb)){
    lb <- rep(-Inf, d) 
  }
  if(missing(ub)){
    ub <- rep(Inf, d) 
  }
  type <- match.arg(type)
  df <- as.vector(df)[1]
  if(isTRUE(all.equal(df, 0)) || isTRUE(all.equal(df, Inf))){
    return(ptmvnorm(q = q, mu = mu, sigma = sigma, lb = lb, ub = ub, B = B, log = log, type = type)) 
  }
  stopifnot(df > 0, length(mu) == ncol(sigma), nrow(sigma) == ncol(sigma), is.logical(log))
  if (is.vector(q)) {
    stopifnot(length(q) == length(mu))
    q <- matrix(q, nrow = 1, ncol = length(q))
  } else {
    stopifnot(ncol(q) == length(mu))
  }
  d <- length(mu)
  if(missing(lb)){
    lb <- rep(-Inf, d) 
  }
  if(missing(ub)){
    ub <- rep(Inf, d) 
  }
  prob <- rep(0, nrow(q))
  
  for(i in 1:nrow(q)){
    if(any(q[i,] > ub) || any(q[i,] < lb)){
      prob[i] <- 0
    } else{
      prob[i] <- switch(type,
                        mc = mvTcdf(l = lb - mu, u = q[i,] - mu, df = df, Sig = sigma, n = B)$prob,
                        qmc = mvTqmc(l = lb - mu, u = q[i,] - mu, df = df, Sig = sigma, n = B)$prob)
    }
  }
  kst <- switch(type,
                mc = mvTcdf(l = lb - mu, u = ub - mu, df = df, Sig = sigma, n = B)$prob,
                qmc = mvTqmc(l = lb - mu, u = ub - mu, df = df, Sig = sigma, n = B)$prob)
  if(log){
    return(log(prob) - log(kst)) 
  } else{
    return(prob/kst)
  }
}

#' Random number generator for the truncated multivariate Student distribution.
#' 
#' This function returns a matrix of draws from a multivariate Student distribution truncated on the interval [\code{lb}, \code{ub}].
#' 
#' @seealso \code{\link{tmvt}}
#' @export
#' @inheritParams tmvt
#' @keywords internal
rtmvt <- function(n, mu, sigma, df, lb, ub){
  if (any(missing(sigma), missing(df))) {
    stop("Arguments missing in function call to `rtmvt`")
  }
  sigma <- as.matrix(sigma)
  df <- as.vector(df)[1]
  if(isTRUE(all.equal(df, 0)) || isTRUE(all.equal(df, Inf))){
    return(rtmvnorm(n = n, mu = mu, sigma = sigma, lb = lb, ub = ub)) 
  }
  if(missing(mu)){
    mu <- rep(0, ncol(sigma)) 
  }
  stopifnot(length(mu) == ncol(sigma), nrow(sigma) == ncol(sigma))
  d <- length(mu)
  if(missing(lb)){
    lb <- rep(-Inf, d) 
  }
  if(missing(ub)){
    ub <- rep(Inf, d) 
  }
  stopifnot(length(lb) == length(ub), length(lb) == d, lb <= ub)
  if(!any((ub - lb) < 1e-10)){
    
    if(n == 1){
      as.vector(mvrandt(l = lb, u = ub, Sig = sigma, df = df, n = n, mu = mu))
    } else{
      t(mvrandt(l = lb, u = ub, Sig = sigma, df = df, n = n, mu = mu))
    }
  } else{
    warning("Some variables have a degenerate distribution.")
    ind <- which((ub - lb) >= 1e-10)
    # check covariance matrix
    stopifnot(isSymmetric(sigma), all(eigen(sigma, only.values = TRUE)$value > 0))
    # compute conditional Gaussian
    schurcomp <- function(sigma, ind) {
      stopifnot(c(length(ind) > 0, ncol(sigma) - length(ind) > 0))
      sigma[ind, ind, drop = FALSE] - sigma[ind, -ind, drop = FALSE] %*%
        solve(sigma[-ind, -ind, drop = FALSE]) %*% sigma[-ind, ind, drop = FALSE]
    }  
    sigmap <- c(df + t(lb[-ind] - mu[-ind]) %*% solve(sigma[-ind, -ind, drop = FALSE]) %*% (lb[-ind] - mu[-ind])) /
      (df + d - length(ind)) * schurcomp(sigma, ind)
    mup <- c(mu[ind] + sigma[ind, -ind, drop = FALSE] %*% solve(sigma[-ind, -ind, drop = FALSE]) %*% (lb[-ind] - mu[-ind]))
    # matrix to store results
    res <- matrix(0, nrow = n, ncol = d)
    res[, -ind] <- rep(lb[-ind], each = n)
    if(n == 1){
      res[, ind] <-  as.vector(mvrandt(l = lb[ind], u = ub[ind], Sig = sigmap, 
                                       df = df + d - length(ind), n = n, mu = mup))
      res <- as.vector(res)
    } else{
      res[, ind] <- t(mvrandt(l = lb[ind], u = ub[ind], Sig = sigmap, df = df + d - length(ind), n = n, mu = mup))
    }
  return(res)
  }
}

#' Distribution function of the multivariate Student distribution for arbitrary limits
#' 
#' This function computes the distribution function of a multivariate normal distribution vector for an arbitrary rectangular region [\code{lb}, \code{ub}].
#' \code{pmvt} computes an estimate and the value is returned along with a relative error and a deterministic upper bound of the distribution function of the multivariate normal distribution.
#' Infinite values for vectors \eqn{u} and \eqn{l} are accepted. The Monte Carlo method uses sample size \eqn{n}: the larger the sample size, the smaller the relative error of the estimator. 
#' @references  Z. I. Botev and P. L'Ecuyer (2015), Efficient probability estimation
#' and simulation of the truncated multivariate Student-t distribution,
#' Proceedings of the 2015 Winter Simulation Conference, pp. 380-391
#' @export
#' @author \code{Matlab} code by Zdravko I. Botev, \code{R} port by Leo Belzile
#' @inheritParams tmvt
#' @export
#' @examples
#' d <- 15; nu <- 30;
#' l <- rep(2, d); u <- rep(Inf, d);
#' sigma <- 0.5 * matrix(1, d, d) + 0.5 * diag(1, d);
#' est <- pmvt(lb = l, ub = u, sigma = sigma, df = nu)
#' # mvtnorm::pmvt(lower = l, upper = u, df = nu, sigma = sigma)
#' \dontrun{
#' d <- 5
#' sigma <- solve(0.5 * diag(d) + matrix(0.5, d, d))
#' # mvtnorm::pmvt(lower = rep(-1,d), upper = rep(Inf, d), df = 10, sigma = sigma)[1]
#' pmvt(lb = rep(-1, d), ub = rep(Inf, d), sigma = sigma, df = 10)
#' }
pmvt  <- function(mu, sigma, df, lb = -Inf, ub = Inf, type = c("mc", "qmc"), log = FALSE, B = 1e4){
  type <- match.arg(type)
  if(any(missing(df), missing(sigma))){
    stop("Missing arguments in `pmvt`")
  }
  sigma <- as.matrix(sigma)
  lb <- rep(lb, length.out = ncol(sigma))
  ub <- rep(ub, length.out = ncol(sigma))
  stopifnot(length(lb) == length(ub))
  if(!missing(mu)){
    stopifnot(length(mu) == length(lb))
    lb <- lb - mu
    ub <- ub - mu
  }
  integ <- switch(type,
                  mc = mvTcdf(l = lb, u = ub, Sig = sigma, df = df, n = B),
                  qmc = mvTqmc(l = lb, u = ub, Sig = sigma, df = df, n = B))
  res <- integ$prob
  attributes(res)$relerr <- integ$relErr
  attributes(res)$upbnd <- integ$upbnd
  res
}
