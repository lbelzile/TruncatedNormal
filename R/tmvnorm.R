# Wrappers around the functions in the TruncatedNormal package, following the R convention
# dtmvnorm, dtmvt,...

#' Multivariate truncated normal distribution
#' 
#' Density, distribution function and random generation for the multivariate truncated normal distribution
#' with mean vector \code{mu}, covariance matrix \code{sigma}, lower truncation limit \code{lb} and upper truncation limit \code{ub}. 
#' The truncation limits can include infinite values. The Monte Carlo (\code{type = "mc"}) uses a sample of size \code{B}, while the
#' quasi Monte Carlo (\code{type = "qmc"}) uses a pointset of size \code{ceiling(n/12)} and estimates the relative error using 12 independent randomized QMC estimators. 
#' 
#' @author Zdravko I. Botev, Leo Belzile (wrappers)
#' @references Z. I. Botev (2017), \emph{The normal law under linear restrictions:
#' simulation and estimation via minimax tilting}, Journal of the Royal
#' Statistical Society, Series B, \bold{79} (1), pp. 1--24.
#' 
#' 
#' @section Usage: \preformatted{
#' dtmvnorm(x, mu, sigma, lb, ub, type = c("mc", "qmc"), log = FALSE, B = 1e4)
#' ptmvnorm(q, mu, sigma, lb, ub, type = c("mc", "qmc"), log = FALSE, B = 1e4)
#' rtmvnorm(n, mu, sigma, lb, ub)}
#' 
#' @name tmvnorm
#' @param n number of observations
#' @param x,q vector of quantiles
#' @param B number of replications for the (quasi)-Monte Carlo scheme
#' @param log logical; if \code{TRUE}, probabilites and density are given on the log scale.
#' @param mu vector of location parameters
#' @param sigma covariance matrix
#' @param lb vector of lower truncation limits
#' @param ub vector of upper truncation limits
#' @param type string, either of \code{mc} or \code{qmc} for Monte Carlo and quasi Monte Carlo, respectively
#' @return \code{dtmvnorm} gives the density, \code{ptmvnorm} and \code{pmvnorm} give the distribution function of respectively the truncated and multivariate Gaussian distribution and \code{rtmvnorm} generate random deviates. 
#' @examples 
#' d <- 4; lb <- rep(0, d)
#' mu <- runif(d)
#' sigma <- matrix(0.5, d, d) + diag(0.5, d)
#' samp <- rtmvnorm(n = 10, mu = mu, sigma = sigma, lb = lb)
#' loglik <- dtmvnorm(samp, mu = mu, sigma = sigma, lb = lb, log = TRUE)
#' cdf <- ptmvnorm(samp, mu = mu, sigma = sigma, lb = lb, log = TRUE, type = "q")

NULL

#' Density function for the truncated multivariate normal distribution
#' 
#' This function returns the (log)-density of a matrix \code{x} of observations lying in the interval [\code{lb}, \code{ub}].
#' 
#' @seealso \code{\link{tmvnorm}}
#' @export
#' @keywords internal
dtmvnorm <- function(x, mu, sigma, lb, ub, log = FALSE, type = c("mc", "qmc"), B = 1e4){
  if (any(missing(x), missing(mu), missing(sigma))) {
    stop("Arguments missing in function call to `dtmvnorm`")
  }
  sigma <- as.matrix(sigma)
  type <- match.arg(type)
  if(missing(mu)){
    mu <- rep(0, length.out = ncol(sigma))
  }
  d <- length(mu)
  stopifnot(d == ncol(sigma), d == ncol(sigma), is.logical(log))
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
  ldens <- as.vector(.dmvnorm_arma(x, mu = as.vector(mu), sigma = as.matrix(sigma), logd = TRUE))
  if(!isTRUE(all.equal(as.vector(c(lb, ub)), c(rep(-Inf, d), rep(Inf, d))))){
    kst <- switch(type,
                  mc = mvNcdf(l = lb - mu, u = ub - mu, Sig = sigma, n = B)$prob,
                  qmc = mvNqmc(l = lb - mu, u = ub - mu, Sig = sigma, n = B)$prob)
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

#' Cumulative distribution function of the truncated multivariate normal distribution.
#' 
#' This function returns the (log)-distribution function of a matrix \code{q} of observations lying in the interval [\code{lb}, \code{ub}].
#' 
#' @seealso \code{\link{tmvnorm}}
#' @export
#' @keywords internal
ptmvnorm <- function(q, mu, sigma, lb, ub, log = FALSE, type = c("mc", "qmc"), B = 1e4){
  if (any(missing(q), missing(sigma))) {
    stop("Arguments missing in function call to `ptmvnorm`")
  }
  sigma <- as.matrix(sigma)
  type <- match.arg(type)
  if(missing(mu)){
    mu <- rep(0, length.out = ncol(sigma))
  }
  d <- length(mu)
  stopifnot(length(mu) == ncol(sigma), nrow(sigma) == ncol(sigma), is.logical(log))
  if (is.vector(q)) {
    stopifnot(length(q) == length(mu))
    q <- matrix(q, nrow = 1, ncol = length(q))
  } else {
    stopifnot(ncol(q) == length(mu))
  }
  if(missing(lb)){
    lb <- rep(-Inf, d) 
  }
  if(missing(ub)){
    ub <- rep(Inf, d) 
  }
  prob <- rep(0, nrow(q))
  
  for(i in 1:nrow(q)){
    if(any(q[i,] > ub) || any(q[i,] < lb)){
      prob[i] <- NA
    } else{
      prob[i] <- switch(type,
                        mc = mvNcdf(l = lb - mu, u = q[i,] - mu, Sig = sigma, n = B)$prob,
                        qmc = mvNqmc(l = lb - mu, u = q[i,] - mu, Sig = sigma, n = B)$prob)
    }
  }
  kst <- switch(type,
                mc = mvNcdf(l = lb - mu, u = ub - mu, Sig = sigma, n = B)$prob,
                qmc = mvNqmc(l = lb - mu, u = ub - mu, Sig = sigma, n = B)$prob)
  if(log){
    return(log(prob) - log(kst)) 
  } else{
    return(prob/kst)
  }
}

#' Random number generator for the truncated multivariate normal distribution.
#' 
#' This function returns a matrix of draws from a multivariate normal distribution truncated on the interval [\code{lb}, \code{ub}].
#' 
#' @seealso \code{\link{tmvnorm}}
#' @export
#' @keywords internal
rtmvnorm <- function(n, mu, sigma, lb, ub){
  if (missing(sigma)) {
    stop("Arguments missing in function call to `rtmvnorm`")
  }
  sigma <- as.matrix(sigma)
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
  if(n == 1){
    as.vector(mvrandn(l = lb, u = ub, Sig = sigma, n = n, mu = mu)) + mu
  } else{
    t(mvrandn(l = lb, u = ub, Sig = sigma, n = n, mu = mu) + mu)
  }
}

#' Distribution function of the multivariate normal distribution for arbitrary limits
#' 
#' This function computes the distribution function of a multivariate normal distribution vector for an arbitrary rectangular region [\code{lb}, \code{ub}].
#' \code{pmvnorm} computes an estimate and the value is returned along with a relative error and a deterministic upper bound of the distribution function of the multivariate normal distribution.
#' Infinite values for vectors \eqn{u} and \eqn{l} are accepted. The Monte Carlo method uses sample size \eqn{n}: the larger the sample size, the smaller the relative error of the estimator. 
#' @author Zdravko I. Botev, Leo Belzile (wrappers)
#' @references Z. I. Botev (2017), \emph{The normal law under linear restrictions:
#' simulation and estimation via minimax tilting}, Journal of the Royal
#' Statistical Society, Series B, \bold{79} (1), pp. 1--24.
#' @inheritParams tmvnorm
#' @seealso \code{\link[mvtnorm]{pmvnorm}}
#' @export
#' @examples 
#' #From mvtnorm
#' mean <- rep(0, 5)
#' lower <- rep(-1, 5)
#' upper <- rep(3, 5)
#' corr <- matrix(0.5, 5, 5) + diag(0.5, 5)
#' prob <- pmvnorm(lb = lower, ub = upper, mu = mean, sigma = corr)
#' stopifnot(pmvnorm(lb = -Inf, ub = 3, mu = 0, sigma = 1) == pnorm(3))
pmvnorm <- function(mu, sigma, lb = -Inf, ub = Inf, B = 1e4, type = c("mc", "qmc"), log = FALSE){
  type <- match.arg(type)
  sigma <- as.matrix(sigma)
  if(missing(sigma)){
    stop("Missing arguments in `ncmvnorm`")
  }
  lb <- rep(lb, length.out = ncol(sigma))
  ub <- rep(ub, length.out = ncol(sigma))
  if(!missing(mu)){
    stopifnot(length(mu) == length(lb))
    lb <- lb - mu
    ub <- ub - mu
  }
  if(missing(sigma)){
    sigma <- diag(1, length(lb)) 
  }
  integ <- switch(type,
                  mc = mvNcdf(l = lb, u = ub, Sig = sigma, n = B),
                  qmc = mvNqmc(l = lb, u = ub, Sig = sigma, n = B))
  res <- integ$prob
  attributes(res)$relerr <- integ$relErr
  attributes(res)$upbnd <- integ$upbnd
  res
}
