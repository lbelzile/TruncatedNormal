#' Truncated univariate normal distribution
#' 
#' The function provides efficient state-of-the-art random number generation of a vector of truncated univariate distribution
#' of the same length as the lower bound vector. The function is vectorized and the vector of means \code{mu} and 
#' of standard deviations \code{sd} are recycled.
#'  
#' If \code{mu} or \code{sd} are not specified they assume the default values of 0 and 1, respectively.
#' @name tnorm
#' @param n number of observations
#' @param p vector or matrix of probabilities
#' @param mu vector of means
#' @param sd vector of standard deviations
#' @param lb vector of lower truncation limits
#' @param ub vector of upper truncation limits
#' @param method string, either of \code{fast} or \code{invtransfo}
#' @return vector or matrix of random variates (\code{rtnorm}) or of quantiles (\code{ptnorm}), depending on the input
#' @examples 
#' rtnorm(n = 10, mu = 2, lb = 1:10, ub = 2:11, method = "fast")
#' qtnorm(runif(10), mu = 2, lb = 1:10, ub = 2:11, sd = 1)
NULL

#' Vectorized random number generation from the univariate truncated normal distribution
#' 
#' If \code{n} is 1, the function returns a vector rather than a matrix. If \code{mu} or \code{sd} are not specified they assume the default values of 0 and 1, respectively.
#' @seealso \code{\link{tnorm}}
#' @export
#' @keywords internal
rtnorm <- function(n, mu, sd, lb, ub, method = c("fast","invtransfo")){
  if(any(missing(lb), missing(ub), length(lb) != length(ub))){
    stop("Invalid input in rtnorm")
    }
  method <- match.arg(method)
  d <- length(lb)
  if(missing(mu)){
    mu <- 0 
  }
  if(missing(sd)){
    sd <- 1 
  }
  mu <- rep(mu, length.out = d)
  sd <- rep(sd, length.out = d)
  #To achieve this, first compute X=norminvp(runif(1),(lb-mu)/sig,(ub-mu)/sig) and then set Z=mu+sig*X
  lb <- (lb - mu) / sd
  ub <- (ub - mu) / sd
  if(d == 1L && n > 1){
    res <- mu + sd * switch(method, 
                            fast = trandn(l = rep(lb, n), u = rep(ub, n)),
                            invtransfo = norminvp(p = runif(n), l = rep(lb, n), u = rep(ub, n)))
  } else{
  res <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    res[i,] <- mu + sd * switch(method, 
                      fast = trandn(l = lb, u = ub),
                      invtransfo = norminvp(p = runif(d), l = lb, u = ub))
  }
  if(d == 1L || n == 1L){
    res <- as.vector(res)
  } 
  }
 res 
}

#' Quantile function using the inversion method
#' 
#' If \code{p} is a matrix, the arguments are recycled.  If \code{mu} or \code{sd} are not specified they assume the default values of 0 and 1, respectively.
#' @seealso \code{\link{tnorm}}
#' @export
#' @keywords internal
qtnorm <- function(p, mu, sd, lb, ub){
  if(any(missing(lb), missing(ub), length(lb) != length(ub))){
    stop("Invalid input in rtnorm")
  }
  d <- length(lb)
  if(missing(mu)){
    mu <- 0 
  }
  if(missing(sd)){
    sd <- 1 
  }
  mu <- rep(mu, length.out = d)
  sd <- rep(sd, length.out = d)
  if(d == 1L && length(p) > 1L){
   p <- as.matrix(p, ncol = 1) 
  }
  if(is.matrix(p)){
   stopifnot(ncol(p) == d)
   n <- nrow(p)  
  } else{
   p <- matrix(p, nrow = 1)
   n <- 1 
  }
  #To achieve this, first compute X=norminvp(runif(1),(lb-mu)/sig,(ub-mu)/sig) and then set Z=mu+sig*X
  res <- matrix(0, ncol = d, nrow = n)
  lb <- (lb - mu) / sd
  ub <- (ub - mu) / sd
  for(i in 1:n){
    res[i,] <- mu + sd * norminvp(p = p[i,], l = lb, u = ub)
  }

    if(d == 1L || n == 1L){
      res <- as.vector(res)
    } 
    res 
  }
