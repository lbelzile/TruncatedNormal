#' Normal quantile function (high precision)
#'
#' Computes with tail-precision the quantile function
#'  of the standard normal distribution at \eqn{0\le p\le 1},
#'  and truncated to the interval \eqn{[l,u]}.
#' Infinite values for vectors \eqn{l} and \eqn{u} are accepted.
#' @param p quantile  at \eqn{0\le p\le 1}
#' @param l lower truncation limit
#' @param u upper truncation limit
#' @keywords internal
#' @author Zdravko I. Botev
#' @references Z. I. Botev (2017), \emph{The Normal Law Under Linear Restrictions:
#' Simulation and Estimation via Minimax Tilting}, Journal of the Royal
#' Statistical Society, Series B, \bold{79} (1), pp. 1--24.
#'
#' @details
#'  Suppose we wish to simulate a random variable \eqn{Z} drawn from \eqn{N(\mu,\sigma^2)} and
#'  conditional on \eqn{l<Z<u} using the inverse transform method.
#'  To achieve this, first compute
#'  \code{X=norminvp(runif(1),(l-mu)/sig,(u-mu)/sig)} and then set
#'  \code{Z=mu+sig*X}
#' @return quantile value of the truncated normal distribution.
#' @note  If you wish to simulate truncated normal variables fast, use \code{\link{trandn}}.
#'  Using \code{norminvp}  is advisable only when needed, for example,
#'  in quasi-Monte Carlo or antithetic sampling, where the inverse transform method
#'  is unavoidable.
#' @export
#' @seealso \code{\link{trandn}}
#' @examples
#'  d <- 150 # simulate via inverse transform method
#'  norminvp(runif(d),l = 1:d, u = rep(Inf, d))
norminvp <- function (p, l, u){
  if ((length(l) != length(p)) | any(l > u) | any(p > 1) |  any(p < 0)) {
    stop("l, u, and p must be the same length with u > l and 0 <= p <= 1")
  }
  x = rep(NaN, length(l))
  K <- ((p < 1) + (p > 0)) == 2L
  if(any(K)){
    x[K] = cases(p = p[K], l = l[K], u = u[K])
  }
  if(sum(K) == length(p)){
   return(x) 
  } 
  I1 <- (p == 1)
  x[I1] = u[I1]
  I2 <- (p == 0)
  x[I2] = l[I2]
  return(x)
}