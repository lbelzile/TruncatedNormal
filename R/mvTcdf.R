#' Truncated multivariate student cumulative distribution
#'
#' Computes an estimator of the probability Pr\eqn{(l<X<u)},
#' where \eqn{X} is a centered multivariate student vector
#' with scale matrix \code{Sig} and degrees of freedom \code{df}.
#' Infinite values for vectors \code{u} and \code{l} are accepted.
#'
#' @details Monte Carlo method uses sample size \code{n}; the larger
#' the \code{n}, the smaller the relative error of the estimator;
#'
#' @inheritParams mvrandt
#' @return a list with components
#' \itemize{
#' \item{\code{prob}: }{estimated value of probability Pr\eqn{(l<X<u)} }
#' \item{\code{relErr}: }{estimated relative error of estimator}
#' \item{\code{upbnd}: }{theoretical upper bound on true Pr\eqn{(l<X<u)} }
#'}
#' @note If you want to estimate Pr\eqn{(l<Y<u)},
#' where \eqn{Y} follows a Student distribution with \code{df} degrees of freedom,
#' location vector \code{m} and scale matrix \code{Sig},
#' then use \code{mvTqmc(Sig, l - m, u - m, nu, n)}.
#'
#' @examples
#'  d <- 15; nu <- 30;
#'  l <- rep(2, d); u <- rep(Inf, d);
#'  Sig <- 0.5 * matrix(1, d, d) + 0.5 * diag(1, d);
#'  est <- mvTcdf(l, u, Sig, nu, n = 1e4)
#'  # mvtnorm::pmvt(lower = l, upper = u, df = nu, sigma = Sig)
#' \dontrun{
#' d <- 5
#' Sig <- solve(0.5*diag(d)+matrix(0.5, d,d))
#' # mvtnorm::pmvt(lower = rep(-1,d), upper = rep(Inf, d), df = 10, sigma = Sig)[1]
#' mvTcdf(rep(-1, d), u = rep(Inf, d), Sig = Sig, df = 10, n=1e4)$prob
#' }
#' @seealso \code{\link{mvTqmc}}, \code{\link{mvrandt}}, \code{\link{mvNqmc}}, \code{\link{mvrandn}}
#' @references Z. I. Botev (2017), \emph{The Normal Law Under Linear Restrictions:
#' Simulation and Estimation via Minimax Tilting}, Journal of the Royal
#' Statistical Society, Series B, \bold{79} (1), pp. 1--24
#' @references  Z. I. Botev and P. L'Ecuyer (2015), Efficient probability estimation
#' and simulation of the truncated multivariate Student-t distribution,
#' Proceedings of the 2015 Winter Simulation Conference, pp. 380-391
#' @export
#' @keywords internal
#' @author \code{Matlab} code by Zdravko Botev, \code{R} port by Leo Belzile
#' @importFrom randtoolbox sobol
mvTcdf <- function(l, u, Sig, df, n = 1e5){
  d <- length(l)
  if (length(u) != d | d != sqrt(length(Sig)) | any(l > u)) {
    stop("The dimensions of l, u, and Sig have to match and u > l")
  }
  if(d == 1L){
    #warning("Univariate problem not handled; using `pt`.")
    return(list(prob = pt(q = u/sqrt(Sig[1]), df = df) - pt(q = l/sqrt(Sig[1]), df = df), err = NA, relErr = NA, upbnd = NA))
  }
  out <- cholperm(Sig, l, u)
  Lfull <- out$L
  D <- diag(Lfull)
  if (any(D < 1e-10)) {
    warning("Method may fail as covariance matrix is singular!")
  }
  L <- Lfull/D - diag(d)
  u <- out$u/D
  l <- out$l/D
  #Starting value
  x0 <- rep(0, 2*d); x0[2*d] <- sqrt(df); x0[d] <- log(x0[2*d])
  solvneq <- nleqslv::nleqslv(x = x0, fn = gradpsiT, L = L, l = l, u = u, nu = df,
                              global = "pwldog", method = "Broyden")
  soln <- solvneq$x
  #fval <- solvneq$fvec
  exitflag <- solvneq$termcd
  if(!(exitflag %in% 1:2) || !all.equal(solvneq$fvec, rep(0, length(x0)))){
    warning('Method may fail as covariance matrix is close to singular!')
  }
  # assign saddlepoint x* and mu*
  soln[d] <- exp(soln[d])
  x <- soln[1:d]
  mu <- soln[(d+1):length(soln)]
  est <- mvtpr(n = n, L = L, l = l, u = u, nu = df, mu = mu)
  # compute psi star
  est$upbnd <- psyT(x = x, L = L, l = l, u = u, nu = df, mu = mu) # compute psi star
  if(est$upbnd < -743){
    warning('Natural log of probability is less than -743, yielding 0 after exponentiation!')
  }
  est$upbnd <- exp(est$upbnd)
  return(est)
}
