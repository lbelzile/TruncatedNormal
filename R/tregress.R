#' Truncated student generator for Bayesian regression simulation
#'
#' Simulates \code{n} random vectors \eqn{X} exactly distributed
#' from the \code{d}-dimensional Student distribution with
#' \code{df=}\eqn{\nu} degrees of freedom, mean zero and scale matrix
#' \code{sigma}, conditional on \eqn{l<X<u},
#'
#' @inheritParams tmvt
#'
#' @author \code{Matlab} code by Zdravko Botev, \code{R} port by Leo Belzile
#' @export
#' @references  Z. I. Botev and P. L'Ecuyer (2015), Efficient probability estimation
#' and simulation of the truncated multivariate Student-t distribution,
#' Proceedings of the 2015 Winter Simulation Conference, pp. 380-391,
#' @return list with components
#' \itemize{
#' \item{\code{R}: } \code{n} vector of scale
#' \item{\code{Z}: } a \code{d} by \code{n} matrix
#' } so that \eqn{\sqrt(\nu)Z/R} follows
#' a truncated Student distribution
#' @examples
#' d <- 5
#' tregress(lb =rep(-2, d), ub = rep(2, d), df = 3, n = 10,
#'   sigma = diag(0.5, d) + matrix(1, d, d))
tregress <- function (n, lb, ub, sigma, df)
{
  d = length(lb)
  if (length(ub) != d | d != sqrt(length(sigma)) | any(lb > ub)) {
    stop("lb, ub, and sigma have to match in dimension with ub>lb")
  }
  out = cholperm(sigma, lb, ub)
  Lfull = out$L
  l = out$l
  u = out$u
  D = diag(Lfull)
  perm = out$perm
  if (any(D < 1e-10)) {
    warning("Method may fail as covariance matrix is singular!")
  }
  L = Lfull/D
  u = u/D
  l = l/D
  L = L - diag(d)
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
  # assigman saddlepoint x* and mu*
  soln[d] <- exp(soln[d])
  x <- soln[1:d];
  mu <- soln[(d+1):length(soln)];
  # compute psi star
  psistar <- psyT(x= x, L = L, l = l, u = u, nu = df, mu = mu);
  # start acceptance rejection sampling
  Z <- matrix(0, nrow = d, ncol = n)
  R <- rep(0, n)
  accept <- 0L; iter <- 0L; nsim <- n
  while(accept < n){ # while # of accepted is less than n
    call <- mvtrnd(n = nsim, L = L, l = l, u = u, nu = df, mu = mu); # simulate n proposals
    idx <-  rexp(nsim) > (psistar - call$p); # acceptance tests
    m <- sum(idx)
    if(m > n - accept){
      m <- n - accept
      idx <- which(idx)[1:m]
    }
    if(m > 0){
      Z[,(accept+1):(accept+m)] <- call$Z[,idx];  # accumulate accepted
      R[(accept+1):(accept+m)] <- call$R[idx];  # accumulate accepted
    }
    accept <- accept + m; # keep track of # of accepted
    iter <- iter + 1L;  # keep track of while loop iterations
    nsim <- min(n, ceiling(nsim/m))
    if(iter == 1e3){ # if iterations are getting large, give warning
      warning('Acceptance prob. smaller than 0.001')
    } else if(iter > 1e4){ # if iterations too large, seek approximation only
      R[,1:accept]
      Z[,1:accept]
      warning('Sample of size smaller than n returned.')
    }
  }
  # # finish sampling; postprocessing
  out = sort(perm, decreasing = FALSE, index.return = TRUE)
  order = out$ix
  Z = Lfull %*% Z
  Z = Z[order, ]
  return(list(R = R, Z = t(Z)))
}
