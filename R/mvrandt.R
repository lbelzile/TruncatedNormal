#' Random number generation from multivariate truncated Student distribution
#'
#' @param l lower bound for truncation (infinite values allowed)
#' @param u upper bound for truncation
#' @param Sig Covariance matrix
#' @param df degrees of freedom
#' @param n sample size
#'
#' @author \code{Matlab} code by Zdravko Botev, \code{R} port by Leo Belzile
#' @export
#' @references  Z. I. Botev and P. L'Ecuyer (2015), Efficient probability estimation
#' and simulation of the truncated multivariate Student-t distribution,
#' Proceedings of the 2015 Winter Simulation Conference, pp. 380-391,
#' @return a \code{d} by \code{n} matrix
#' @importFrom nleqslv nleqslv
#' @examples
#' \dontrun{
#' d <- 60L; n <- 1e3;
#' Sig <- 0.9 * matrix(1, d, d) + 0.1 * diag(d);
#' l <- (1:d)/d * 4; u <- l+2; df <- 10;
#' X <- mvrandt(l,u,Sig,df,n)
#' stopifnot(all(X>l))
#' stopifnot(all(X<u))
#' }
mvrandt <- function (l, u, Sig, df, n)
{
  d = length(l)
  if (length(u) != d | d != sqrt(length(Sig)) | any(l > u)) {
    stop("l, u, and Sig have to match in dimension with u>l")
  }
  out = cholperm(Sig, l, u)
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
  # assign saddlepoint x* and mu*
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
  return(t(sqrt(df)*t(Z)/R))
}
