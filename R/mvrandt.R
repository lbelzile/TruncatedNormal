#' Random number generation from multivariate truncated Student distribution
#'
#' @param l lower bound for truncation (infinite values allowed)
#' @param u upper bound for truncation
#' @param Sig covariance matrix
#' @param df degrees of freedom
#' @param n sample size
#' @param mu location parameter
#'
#' @author \code{Matlab} code by Zdravko Botev, \code{R} port by Leo Belzile
#' @export
#' @references  Z. I. Botev and P. L'Ecuyer (2015), Efficient probability estimation
#' and simulation of the truncated multivariate Student-t distribution,
#' Proceedings of the 2015 Winter Simulation Conference, pp. 380-391,
#' @return a \code{d} by \code{n} matrix
#' @importFrom nleqslv nleqslv
#' @importFrom stats "pt" "qt"
#' @keywords internal
#' @examples
#' \dontrun{
#' d <- 60L; n <- 1e3;
#' Sig <- 0.9 * matrix(1, d, d) + 0.1 * diag(d);
#' l <- (1:d)/d * 4; u <- l+2; df <- 10;
#' X <- mvrandt(l,u,Sig,df,n)
#' stopifnot(all(X>l))
#' stopifnot(all(X<u))
#' }
mvrandt <- function (l, u, Sig, df, n, mu = NULL)
{
  d = length(l)
  if (length(u) != d | d != sqrt(length(Sig)) | any(l > u)) {
    stop("l, u, and Sig have to match in dimension with u>l")
  }
  if(!is.null(mu)){
    l <- l - mu
    u <- u - mu
  }
  if (d == 1){ #Univariate case, via inverse CDF method
    std.dev <- sqrt(Sig[1])
    #Inverse CDF method
    if(l > 0){
      if(is.null(mu)){
        return( std.dev * (-qt( pt(l/std.dev, df = df, lower.tail = FALSE) - runif(n) * 
                                  (pt(l/std.dev, df = df, lower.tail = FALSE) - pt(u/std.dev, 
                                                                                   df = df, lower.tail = FALSE)), df = df))) 
      } else{
        return( std.dev * (-qt( pt(l/std.dev, df = df, lower.tail = FALSE) - runif(n) * 
                                (pt(l/std.dev, df = df, lower.tail = FALSE) - pt(u/std.dev, 
                                                                                 df = df, lower.tail = FALSE)), df = df)) + mu)
      }
    } else{
      if(is.null(mu)){
        return(std.dev * (qt(pt(l/std.dev, df = df) + runif(n) * 
                             (pt(u/std.dev, df = df) - pt(l/std.dev, df = df)), df = df)))         
      } else{
        return(std.dev * (qt(pt(l/std.dev, df = df) + runif(n) * 
                               (pt(u/std.dev, df = df) - pt(l/std.dev, df = df)), df = df)) + mu)     
      }
    }
  }
  
  out <- cholperm(Sig, l, u)
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
                              global = "pwldog", method = "Broyden", control = list(maxit = 500L))
  
  soln <- solvneq$x
  #fval <- solvneq$fvec
  exitflag <- solvneq$termcd
  if(!(exitflag %in% c(1,2)) || !isTRUE(all.equal(solvneq$fvec, rep(0, length(x0)), tolerance = 1e-6))){
    warning('Did not find a solution to the nonlinear system in `mvrandt`!')
  }
  # assign saddlepoint x* and mu*
  soln[d] <- exp(soln[d])
  x <- soln[1:d];
  #TODO check that the solution lies in convex set
  #Problem is: we did not generate Zd
  # mid <- sqrt(df) * L[,-d] %*% x[-d]
  # if(any(x[d] * l[-d] > mid, x[d] * u[-d] < mid)){
  #   warning("Solution of nonlinear system does not satisfy convex optimization problem constraints")
  # }
  muV <- soln[(d+1):length(soln)];
  # compute psi star
  psistar <- psyT(x= x, L = L, l = l, u = u, nu = df, mu = muV);
  # start acceptance rejection sampling
  Z <- matrix(0, nrow = d, ncol = n)
  R <- rep(0, n)
  accept <- 0L; iter <- 0L; nsim <- n; ntotsim <- 0L
  while(accept < n){ # while # of accepted is less than n
    call <- mvtrnd(n = nsim, L = L, l = l, u = u, nu = df, mu = muV); # simulate n proposals
    ntotsim <- ntotsim + nsim
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
    if((ntotsim > 1e4) && (accept / ntotsim < 1e-3)){ # if iterations are getting large, give warning
      warning('Acceptance probability smaller than 0.001')
    } else if(iter > 1e5){ # if iterations too large, seek approximation only
      if(accept == 0){
        stop("Could not sample from truncated Student - check input")
      } else if(accept > 1){
        R <- R[1:accept]
        Z <- Z[,1:accept]
        warning('Sample of size smaller than n returned.')
      }
    }
  }
  # # finish sampling; postprocessing
  out = sort(perm, decreasing = FALSE, index.return = TRUE)
  order = out$ix
  Z = Lfull %*% Z
  Z = Z[order, ]
  #Add back mean only if non-zero
  if(!is.null(mu)){
    return(t(sqrt(df)*t(Z)/R) + mu)
  } else{
    return(t(sqrt(df)*t(Z)/R)) 
  }
}
