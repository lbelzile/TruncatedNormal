#' Truncated multivariate normal generator
#'
#' Simulate \eqn{n}  independent and identically distributed random vectors
#'  from the \eqn{d}-dimensional \eqn{N(0,\Sigma)} distribution
#'  (zero-mean normal with covariance \eqn{\Sigma}) conditional on \eqn{l<X<u}.
#' Infinite values for \eqn{l} and \eqn{u} are accepted.
#' @param l lower truncation limit
#' @param u upper truncation limit
#' @param Sig covariance matrix of \eqn{N(0,\Sigma)}
#' @param n number of simulated vectors
#' @details
#' \itemize{
#' \item{Bivariate normal:}{
#' Suppose we wish to simulate a bivariate \eqn{X} from \eqn{N(\mu,\Sigma)}, conditional on
#' \eqn{X_1-X_2<-6}. We can recast this as the problem of simulation
#' of \eqn{Y} from \eqn{N(0,A\Sigma A^\top)} (for an appropriate matrix \eqn{A})
#' conditional on \eqn{l-A\mu < Y < u-A\mu} and then setting \eqn{X=\mu+A^{-1}Y}.
#'    See the example code below.}
#'  \item{Exact posterior simulation for Probit regression:}{Consider the
#'      Bayesian Probit Regression model applied to the \code{\link{lupus}} dataset.
#'      Let the prior for the regression coefficients \eqn{\beta} be \eqn{N(0,\nu^2 I)}. Then, to simulate from the Bayesian
#'      posterior exactly, we first simulate
#'      \eqn{Z} from \eqn{N(0,\Sigma)}, where  \eqn{\Sigma=I+\nu^2 X X^\top,}
#'      conditional on \eqn{Z\ge 0}. Then, we simulate the posterior regression coefficients, \eqn{\beta}, of the Probit regression
#'      by drawing \eqn{(\beta|Z)} from \eqn{N(C X^\top Z,C)}, where \eqn{C^{-1}=I/\nu^2+X^\top X}.
#' See the example code below.}
#' }
#' @return a \eqn{d} by \eqn{n} matrix storing the random vectors, \eqn{X}, drawn from \eqn{N(0,\Sigma)}, conditional on \eqn{l<X<u};
#' @note The algorithm may not work or be very inefficient if \eqn{\Sigma} is close to being rank deficient.
#' @seealso \code{\link{mvNqmc}}, \code{\link{mvNcdf}}
#' @export
#' @examples
#'  # Bivariate example.
#'
#'  Sig <- matrix(c(1,0.9,0.9,1), 2, 2);
#'  mu <- c(-3,0); l <- c(-Inf,-Inf); u <- c(-6,Inf);
#'  A <- matrix(c(1,0,-1,1),2,2);
#'  n <- 1e3; # number of sampled vectors
#'  Y <-mvrandn(l - A %*% mu, u - A %*% mu, A %*% Sig %*% t(A), n);
#'  X <- rep(mu, n) + solve(A, diag(2)) %*% Y;
#'  # now apply the inverse map as explained above
#'  plot(X[1,], X[2,]) # provide a scatterplot of exactly simulated points
#' # Exact Bayesian Posterior Simulation Example.
#'
#' data("lupus"); # load lupus data
#' Y = lupus[,1]; # response data
#' X = lupus[,-1]  # construct design matrix
#' m=dim(X)[1]; d=dim(X)[2]; # dimensions of problem
#'  X=diag(2*Y-1) %*%X; # incorporate response into design matrix
#'  nu=sqrt(10000); # prior scale parameter
#'  C=solve(diag(d)/nu^2+t(X)\%*\%X);
#'  L=t(chol(t(C))); # lower Cholesky decomposition
#'  Sig=diag(m)+nu^2*X %*% t(X); # this is covariance of Z given beta
#'  l=rep(0,m);u=rep(Inf,m);
#'  est=mvNcdf(l,u,Sig,1e3);
#'  # estimate acceptance probability of Crude Monte Carlo
#'  print(est$upbnd/est$prob)
#'  # estimate the reciprocal of acceptance probability
#'  n=1e4 # number of iid variables
#'  z=mvrandn(l,u,Sig,n);
#'  # sample exactly from auxiliary distribution
#'  beta=L %*% matrix(rnorm(d*n),d,n)+C %*% t(X) %*% z;
#'  # simulate beta given Z and plot boxplots of marginals
#'  boxplot(t(beta))
#'  # plot the boxplots of the marginal
#'  # distribution of the coefficients in beta
#'  print(rowMeans(beta)) # output the posterior means
mvrandn <-  function(l, u, Sig, n){
    d <- length(l); # basic input check
    if(length(u) != d | d != sqrt(length(Sig)) | any(l > u)){
      stop('l, u, and Sig have to match in dimension with u>l')
    }
    if(d == 1){
      std.dev <- sqrt(Sig[1]) #if Sigma not declared as matrix
      return(std.dev * trandn(rep(l/std.dev, n), rep(u/std.dev, n)))
    }
    # Cholesky decomposition of matrix
    out <- cholperm(Sig,l,u);
    Lfull=out$L;l=out$l;u=out$u;D=diag(Lfull);perm=out$perm;
    if (any(D<1e-10)){
      warning('Method may fail as covariance matrix is singular!')
    }
    L=Lfull/D;u=u/D;l=l/D; # rescale
    L=L-diag(d); # remove diagonal
    # find optimal tilting parameter via non-linear equation solver
    xmu<-nleq(l,u,L) # nonlinear equation solver
    x=xmu[1:(d-1)];mu=xmu[d:(2*d-2)]; # assign saddlepoint x* and mu*
    # compute psi star
    psistar=psy(x,L,l,u,mu);
    # start acceptance rejection sampling
    iter=0; rv=c();
    repeat{
      out=mvnrnd(n,L,l,u,mu);logpr=out$logpr;Z=out$Z; # simulate n proposals
      idx=-log(runif(n))>(psistar-logpr); # acceptance tests
      rv=cbind(rv,Z[,idx]);  # accumulate accepted
      accept=dim(rv)[2]; # keep track of # of accepted
      iter=iter+1;  # keep track of while loop iterations
      if (iter==1e3){ # if iterations are getting large, give warning
        warning('Acceptance prob. smaller than 0.001')
      } else if (iter>1e4){ # if iterations too large, seek approximation only
        accept=n;rv=cbind(rv,Z); # add the approximate samples
        warning('Sample is only approximately distributed.')
      }
      if (accept>=n){# if # of accepted is less than n
        break
      }
    }
    # finish sampling; postprocessing
    out=sort(perm,decreasing = FALSE,index.return = TRUE);order=out$ix;
    rv=rv[,1:n]; # cut-down the array to desired n samples
    rv=Lfull%*%rv; # reverse scaling of L
    rv=rv[order,]; # reverse the Cholesky permutation
    return(rv)
  }
