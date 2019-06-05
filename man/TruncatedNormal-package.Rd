\name{TruncatedNormal-package}
\alias{TruncatedNormal-package}
\alias{TruncatedNormal}
\docType{package}
\title{
 Truncated Normal Distribution Toolbox
}
\description{
The routines include:
\itemize{
\item generator of \bold{independent and identically distributed} random vectors from the truncated univariate and multivariate distributions;
\item (Quasi-) Monte Carlo estimator and a \bold{deterministic upper bound}  of the cumulative distribution function of the multivariate normal and Student distributions;
\item algorithm for the accurate computation of the quantile function of the normal distribution in the extremes of its tails.
}
}
\details{
\describe{
\item{\code{\link{cholperm}}}{computes the Cholesky decomposition of a covariance matrix and returns a list containing the permuted bounds for integration. The prioritization of the variables follow the rule proposed in Gibson, Glasbey and Elston}
 \item{\code{\link{mvNcdf}}}{uses a Monte Carlo sample of size \eqn{n} to estimate the cumulative
 distribution function, Pr\eqn{( l < X < u)}, of the \eqn{d}-dimensional
 multivariate normal with zero-mean and covariance \eqn{\Sigma}, that is,
 \eqn{X} has \eqn{N(0,\Sigma)} distribution;}
\item{\code{\link{mvNqmc}}}{ provides a Quasi Monte Carlo algorithm for medium dimensions
 (say, \eqn{d<50}), in addition to the faster Monte Carlo algorithm in \code{\link{mvNcdf}};}
  \item{\code{\link{mvTqmc}}}{ provides a Quasi Monte Carlo algorithm for medium dimensions
 (say, \eqn{d<50}), in addition to the faster Monte Carlo algorithm in \code{\link{mvTcdf}};}
\item{\code{\link{mvrandn}}}{
simulates \eqn{n} \bold{identically and independently distributed} random vectors \eqn{X} from \eqn{N(0,\Sigma)}, conditional on \eqn{l<X<u};}
\item{\code{\link{mvrandt}}}{
simulates \eqn{n} \bold{identically and independently distributed} random vectors \eqn{X} from a Student distribution with scale matrix \eqn{\Sigma} and \eqn{\nu} degrees of freedom, conditional on \eqn{l<X<u};}
\item{\code{\link{norminvp}}}{ computes the quantile function at \eqn{0\le p\le 1} of the univariate \eqn{N(0,1)} distribution
 truncated to \eqn{[l,u]}, and with high precision in the tails;}
\item{\code{\link{trandn}}}{ is a fast random number generator from the univariate \eqn{N(0,1)}
 distribution truncated to \eqn{[l,u]}.}
 \item{\code{\link{mvTcdf}}}{ uses a Monte Carlo sample of size \eqn{n} to estimate the cumulative distribution function, Pr\eqn{(l < X < u)}, of the \eqn{d}-dimensional multivariate Student with zero-mean and covariance \eqn{\Sigma} and degrees of freedom \eqn{\nu}, that is, \eqn{X \sim t_\nu(0, \Sigma)}{X ~  t[\nu](0, \Sigma)}.}
\item{\code{\link{mvTqmc}}}{ provides a Quasi Monte Carlo algorithm for medium dimensions
(say, \eqn{d < 50}), in addition to the faster Monte Carlo algorithm in \code{\link{mvTcdf}}.}
\item{\code{\link{mvrandt}}}{ simulates \eqn{n} random vectors \eqn{X \sim t_\nu(0, \Sigma)}{X ~  t[\nu](0, \Sigma)}, conditional on \eqn{l < X <u}.}
\item{\code{\link{tregress}}}{ simulates \eqn{n} pairs, \eqn{Z, R}, so that \eqn{\sqrt{\nu}Z/R \sim t_\nu(0, \Sigma)}{\sqrt{\nu}Z/R ~  t[\nu]((0, \Sigma)}, conditional on \eqn{l < X <u}.}
}
}
\author{
Z. I. Botev, email: \email{botev@unsw.edu.au} and web page: \code{https://web.maths.unsw.edu.au/~zdravkobotev/}
}
\references{
\itemize{
\item Z. I. Botev (2017), \emph{The Normal Law Under Linear Restrictions:
Simulation and Estimation via Minimax Tilting}, Journal of the Royal
 Statistical Society, Series B, \bold{79} (1), pp. 1--24.
 \item Z. I. Botev and P. L'Ecuyer (2015), \emph{Efficient Estimation
 and Simulation of the Truncated Multivariate Student-t Distribution}, Proceedings of the 2015 Winter Simulation Conference,
   Huntington Beach, CA, USA
 \item Gibson G. J., Glasbey C. A., Elston D. A. (1994),
  \emph{Monte Carlo evaluation of multivariate normal integrals and sensitivity to variate ordering},
  In: Advances in Numerical Methods and Applications, pages 120--126
 }
}
