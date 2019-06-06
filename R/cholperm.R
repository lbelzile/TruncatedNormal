#' Cholesky decomposition for Gaussian distribution function with permutation
#' 
#' This function computes the Cholesky decomposition of a covariance matrix
#' \code{Sigma} and returns a list containing the permuted bounds for integration. 
#' The prioritization of the variables follows either the rule proposed in Gibson, Glasbey and Elston (1994),
#' reorder variables to have outermost variables with smallest expected values. The alternative is the scheme proposed
#' in Genz and Bretz (2009) that minimizes the variance of the truncated Normal variates.
#' 
#' The list contains an integer vector \code{perm} with the indices of the permutation, which is such that
#' \code{Sigma(perm, perm) == L \%*\% t(L)}.
#' The permutation scheme is described in Genz and Bretz (2009) in Section 4.1.3, p.37.
#' @param Sigma \code{d} by \code{d} covariance matrix
#' @param l \code{d} vector of lower bounds
#' @param u \code{d} vector of upper bounds
#' @param method string indicating which method to use. Default to \code{"GGE"}
#' @return a list with components
#' \itemize{
#' \item{\code{L}: }{Cholesky root}
#' \item{\code{l}: }{permuted vector of lower bounds}
#' \item{\code{u}: }{permuted vector of upper bounds}
#' \item{\code{perm}: }{vector of integers with ordering of permutation}
#' }
#' @export
#' @references Genz, A. and Bretz, F. (2009). Computations of Multivariate Normal and t Probabilities, volume 105. Springer, Dordrecht.
#' @references Gibson G.J., Glasbey C.A. and D.A. Elton (1994).  Monte Carlo evaluation of multivariate normal integrals and sensitivity to variate ordering. In: Dimon et al., Advances in Numerical Methods and Applications, WSP, pp. 120-126.
cholperm <- function(Sigma, l, u, method = c("GGE", "GB")){
  method <- match.arg(method, c("GGE","GB"))[1]
  if(method == "GGE"){
    .cholpermGGE(Sigma = Sigma, l = l, u = u)
  } else{
    .cholpermGB(Sigma = Sigma, l = l, u = u)
  }
}