#' Fast truncated normal generator
#'
#' Efficient state-of-the-art generator of a vector of \code{length(l)=length(u)}
#'  from the standard multivariate normal distribution truncated over the region \eqn{[l,u]}.
#'  Infinite values for \code{u} and \code{l} are accepted.
#'
#' @param l lower truncation limit
#' @param u upper truncation limit
#' @details
#'  Suppose we wish to simulate a random variable \eqn{Z} drawn from \eqn{N(\mu,\sigma^2)} and
#'  conditional on \eqn{l<Z<u} using the inverse transform method.
#'  To achieve this, first compute
#'  \code{X=norminvp(runif(1),(l-mu)/sig,(u-mu)/sig)} and then set
#'  \code{Z=mu+sig*X}
#' @return random variable drawn from the truncated normal distribution
#' @note  Use \code{\link{norminvp}} for the (slower) inverse transform method of simulating truncated normal variables.
#' @seealso \code{\link{norminvp}}
#' @export
#' @keywords internal
#' @author Zdravko I. Botev
#' @references Z. I. Botev (2017), \emph{The Normal Law Under Linear Restrictions:
#' Simulation and Estimation via Minimax Tilting}, Journal of the Royal
#' Statistical Society, Series B, \bold{79} (1), pp. 1--24.
#'
#' @examples
#' trandn(l = 1,u = Inf)
#' trandn(l = rep(1, 10), u = rep(Inf, 10))
trandn <-  function(l, u){
     if(any(l>u)){
      stop('Truncation limits have to be vectors of the same length with l<u')
    }
    x=rep(0,length(l));
    a=.4; # threshold for switching between methods
    # threshold can be tuned for maximum speed for each Matlab version
    # three cases to consider:
    # case 1: a<l<u
    I <- l>a;
    if (any(I)){
      x[I] <- ntail(l[I], u[I]);
    }
    # case 2: l<u<-a
    J <- u < (-a);
    if (any(J)){
      x[J] <- -ntail(-u[J],-l[J]);
    }
    # case 3: otherwise use inverse transform or accept-reject
    L <- !(I|J);
    if (any(L)){
      x[L] <- tn(l[L],u[L]);
    }
    return(x)
  }
