gradpsiT <- function(y, L, l, u, nu)
{
  # implements gradient of psi(x) to find optimal exponential twisting;
  # assumes scaled 'L' with zero diagonal;
  d <- length(u)
  co <- x <- mu <- rep(0, d)
  x[-d] <- y[1:(d-1)]
  r <- exp(y[d])
  mu[-d] <- y[(d+1):(2*d-1)]
  eta <- y[2*d]
  l <- l/sqrt(nu);
  u <- u/sqrt(nu);
  # compute now ~l and ~u
  co[-1] <- L[2:d,] %*% x;
  lt <- r * l - mu - co;
  ut <- r * u - mu - co;
  # compute gradients avoiding catastrophic cancellation
  w <- lnNpr(lt, ut);
  pl <- exp(-0.5 * lt^2 - w) / sqrt(2 * pi);
  pu <- exp(-0.5 * ut^2 - w) / sqrt(2 * pi);
  P <- pl - pu;
  # output the gradient
  dfdx <- - mu[1:(d-1)] + t(t(P) %*% L[,1:(d-1)]);
  dfdm <-  mu - x + P;
  l[is.infinite(l)] <- 0;
  u[is.infinite(u)] <- 0
  dfdr <- (nu-1) / r - eta + sum(u * pu - l * pl);
  dfde <- eta - r + exp(dnorm(eta, log = TRUE) - lnNpr(-eta, Inf));
  c(dfdx, dfdr, dfdm[1:(d-1)], dfde)
}
