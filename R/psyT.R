psyT <- function(x, L, l, u, nu, mu)
{
  #implements psi(x,mu); assumes scaled 'L' without diagonal;
  d <- length(u);
  r <- x[d]
  eta <- mu[d]
  x[d] <- 0; mu[d] <- 0;
  l <- l/sqrt(nu);
  u <- u/sqrt(nu);
  # compute now ~l and ~u
  c <- L %*% x;
  l <- r * l - mu - c;
  u <- r * u - mu - c;
  sum(lnNpr(l, u)) + 0.5 * sum(mu * mu) - sum(x * mu) +
    log(2 * pi) / 2 - lgamma(nu / 2) - (0.5 * nu - 1) * log(2) +
    0.5 * sum(eta * eta) - r * eta + (nu - 1) * log(r) +
    lnNpr(-eta, Inf);
}
