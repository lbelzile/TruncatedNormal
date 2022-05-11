library(TruncatedNormal)
lb <- c(0, 0)
ub <- c(740.0, 76.2)
mu <- c(344.31293403, 62.6937066)
sigma <- matrix(c(36407.0005966, -1167.50805662, -1167.50805662, 290.76915744), 2, 2)
df1 <- 3
df2 <- 300
x = c(100, 50)
set.seed(1234)

## "Truncated normal DF (MC versus QMC) give similar answers" 

  tinytest::expect_equal(ptmvnorm(x, mu=mu, sigma=sigma, lb=lb, ub=ub, log=FALSE, type="qmc"),
                         ptmvnorm(x, mu=mu, sigma=sigma, lb=lb, ub=ub, log=FALSE, type="mc"), 
                         tolerance = 1e-4)
  tinytest::expect_equal(ptmvt(x, B = 1e5, sigma=sigma, df = 2, lb=lb, ub=ub, log=FALSE, type="mc"),
                         ptmvt(x, B = 1e5, sigma=sigma, lb=lb, ub=ub, df = 2, log=FALSE, type="qmc"), 
                         tolerance = 3e-3)


## "Student with large df gives same answer as normal"
  
  tinytest::expect_equal(ptmvnorm(x, B = 1e6, sigma=sigma, lb=lb, ub=ub, log=FALSE, type="qmc"),
                         ptmvt(x, B = 1e6, sigma=sigma, lb=lb, ub=ub, df = 300, log=FALSE, type="qmc"), 
                         tolerance = 2e-3)



mean <- rep(0, 5)
lower <- rep(-1, 5)
upper <- rep(3, 5)
corr <- matrix(0.5, 5, 5) + diag(0.5, 5)
prob <- pmvnorm(lb = lower, ub = upper, mu = mean, sigma = corr)


## "Univariate probabilities",{
  
  tinytest::expect_equivalent(pmvnorm(lb = -Inf, ub = 3, mu = 2, sigma = 1), pnorm(3, mean = 2))
  tinytest::expect_equivalent(pmvt(lb = -Inf, ub = 3, df = 2, mu = 0, sigma = 1), pt(3, 2))


mean_tnorm <- function(lb, ub, mu, sigma){
  # Mean of independent univariate truncated Normal distribution - vectorized
  stopifnot(length(lb) == length(ub))
  mu <- rep(mu, length.out = length(lb))
  sigma <- rep(sigma, length.out = length(lb))
  mu + sigma*(dnorm((lb-mu)/sigma)-dnorm((ub-mu)/sigma))/
    (pnorm((ub-mu)/sigma) - pnorm((lb-mu)/sigma))
}

mean_tt <- function(lb, ub, df){
  ((df + lb^2)^(-(df-1)/2) -  (df + ub^2)^(-(df-1)/2)) * exp(lgamma(0.5*(df-1)) + 0.5*df*log(df) - lgamma(0.5*df) - lgamma(0.5))/(2*(pt(ub, df) - pt(lb, df)))
}

B <- 1e5
D <- 10
muV <- 1:D
Smat <- diag(0.5, D) + matrix(0.5, D, D)
## "Expectation of (truncated) elliptical distributions"
  
  tinytest::expect_equal(colMeans(rtmvnorm(n = B, sigma = Smat)),
                         rep(0, D), tolerance = 5/sqrt(B))
  tinytest::expect_equal(colMeans(rtmvnorm(n = B, mu = muV, sigma = Smat)),
                         muV, tolerance = 5/sqrt(B))
  tinytest::expect_equal(colMeans(rtmvnorm(n = B, lb = rep(0, D), ub = rep(2*D, D), mu = muV, sigma = diag(D, D))),
                         mean_tnorm(lb = rep(0, D), ub = rep(2*D, D), mu = muV, sigma = sqrt(D)),
                         tolerance = 5/sqrt(B))
  tinytest::expect_equal(colMeans(rtmvt(n = B, lb = (1:D)/D, df = 3, ub = 2*(1:D)/D, mu = rep(0,D), sigma = diag(1, D))),
                         mean_tt(lb = (1:D)/D, ub = 2*(1:D)/D, df = 3),
                         tolerance = 5/sqrt(B))
  


lb <- rnorm(n = D, mean = 0, sd = 10)
ub <- lb + rgamma(D, shape = 4, rate = 1)
## "Bounds of simulated variables"
  
  tinytest::expect_true(isTRUE(all(apply(rtmvnorm(n = 1e4, lb = lb, ub = ub, mu = muV, 100*Smat), 2, min) > lb)))
  tinytest::expect_true(isTRUE(all(apply(rtmvnorm(n = 1e4, lb = lb, ub = ub, mu = muV, 100*Smat), 2, min) < ub)))
  tinytest::expect_true(isTRUE(all(apply(rtmvt(n = 1e4, df = 3, lb = lb, ub = ub, mu = muV, 100*Smat), 2, min) > lb)))
  tinytest::expect_true(isTRUE(all(apply(rtmvt(n = 1e4, df = 3, lb = lb, ub = ub, mu = muV, 100*Smat), 2, min) < ub)))
  


## "Bounds on distribution function beyond truncation points"
  
  tinytest::expect_equal(ptmvnorm(q = ub + runif(D), lb = lb, ub = ub, mu = muV, Smat), 1)
  tinytest::expect_equal(ptmvt(q = ub + runif(D), df = 2, lb = lb, ub = ub, mu = muV, Smat), 1)
  tinytest::expect_equal(ptmvnorm(q = lb + c(-1, runif(D-1)), lb = lb, ub = ub, mu = muV/D, Smat), 0)
  tinytest::expect_equal(ptmvt(q = lb + c(-1, runif(0, D -1)), df = 2, lb = lb, ub = ub, mu = muV/D, Smat), 0)


pt <- rnorm(D)
## "Untruncated density agrees with that in the mvtnorm package"
  
  tinytest::expect_equal(mvtnorm::dmvnorm(x = pt, mean = muV, sigma = Smat, log = TRUE), 
                         TruncatedNormal::dtmvnorm(x = pt, mu = muV, lb = rep(-Inf, D), ub = rep(Inf, D), sigma = Smat, log = TRUE))
  tinytest::expect_equal(mvtnorm::dmvt(x = pt, df = 2, delta = muV, sigma = Smat, log = TRUE), 
                         TruncatedNormal::dtmvt(x = pt, df = 2, mu = muV, lb = rep(-Inf, D), ub = rep(Inf, D), sigma = Smat, log = TRUE))



pt <- rnorm(D)
## "Untruncated density agrees with that in the mvtnorm package"
  
  tinytest::expect_equal(mvtnorm::dmvnorm(x = pt, mean = muV, sigma = Smat, log = TRUE), 
                         TruncatedNormal::dtmvnorm(x = pt, mu = muV, lb = rep(-Inf, D), ub = rep(Inf, D), sigma = Smat, log = TRUE))
  tinytest::expect_equal(mvtnorm::dmvt(x = pt, df = 2, delta = muV, sigma = Smat, log = TRUE), 
                         TruncatedNormal::dtmvt(x = pt, df = 2, mu = muV, lb = rep(-Inf, D), ub = rep(Inf, D), sigma = Smat, log = TRUE))


d <- 15
sigma <- 0.5 * (diag(d) + matrix(1, d, d))
## "Known probability"
  
  tinytest::expect_equivalent((d+1)*pmvnorm(sigma = sigma, lb = rep(0, d), type = "qmc", B = B),
  1, tolerance = 1/sqrt(B))
  
