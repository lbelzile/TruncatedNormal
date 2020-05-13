gradpsi <-  function(y, L, l, u){ 
    # implements grad_psi(x) to find optimal exponential twisting;
    # assumes scaled 'L' with zero diagonal;
    d <- length(u);
    cv <- rep(0,d);
    x <- cv; 
    mu <- cv
    x[1:(d-1)] <- y[1:(d-1)];
    mu[1:(d-1)] <- y[d:(2*d-2)]
    # compute now ~l and ~u
    cv[-1] <- L[-1,] %*% x;
    lt <- l - mu - cv;
    ut <- u - mu - cv;
    # compute gradients avoiding catastrophic cancellation
    w <- lnNpr(lt, ut);
    pl <- exp(-0.5*lt^2-w)/sqrt(2*pi);
    pu <- exp(-0.5*ut^2-w)/sqrt(2*pi)
    P <- pl-pu;
    # output the gradient
    dfdx <- -mu[-d] + as.vector(crossprod(P, L[,-d]))
    dfdm <- mu - x + P
    grad <- c(dfdx, dfdm[-d])
    # here compute Jacobian matrix
    grad
  }

jacpsi <-  function(y, L, l, u){ 
  # implements grad_psi(x) to find optimal exponential twisting;
  # assume scaled 'L' with zero diagonal;
  d <- length(u);
  cv <- rep(0,d);
  x <- cv; 
  mu <- cv
  x[1:(d-1)] <- y[1:(d-1)];
  mu[1:(d-1)] <- y[d:(2*d-2)]
  # compute now ~l and ~u
  cv[-1] <- L[-1,] %*% x;
  lt <- l - mu - cv;
  ut <- u - mu - cv;
  # compute gradients avoiding catastrophic cancellation
  w <- lnNpr(lt, ut);
  pl <- exp(-0.5*lt^2-w)/sqrt(2*pi);
  pu <- exp(-0.5*ut^2-w)/sqrt(2*pi)
  P <- pl-pu;
 
  # here compute Jacobian matrix
  lt[is.infinite(lt)] <- 0
  ut[is.infinite(ut)] <- 0
  dP <- (-P^2) + lt * pl - ut * pu # dPdm
  DL <- rep(dP,1,d) * L
  mx <- -diag(d) + DL
  xx <- crossprod(L, DL)
  mx <- mx[1:(d-1),1:(d-1)]
  xx <- xx[1:(d-1),1:(d-1)]
  if (d>2){
    Jac <- rbind(cbind(xx, t(mx)), cbind(mx, diag(1+dP[1:(d-1)])))
  } else {
    Jac <- rbind(cbind(xx, t(mx)), cbind(mx, 1 + dP[1:(d-1)]))
  }
  Jac
}
