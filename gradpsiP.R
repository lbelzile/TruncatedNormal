gradpsiP <- function(y, M, l, u)
{
  # implements grad_psi(x) (precision matrix) to find 
  # optimal exponential twisting;
  # y = (x,mu);
  # assumes scaled bounds T_jj*l_j and T_jj*u_j;
  # assume M = T;
  d <- length(u);
  nu <- rep(0,d);
  x <- nu; 
  mu <- nu;
  x[2:d] <- y[1:(d-1)];
  mu[2:d] <- y[d:(2*d-2)]
  #print(mu);
  D = diag(M);
  for(k in 1:d){
    M[k,k] = 0;
  }
  
  # compute now ~l and ~u
  nu[-d] <- M[-d,] %*% x;
  lt <- l - mu + nu;
  ut <- u - mu + nu;
  # compute gradients avoiding catastrophic cancellation
  w <- lnNpr(lt, ut);
  pl <- exp(-0.5*lt^2-w)/sqrt(2*pi);
  pu <- exp(-0.5*ut^2-w)/sqrt(2*pi);
  P <- pl-pu;
  # output the gradient
  dfdx <- -mu[-1] - as.vector(crossprod(M[,-1],P)); #as.vector(crossprod(P, M[,-1]);
  dfdm <- mu/D - x + P;
  grad <- c(dfdx, dfdm[-1])
  # here compute Jacobian matrix
  grad
}


jacpsiP <-  function(y, M, l, u){ 
  # implements grad_psi(x) to find optimal exponential twisting;
  # y = (x,mu);
  # assumes scaled bounds T_jj*l_j and T_jj*u_j;
  # assume M = T
  d <- length(u);
  nu <- rep(0,d);
  x <- nu; 
  mu <- nu;
  x[2:d] <- y[1:(d-1)];
  mu[2:d] <- y[d:(2*d-2)]
  
  D = diag(M);
  for(k in 1:d){
    M[k,k] = 0;
  }
  
  # compute now ~l and ~u
  nu[-d] <- M[-d,] %*% x;
  lt <- l - mu + nu;
  ut <- u - mu + nu;
  # compute gradients avoiding catastrophic cancellation
  w <- lnNpr(lt, ut);
  pl <- exp(-0.5*lt^2-w)/sqrt(2*pi);
  pu <- exp(-0.5*ut^2-w)/sqrt(2*pi)
  P <- pl-pu;
  
  # here compute Jacobian matrix
  lt[is.infinite(lt)] <- 0 
  ut[is.infinite(ut)] <- 0
  
  dP <- -(P^2) + lt * pl - ut * pu # dPdm
  DM <- rep(dP,1,d) * M
  mx <- -diag(1/D) - DM
  xx <- crossprod(t(M), DM)
  mx <- mx[2:d,2:d]
  xx <- xx[2:d,2:d]
  if (d>2){
    Jac <- rbind(cbind(xx, t(mx)), cbind(mx, diag((1/D)[2:d] + dP[2:d])))
  } else {
    Jac <- rbind(cbind(xx, t(mx)), cbind(mx, (1/D)[2:d] + dP[2:d]))
  }
  Jac
}
