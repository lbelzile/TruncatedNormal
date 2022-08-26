psyP <- function(x, M, l, u, mu)
  {
  # psi(x,mu) function adapted to back-substituted precision matrix
  # assume M = T with diagonal
  # assume scaled truncation bounds
  d = length(l);
  nu = rep(0,d);
  x[1] = 0;
  mu[1] = 0;
  
  D = diag(M);
  for(k in 1:d){
    M[k,k] = 0;
  }
  
  nu[-d] = M[-d,] %*% x;
  lt = l + nu - mu;
  ut = u + nu - mu;
  
  return(sum(lnNpr(lt, ut) + 0.5*(mu^2)/D - x*mu));
  }