precisionSOV <- 
  function(SigInv, l, u){
    if(any(l>u)){
      stop("Lower bound needs to be smaller than upper bound");
    }
    
    if(length(l)!=length(u)){
      stop("Truncation limits have to be of the same length");
    }
    
    Tchol = t(chol(SigInv));
    d = length(l);
    X = rep(0,d);
    rem = 0;
    
    l[d] = Tchol[d,d]*l[d];
    u[d] = Tchol[d,d]*u[d];
    X[d] = trandn(l[d],u[d]);
    
    for(j in ((d-1):1)){
      l[j] = l[j]*Tchol[j,j];
      u[j] = u[j]*Tchol[j,j];
      
      for(i in ((j+1):d)){
        rem = rem + (Tchol[i,j]*(X[i] - eta(i,d,Tchol,X[i])))/Tchol[i,i];
      }
      
      l[j] = l[j] + rem;
      u[j] = u[j] + rem;
      
      X[j] = trandn(l[j], u[j]);
    }
    est = 1;
    for(k in (1:d)){
      est = est * (pnorm(u[k]) - pnorm(l[k]));
    }
    return(est);
  }

eta <-
  function(j, d, M, x){
    eta = 0;
    if(j == d){
      return(0);
    }
    for(i in ((j+1):d)){
      eta = eta + M[i,j]/M[i,i];
    }
    eta = (eta * x)/M[j,j];
    return(eta);
  }
