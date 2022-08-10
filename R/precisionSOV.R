precisionSOV <- 
  function(SigInv, l, u){
    if(any(l>u)){
      stop("Lower bound needs to be smaller than upper bound");
    }
    
    if(length(l)!=length(u)){
      stop("Truncation limits have to be of the same length");
    }
    
    Tchol = chol(SigInv);
    
    d = length(l);
    X = rep(0,d);
    
    X[d] = mvrandn(l[d],u[d],Tchol[d,d]^(-2),1, mu = 0);
    
    for(j in ((d-1):1)){
      nu = sum(-(Tchol[j, (j+1):d] %*% X[(j+1):d])/Tchol[j,j]);
      
      
      l[j] = (l[j] - nu);
      u[j] = (u[j] - nu);
      
      
      l[j] = l[j] - nu;
      u[j] = u[j] - nu;
      
      if(j != 1){
        X[j] = mvrandn(l[j], u[j], Tchol[j,j]^(-2), 1, mu = 0);
       
      }
    }
    est = 1;
    
    u = u*diag(SigInv);
    l = l*diag(SigInv);
   
    for(k in (1:d)){
      est = est * (pnorm(u[k]) - pnorm(l[k]));
    }
    return(est);
  }


