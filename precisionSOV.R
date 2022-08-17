precisionSOV <- 
  function(SigInv, l, u, n = 10000){
    if(any(l>u)){
      stop("Lower bound needs to be smaller than upper bound");
    }
    
    if(length(l)!=length(u)){
      stop("Truncation limits have to be of the same length");
    }
    
    Tchol = chol(SigInv);
    
    d = length(l);
    
    for(k in 1:d){
      l[k] = Tchol[k,k]*l[k];
      u[k] = Tchol[k,k]*u[k];
    }
    a = l;
    b = u;
    acc = 0;
    
    for(m in 1:n){
      
      x = rep(0,d);
      nu = rep(0,d);
      l = a;
      u = b;
    
      x[d] = trandn(l[d],u[d])/Tchol[d,d];
      
      for(j in ((d-1):1)){
        for(i in (j+1):d){
          nu[j] = nu[j] + Tchol[i,j]*x[i];
        }
        
        l[j] = l[j] + nu[j];
        u[j] = u[j] + nu[j];
        
        if(j != 1){
          x[j] = (trandn(l[j],u[j])-nu[j])/Tchol[j,j];
        }
      }
      est = 1;
      for(k in (1:d)){
        est = est * (pnorm(u[k]) - pnorm(l[k]));
      }
      
      acc = acc + est;
    }
    return(acc/n);
  }


