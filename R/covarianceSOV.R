covarianceSOV <-
  function(Sig, l, u){
    if(any(l>u)){
      stop("Lower bound needs to be smaller than upper bound");
    }
    
    if(length(l)!=length(u)){
      stop("Truncation limits have to be of the same length");
    }
    
    Lchol = t(chol(Sig));
    d = length(l);
    X = rep(0,d);
    rem = 0;
    
    l[1] = l[1]/Lchol[1,1];
    u[1] = u[1]/Lchol[1,1];
    
    X[1] = trandn(l[1], u[1]);
    
    for(j in (2:d)){
      for(i in (1:(j-1))){
        rem = rem + Lchol[j,i]*X[i];
      }
      l[j] = (l[j] - rem)/Lchol[j,j];
      u[j] = (u[j] - rem)/Lchol[j,j];
      
      X[j] = trandn(l[j],u[j]);
    }
    est = 1;
    for(k in (1:d)){
      est = est * (pnorm(u[k]) - pnorm(l[k]));
    }
    return(est);
  }