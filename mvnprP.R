mvnprP <-
  function(n,M,l,u,mu){
    # computes P(l<X<u), where X is normal with
    # 'Prec(X)= M*M' and zero mean vector;
    # exponential tilting uses parameter 'mu';
    # Monte Carlo uses 'n' samples;
    # assumes scaled truncation bounds
    # assumes M = T
    d=length(l); # Initialization
    
    D = diag(M);
    for(k in 1:d){
      M[k,k] = 0;
    }
    Z=matrix(0,d,n); # create array for variables
    p=0;
    for (k in d:2){
      # compute matrix multiplication M*Z
      col= M[k, k:d] %*% Z[k:d,];
      # compute limits of truncation
      tl=l[k]-mu[k]+col;
      tu=u[k]-mu[k]+col;
      #simulate N(mu+nu,T_jj^-2) conditional on [tl,tu]
      Z[k,]=mu[k]+(trandn(tl,tu) - col)/D[k];
      # update likelihood ratio
      p = p+lnNpr(tl,tu)+.5*(mu[k]^2)/D[k]-mu[k]*Z[k,];
    }
    # deal with final Z(1) which need not be simulated
    col=M[1,]%*%Z;tl=l[1]+col;tu=u[1]+col;
    p=p+lnNpr(tl,tu); # update LR corresponding to Z(1)

    p=exp(p); # now switch back from logarithmic scale
    prob=mean(p);relErr=sd(p)/sqrt(n)/prob; # relative error
    est=list(prob=prob,relErr=relErr)
  }
