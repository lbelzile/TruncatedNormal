cases <-
  function(p,l,u){
    x=rep(NaN,length(l))
    a=35; # threshold for switching between erfcinv and newton method
    # case 1: a<l<u
    I=l>a;
    if (any(I)){
      tl=l[I]; tu=u[I]; tp=p[I]; x[I]=normq(tp,tl,tu);
    }
    # case 2: l<u<-a
    J=(u<(-a));
    if (any(J)){
      tl=-u[J]; tu=-l[J]; tp=p[J]; x[J]=-normq(1-tp,tl,tu);
    }
    # case 3: otherwise use erfcinv
    K=!(I|J);
    if  (any(K)){
      tl=l[K]; tu=u[K]; tp=p[K]; 
      x[K]=Phinv(p = tp, l = tl, u = tu);
    }
    return(x)
  }
