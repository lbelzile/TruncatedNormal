cases <-
  function(p,l,u){
    x=rep(NaN,length(l))
    a=35; # threshold for switching between erfcinv and newton method
    # case 1: a<l<u
    I1=l>a;
    if (any(I1)){
      tl=l[I1]; tu=u[I1]; tp=p[I1]; x[I1]=normq(tp,tl,tu);
    }
    # case 2: l<u<-a
    J1=(u<(-a));
    if (any(J1)){
      tl=-u[J1]; tu=-l[J1]; tp=p[J1]; x[J1]=-normq(1-tp,tl,tu);
    }
    # case 3: otherwise use erfcinv
    K=!(I1|J1);
    if  (any(K)){
      tl=l[K]; tu=u[K]; tp=p[K]; 
      x[K]=Phinv(p = tp, l = tl, u = tu);
    }
    return(x)
  }
