#Holds functions for check loss (aggregate)

require(MASS)
require(far)
require(esaBcv)
require(FACTMLE)
require(POET)
require(dpglasso)


rmt.est<- function(K,S,m,type){
  
  n<- dim(S)[1]
  spd<- eigen(S,symmetric = TRUE)
  rho<- n/(m-1)
  
  if(type==0){
    
    l.hat<- spd$values[1:K]
    l0.hat<-spd$values[K+1]
    zeta<-matrix(0,K,1)
    for(k in 1:K){
      
      temp<-l.hat[k]/l0.hat
      a<- rho/(temp-1)^2
      b<- rho/(temp-1)
      zeta[k]<- sqrt((1-a)/(1+b))
    }
  }
  
  if(type==1){
    l.tilde<-spd$values
    l0.tilde<- ((n-K)^{-1})*sum(l.tilde[(K+1):n])
    
    out<- eigen.est(l.tilde,l0.tilde,rho,K)
    l0.hat<-out$l0.hat
    l.hat<-out$l.hat
    
    zeta<-matrix(0,K,1)
    for(k in 1:K){
      
      temp<-l.hat[k]/l0.hat
      a<- rho/(temp-1)^2
      b<- rho/(temp-1)
      zeta[k]<- sqrt((1-a)/(1+b))
    }
  }
  if(type==2){#for cv
    
    l.hat<- spd$values[1:K]
    l0.hat<-spd$values[K+1]
    zeta<-matrix(1,K,1)
  }
  
  return(list("l0.hat"=l0.hat,"l.hat"=l.hat,"zeta"=zeta,
              "pj"=spd$vectors,"K"=K))
  
}
S.est.naive.1<-function(S,K){
  
  n<- dim(S)[1]
  spd<- eigen(S)
  l.tilde<-spd$values
  l0.tilde<- ((n-K$numOfSpikes)^{-1})*sum(l.tilde[(K$numOfSpikes+1):n])#K$sigma.sq
  eigvec<- spd$vectors
  
  # lam.mat<-diag(c(l.tilde[1:K$numOfSpikes],rep(l0.tilde,n-K$numOfSpikes)))
  # S.naive<- pj%*%lam.mat%*%t(pj)
  
  S.naive<- matrix(0,n,n)
  temp<- S.naive
  for(k in 1:K$numOfSpikes){

    temp<- eigvec[,k]%*%t(eigvec[,k])+temp
    S.naive<- l.tilde[k]*(eigvec[,k]%*%t(eigvec[,k]))+S.naive

  }
  S.naive<- S.naive+l0.tilde*(diag(n)-temp)
  #S.naive<- S.naive+l0.tilde*diag(n)
  
  return(S.naive)
}

eigen.est<-function(l,l.0,rho,K){
  
  x<- matrix(0,K,1)
  b<-l.0
  
  for(k in 1:K){
    a<- l[k]
    temp.1<-b+a-rho*b
    x[k]<- (temp.1+sqrt((temp.1^2)-4*a*b))/(2)
  }
  
  temp.2<- ((x/l.0)-1)^{-1}
  xi<-K+sum(temp.2[1:K])
  
  l0.hat<- l.0*(1+rho*xi/(n-K))
  l.hat<-matrix(0,K,1)
  for(k in 1:K){
    
    a<- l[k]/l0.hat
    temp.3<- a+1-rho
    if (temp.3^2-4*a<0){
      #l.hat[k]<- x[k]
      l.hat[k]<-(l0.hat/2)*(temp.3)
      print(k)
      print('imaginary roots')
    }
    if (temp.3^2-4*a>=0){
      
      l.hat[k]<- (l0.hat/2)*(temp.3+sqrt(temp.3^2-4*a))
    }
  }
  return(list("l0.hat"=l0.hat,"l.hat"=l.hat))
}

q.est.agg<- function(rmt,A,X,b.tilde,tau,beta,eta,m,m.0,type){
  
  q<- dim(A)[1]
  f<- rep(1,q)
  if (type == 0){
    
    f<- f.agg.est(rmt,A,beta,tau,m)
    
  }
  G1<- Gfun.est(A,1,-1,beta,tau,rmt,m)
  G2<- Gfun.est(A,1,0,beta,tau,rmt,m)
  G3<- Gfun.est(A,0,1,0,tau,rmt,m)
 
  term1<- A%*%eta
  if(m>1){
    term2<- f*(G1%*%A%*%(colMeans(X)-eta))
  }
  if(m==1){
    term2<- f*(G1%*%A%*%(X-eta))
  }
  term3<- sqrt(m^{-1}*diag(G2)+m.0^{-1}*diag(G3))
  
  q<- term1+term2+term3*qnorm(b.tilde)
  
  return(list("q"=q,"f"=f,"G1"=G1,"G2"=G2,"G3"=G3))
  
}
q.Bayes.agg<- function(Sigma,A,X,b.tilde,tau,beta,eta,m,m.0){
  
  n<- dim(Sigma)[1]
  spd<- eigen(Sigma,symmetric = TRUE)
  V<- spd$vectors
  lam<- spd$values
  
  W<- chol2inv(chol(A%*%Sigma%*%t(A)))
  
  temp<- diag((lam)^(beta))#diag((lam)^(-beta))
  WW<- A%*%(V%*%temp%*%t(V))%*%t(A)
  
  C<- m*W+(1/tau)*(chol2inv(chol(WW)))#m*W+(1/tau)*(A%*%WW%*%t(A))
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- A%*%eta
  if(m>1){
    term2<- h1%*%A%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%A%*%(X-eta)
  }
  term3<- sqrt(diag(B)+m.0^(-1)*diag(A%*%Sigma%*%t(A)))
  
  q<- term1+term2+term3*qnorm(b.tilde)
  return(q)
  
}
q.naive.1.agg<- function(S,A,X,b.tilde,tau,beta,eta,m,m.0){
  
  n<- dim(S)[1]
  spd<- eigen(S)
  V<- spd$vectors
  lam<- spd$values
  
  W<- chol2inv(chol(A%*%S%*%t(A)))
  
  temp<- diag((lam)^(beta))#diag((lam)^(-beta))
  WW<- A%*%(V%*%temp%*%t(V))%*%t(A)
  
  C<- m*W+(1/tau)*(chol2inv(chol(WW)))#m*W+(1/tau)*(A%*%WW%*%t(A))
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- A%*%eta
  if(m>1){
    term2<- h1%*%A%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%A%*%(X-eta)
  }
  term3<- sqrt(diag(B)+m.0^(-1)*diag(A%*%S%*%t(A)))
  
  q<- term1+term2+term3*qnorm(b.tilde)
  return(q)
  
}
q.Bcv.agg<-function(S,A,X,b.tilde,tau,beta,eta,m,m.0){
  
  n<- dim(S)[1]
  spd<- eigen(S,symmetric = TRUE)
  V<- spd$vectors
  lam<- abs(spd$values)
  
  W<- chol2inv(chol(A%*%S%*%t(A)))
  
  temp<- diag((lam)^(beta))#diag((lam)^(-beta))
  WW<- A%*%(V%*%temp%*%t(V))%*%t(A)
  
  C<- m*W+(1/tau)*(chol2inv(chol(WW)))#m*W+(1/tau)*(A%*%WW%*%t(A))
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- A%*%eta
  if(m>1){
    term2<- h1%*%A%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%A%*%(X-eta)
  }
  term3<- sqrt(diag(B)+m.0^(-1)*diag(A%*%S%*%t(A)))
  
  q<- term1+term2+term3*qnorm(b.tilde)
  return(q)
  
}
q.factmle.agg<-function(S,A,X,b.tilde,tau,beta,eta,m,m.0){
  
  n<- dim(S)[1]
  spd<- eigen(S,symmetric = TRUE)
  V<- spd$vectors
  lam<- spd$values
  
  W<- chol2inv(chol(A%*%S%*%t(A)))
  
  temp<- diag((lam)^(beta))#diag((lam)^(-beta))
  WW<- A%*%(V%*%temp%*%t(V))%*%t(A)
  
  C<- m*W+(1/tau)*(chol2inv(chol(WW)))#m*W+(1/tau)*(A%*%WW%*%t(A))
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- A%*%eta
  if(m>1){
    term2<- h1%*%A%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%A%*%(X-eta)
  }
  term3<- sqrt(diag(B)+m.0^(-1)*diag(A%*%S%*%t(A)))
  
  q<- term1+term2+term3*qnorm(b.tilde)
  return(q)
  
}
q.poet.agg<-function(S,A,X,b.tilde,tau,beta,eta,m,m.0){
  
  n<- dim(S)[1]
  spd<- eigen(S,symmetric = TRUE)
  V<- spd$vectors
  lam<- abs(spd$values)
  
  W<- chol2inv(chol(A%*%S%*%t(A)))
  
  temp<- diag((lam)^(beta))#diag((lam)^(-beta))
  WW<- A%*%(V%*%temp%*%t(V))%*%t(A)
  
  C<- m*W+(1/tau)*(chol2inv(chol(WW)))#m*W+(1/tau)*(A%*%WW%*%t(A))
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- A%*%eta
  if(m>1){
    term2<- h1%*%A%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%A%*%(X-eta)
  }
  term3<- sqrt(diag(B)+m.0^(-1)*diag(A%*%S%*%t(A)))
  
  q<- term1+term2+term3*qnorm(b.tilde)
  return(q)
  
}

hfun.est<-function(r,alpha,beta,rmt,m,tau){
  
  pj<- rmt$pj
  l0.hat<- rmt$l0.hat
  l.hat<-rmt$l.hat
  zeta<-rmt$zeta
  K<-rmt$K
  
  n<- dim(pj)[1]
  h<-matrix(0,n,n)
  
  #if(beta<=1){
  
  a.1<- (l0.hat^(r+alpha))/((1+(1/(m*tau))*l0.hat^{1-beta})^r)
  
  for(k in 1:K){
    
    b.1<- (l.hat[k]^(r+alpha))/((1+(1/(m*tau))*l.hat[k]^{1-beta})^r)
    h<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+h
    
  }
  h<-h+a.1*diag(n)
  # }
  # if(beta>1){
  #   
  #   a.1<- (l0.hat^(r*beta+alpha))/((1/(m*tau)+l0.hat^(beta-1))^r)
  #   
  #   for(k in 1:K){
  #     
  #     b.1<- (l.hat[k]^(r*beta+alpha))/((1/(m*tau)+l.hat[k]^(beta-1))^r)
  #     h<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+h
  #     
  #   }
  #   h<-h+a.1*diag(n) 
  #   
  #   
  # }
  
  return(h) 
}
f.agg.est<-function(rmt,A,beta,tau,m){
  
  pj<- rmt$pj
  l0.hat<- rmt$l0.hat
  l.hat<-rmt$l.hat
  zeta<-rmt$zeta
  K<- rmt$K
  
  n<- dim(pj)[1]
  q<-dim(A)[1]
  
  G1<- Gfun.est(A,1,-1,beta,tau,rmt,m)
  g<-g.est(rmt,beta,tau,m)$G
  H<- G1%*%A%*%g%*%t(G1%*%A)
  
  numerator<- diag(H)
  
  g<- l0.hat+(m*tau)*l0.hat^beta
  M<-matrix(0,q,q)
  
  a.1<-1/(1+(1/(m*tau))*l0.hat^{1-beta})
  
  for(k in 1:K){
    
    b.1<- 1/(1+(1/(m*tau))*l.hat[k]^{1-beta})
    M<- (1/zeta[k]^4)*((b.1-a.1)^2)*(A%*%pj[,k]%*%t(A%*%pj[,k]))+M
  }
  
  R<- H+g*M
  denominator<- diag(R)
  return(numerator/denominator)
  
}
g.est<-function(rmt,beta,tau,m){
  
  pj<- rmt$pj
  l0.hat<- rmt$l0.hat
  l.hat<-rmt$l.hat
  zeta<-rmt$zeta
  K<- rmt$K
  
  n<- dim(pj)[1]
  G<-matrix(0,n,n)
  a.1<- l0.hat+(m*tau)*l0.hat^beta 
  for(k in 1:K){
    
    b.1<- l.hat[k]+(m*tau*l.hat[k]^beta)
    G<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+G
    
  }
  G<-G+a.1*diag(n)
  
  Ginv<-matrix(0,n,n)
  a.2<- 1/(l0.hat+(m*tau)*l0.hat^beta) 
  for(k in 1:K){
    
    b.2<- 1/(l.hat[k]+(m*tau*l.hat[k]^beta))
    Ginv<- (1/zeta[k]^2)*(b.2-a.2)*(pj[,k]%*%t(pj[,k]))+Ginv
    
  }
  Ginv<-Ginv+a.2*diag(n)
  return(list("G"=G,"Ginv"=Ginv))
  
}
g.agg.est<-function(rmt,A,beta,tau,m){
  
  pj<- rmt$pj
  l0.hat<- rmt$l0.hat
  l.hat<-rmt$l.hat
  zeta<-rmt$zeta
  K<- rmt$K
  
  n<- dim(pj)[1]
  G<-matrix(0,n,n)
  a.1<- l0.hat+(m*tau)*l0.hat^beta 
  for(k in 1:K){
    
    b.1<- l.hat[k]+(m*tau*l.hat[k]^beta)
    G<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+G
    
  }
  G<-A%*%(G+a.1*diag(n))%*%t(A)
  
  Ginv<-chol2inv(chol(G))
  return(list("G"=G,"Ginv"=Ginv))
  
}
Gfun.est<- function(A,r,alpha,beta,tau,rmt,m){
  
  H33<- A%*%hfun.est(0,1,0,rmt,m,tau)%*%t(A)
  
  if(alpha!=0){
    spd<- eigen(H33)
    eigvec<- spd$vectors
    temp<- diag((spd$values)^alpha)
    H3<- eigvec%*%temp%*%t(eigvec)
  }
  if(alpha==0){
    
    H3<- diag(dim(H33)[1])
  }
  H4<- diag(dim(H33)[1])
  if(r!=0){
    H1<- A%*%hfun.est(0,beta,0,rmt,m,tau)%*%t(A)
    H2<- A%*%(tau*hfun.est(0,beta,0,rmt,m,tau)+hfun.est(0,1,0,rmt,m,tau))%*%t(A)
    H2<- chol2inv(chol(H2))
    
    H4<- H1%*%H2%*%H33
    spd<- eigen(H4)
    eigvec<- spd$vectors
    
    temp<- diag((spd$values)^r)
    H4<- eigvec%*%temp%*%t(eigvec)
  }
  
  Gfun<- (tau^r)*(H4%*%H3)
  
  return(Gfun)
  
}
risk.agg.est<- function(q,theta,sigma.Y,b,h,type){
  
  if(type == 1){
  n<- length(sigma.Y)
  z<- (q-theta)/sigma.Y
  a.1<- (q-theta)*pnorm(z)
  a.2<- sigma.Y*dnorm(z)
  temp<- b*(theta-q)+(b+h)*(a.1+a.2)
  return(sum(temp))
  }
  if(type==0){
    
    z<- (q-theta)/sigma.Y
    temp = z^2
    return(sum(temp))
  }
  
}
taubeta.q.est<-function(grid.val,rmt,m,X){
  
  cv.taubeta<-matrix(0,nrow(grid.val),1)
  if (m>1){
    
    X<-colMeans(X)
  }
  for (t in 1:nrow(grid.val)){
    
    ttau<- grid.val[t,1]
    bbeta<-grid.val[t,2]
    
    out.g<-g.est(out.rmt,bbeta,ttau,m)
    S.q<- out.g$G
    Sinv.q<- out.g$Ginv
    cv.taubeta[t]<- -0.5*log(det(S.q))-0.5*t(X)%*%Sinv.q%*%X
    
  }
  o<- order(-cv.taubeta)
  return(colMeans(grid.val[o[1:3],]))

}
taubeta.est<-function(grid.val,S.val,m,X){
  
  cv.taubeta<-matrix(0,nrow(grid.val),1)
  if (m>1){
    
    X<-colMeans(X)
  }
  for (t in 1:nrow(grid.val)){
    
    ttau<- grid.val[t,1]
    bbeta<-grid.val[t,2]
    
    spd<-eigen(S.val)
    lam<-spd$values
    pp<-spd$vectors
    temp<- diag(lam^bbeta)
    
    S.1<- S.val+ttau*pp%*%temp%*%t(pp)
    Sinv.1<- chol2inv(chol(S.1))
    cv.taubeta[t]<- -0.5*log(det(S.1))-0.5*t(X)%*%Sinv.1%*%X
    
  }
  o<- order(-cv.taubeta)
  return(colMeans(grid.val[o[1:3],]))
  
}


#--------------------------------------
