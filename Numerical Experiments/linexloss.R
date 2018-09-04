# Holds functions for Linex loss
require(MASS)
require(far)
require(esaBcv)
require(FACTMLE)
require(POET)
require(dpglasso)

rmt.est<- function(K,S,m,type){
  
  # Take S and adjust eigen values and eigen vectors. To be used in q.est
  
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
  
  return(list("l0.hat"=l0.hat,"l.hat"=l.hat,"zeta"=zeta,
              "pj"=spd$vectors,"K"=K))
  
}
S.est.naive.1<-function(S,K){
  
  n<- dim(S)[1]
  spd<- eigen(S)
  l.tilde<-spd$values
  l0.tilde<- ((n-K$numOfSpikes)^{-1})*sum(l.tilde[(K$numOfSpikes+1):n])#K$sigma.sq#((n-K)^{-1})*sum(l.tilde[(K+1):n])
  pj<- spd$vectors
  
  S.naive<- matrix(0,n,n)
  temp<- S.naive
  for(k in 1:K$numOfSpikes){
    
    temp<- pj[,k]%*%t(pj[,k])+temp
    S.naive<- l.tilde[k]*(pj[,k]%*%t(pj[,k]))+S.naive
    
  }
  S.naive<- S.naive+l0.tilde*(diag(n)-temp)
  return(S.naive)
}

eigen.est<-function(l,l.0,rho,K){
  
  x<- matrix(0,K,1)
  b<-l.0
  
  for(k in 1:K){
    a<- l[k]
    temp.1<-b+a-rho*b
    x[k]<- (temp.1+sqrt(temp.1^2-4*a*b))/(2)
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
      print('imaginary roots')
    }
    if (temp.3^2-4*a>=0){
      
      l.hat[k]<- (l0.hat/2)*(temp.3+sqrt(temp.3^2-4*a))
    }
  }
  return(list("l0.hat"=l0.hat,"l.hat"=l.hat))
}

q.Bayes<- function(Sigma,X,a,tau,beta,eta,m,m.0){
  
  n<- dim(Sigma)[1]
  spd<- eigen(Sigma,symmetric = TRUE)
  V<- spd$vectors
  lam<- spd$values
  
  W<- chol2inv(chol(Sigma))
  
  temp<- diag((lam)^(-beta))
  WW<- V%*%temp%*%t(V)
  
  C<- m*W+(1/tau)*WW
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- eta
  if(m>1){
    term2<- h1%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%(X-eta)
  }
  term3<- diag(B)+m.0^(-1)*diag(Sigma)
  
  q<- term1+term2-term3*(a/2)
  return(q)
  
}
q.naive.1<- function(S.naive,X,a,tau,beta,eta,m,m.0){
  
  n<- dim(S.naive)[1]
  spd<- eigen(S.naive,symmetric = TRUE)
  V<- spd$vectors
  lam<- spd$values
  
  W<- chol2inv(chol(S.naive))
  
  temp<- diag((lam)^(-beta))
  WW<- V%*%temp%*%t(V)
  
  C<- m*W+(1/tau)*WW
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- eta
  if(m>1){
    term2<- h1%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%(X-eta)
  }
  term3<- diag(B)+m.0^(-1)*diag(S.naive)
  
  q<- term1+term2-term3*(a/2)
  return(q)
  
}
q.est<- function(rmt,X,a,tau,beta,eta,m,m.0,type){
  
  #Use rmt.est
  
  f<- rep(1,n)
  if (type == 0){
    
    f<- f.est(rmt,beta,tau,m)
    
  }
  
  h1<- hfun.est(1,-1,beta,rmt,m,tau)
  h2<- hfun.est(1,0,beta,rmt,m,tau)
  h3<- hfun.est(0,1,0,rmt,m,tau)
  term1<- eta
  if(m>1){
    term2<- f*(h1%*%(colMeans(X)-eta))
  }
  if(m==1){
    term2<- f*(h1%*%(X-eta))
  }
  term3<- m^{-1}*diag(h2)+m.0^{-1}*diag(h3)
  
  q<- term1+term2-term3*(a/2)
  return(list("q"=q,"f"=f,"h1"=h1,"h2"=h2,"h3"=h3))
  
}
q.Bcv<-function(S.Bcv,X,a,tau,beta,eta,m,m.0){
  
  n<- dim(S.Bcv)[1]
  spd<- eigen(S.Bcv,symmetric = TRUE)
  V<- spd$vectors
  lam<- spd$values
  
  W<- chol2inv(chol(S.Bcv))
  
  temp<- diag((lam)^(-beta))
  WW<- V%*%temp%*%t(V)
  
  C<- m*W+(1/tau)*WW
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- eta
  if(m>1){
    term2<- h1%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%(X-eta)
  }
  term3<- diag(B)+m.0^(-1)*diag(S.Bcv)
  
  q<- term1+term2-term3*(a/2)
  return(q)
  
}
q.factmle<-function(S,X,a,tau,beta,eta,m,m.0){
  
  n<- dim(S)[1]
  spd<- eigen(S,symmetric = TRUE)
  V<- spd$vectors
  lam<- spd$values
  
  W<- chol2inv(chol(S))
  
  temp<- diag((lam)^(-beta))
  WW<- V%*%temp%*%t(V)
  
  C<- m*W+(1/tau)*WW
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- eta
  if(m>1){
    term2<- h1%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%(X-eta)
  }
  term3<- diag(B)+m.0^(-1)*diag(S)
  
  q<- term1+term2-term3*(a/2)
  return(q)
  
}
q.poet<-function(S,X,a,tau,beta,eta,m,m.0){
  
  n<- dim(S)[1]
  spd<- eigen(S,symmetric = TRUE)
  V<- spd$vectors
  lam<- spd$values
  
  W<- chol2inv(chol(S))
  
  temp<- diag((lam)^(-beta))
  WW<- V%*%temp%*%t(V)
  
  C<- m*W+(1/tau)*WW
  B<- chol2inv(chol(C))
  
  h1<- m*(B%*%W)
  term1<- eta
  if(m>1){
    term2<- h1%*%(colMeans(X)-eta)
  }
  if(m==1){
    term2<- h1%*%(X-eta)
  }
  term3<- diag(B)+m.0^(-1)*diag(S)
  
  q<- term1+term2-term3*(a/2)
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
  
  a.1<- (l0.hat^(r+alpha))/((1+(1/(m*tau))*l0.hat^{1-beta})^r)
  
  for(k in 1:K){
    
    b.1<- (l.hat[k]^(r+alpha))/((1+(1/(m*tau))*l.hat[k]^{1-beta})^r)
    h<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+h
    
  }
  h<-h+a.1*diag(n)
  
  return(h) 
}
f.est<-function(rmt,beta,tau,m){
  
  pj<- rmt$pj
  l0.hat<- rmt$l0.hat
  l.hat<-rmt$l.hat
  zeta<-rmt$zeta
  K<- rmt$K
  
  n<- dim(pj)[1]
  H<-matrix(0,n,n)
  
  a.1<- ((m*tau)*l0.hat^beta)/(1+(1/(m*tau))*l0.hat^{1-beta}) 
  
  for(k in 1:K){
    
    b.1<- (m*tau*l.hat[k]^beta)/(1+(1/(m*tau))*l.hat[k]^{1-beta})
    H<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+H
    
  }
  H<-H+a.1*diag(n)
  numerator<- diag(H)
  
  g<- l0.hat+(m*tau)*l0.hat^beta
  M<-matrix(0,n,n)
  
  a.1<-1/(1+(1/(m*tau))*l0.hat^{1-beta})
  
  for(k in 1:K){
    
    b.1<- 1/(1+(1/(m*tau))*l.hat[k]^{1-beta})
    M<- (1/zeta[k]^4)*((b.1-a.1)^2)*(pj[,k]%*%t(pj[,k]))+M
  }
  
  R<- H+g*M
  #denominator<- sqrt(diag(R))
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

risk.est<- function(q,theta,sigma.Y,a,b){
  
  z.1<- a*(q-theta)
  z.2<- exp(z.1+0.5*(a^2)*sigma.Y^2)
  temp<- b*(z.2-z.1-1)
  return(sum(temp))
  
}
taubeta.q.est<-function(grid.val,rmt,m,X){
  
  cv.taubeta<-matrix(0,nrow(grid.val),1)
  if (m>1){
    
    X<-colMeans(X)
  }
  for (t in 1:nrow(grid.val)){
    
    ttau<- grid.val[t,1]
    bbeta<-grid.val[t,2]
    
    S.q<- g.est(out.rmt,bbeta,ttau,m)$G
    Sinv.q<- g.est(out.rmt,bbeta,ttau,m)$Ginv
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
    
    S.1<- S.val+ttau*(pp%*%temp%*%t(pp))
    Sinv.1<- chol2inv(chol(S.1))
    cv.taubeta[t]<- -0.5*log(det(S.1))-0.5*t(X)%*%Sinv.1%*%X
    
  }
  o<- order(-cv.taubeta)
  return(colMeans(grid.val[o[1:3],]))
  
}