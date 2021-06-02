
# Experiment 1 - scenario 2

source('linexloss.R')
source('spikelib.R')

require(ggplot2)
require(gridExtra)
require(Matrix)
library(foreach)
library(doParallel)

gen.spikedcov<-function(n,mz,m.0,K,l0,tau,beta){
  
  #1. Generate the n x n spiked covariance matrix 
  spikes<-seq(25,5,length.out = K)

  set.seed(1)
  A<- matrix(runif(n^2,1,2),n,n)
  p<- orthonormalization(A, basis=TRUE, norm=TRUE)
  Sigma<- matrix(0,n,n)
  temp<- Sigma
  for(k in 1:K){

    temp<- p[,k]%*%t(p[,k])+temp
    Sigma<- spikes[k]*(p[,k]%*%t(p[,k]))+Sigma

  }
  Sigma<- Sigma+l0*(diag(n)-temp)
  # set.seed(1)
  # B<- matrix(rnorm(n*K,0,0.5),n,K)
  # Sigma<-B%*%t(B)+l0*diag(n)
  
  
  spd<-eigen(Sigma)
  p<-spd$vectors
  temp<- diag((spd$values)^beta)
  Sigma.theta<- tau*(p%*%temp%*%t(p))
  SigmaY<- (1/m.0)*Sigma
  
  # set.seed(1)
  # Z<- mvrnorm(mz,rep(0,n),Sigma)
  # S<- cov(Z)*(mz/(mz-1))
  # 
  # return(list("Sigma"=Sigma,"S"=S,"Sigma.Y"=SigmaY,
  #             "Sigma.theta"=Sigma.theta,"W"=Z))
  
  return(list("Sigma"=Sigma,"Sigma.Y"=SigmaY,
              "Sigma.theta"=Sigma.theta))
  
}
gen.data<-function(eta,m,Sigma,Sigma.theta,reps){
  
  set.seed(reps)
  theta<- mvrnorm(1,eta,Sigma.theta)
  set.seed(reps^2)
  X<- mvrnorm(m,theta,Sigma)
  return(list("theta"=theta,"X"=X))
}


#--- Fixed Parameters ---------------------------------  

n<- 200 # dimensionality
K.true<- 10
l0<-1
set.seed(1)
a<- runif(n,-2,-1)#check loss parameters
b<- rep(1,n)#check loss parameters
m<-1 # sample sizes
Mz<- seq(15,50,5)
m.0<- 1 # sample size for future data Y
reps<- 500
tau<-1
beta<-1.75
eta<-rep(0,n)
tau.grid<-c(0.6,0.7,0.8,0.9,1)
beta.grid<-c(1.5,1.6,1.7,1.8,1.9,2)


#----Initialization ------------------------
taubeta.grid<- cbind(rep(tau.grid,each = length(beta.grid)),
                     rep(beta.grid,length(beta.grid)))

taubeta.estimated<-array(rep(0,length(Mz)*5*2),
                         c(length(Mz),5,2))

taubeta.std<-taubeta.estimated

factors.estimated<-matrix(0,length(Mz),5)
factors.std<-factors.estimated
shrink.f<-matrix(0,n,reps)
ineff.q<-matrix(0,length(Mz),1)
ineff.naive<-ineff.q
ineff.Bcv<-ineff.q
ineff.poet<-ineff.q
ineff.fact<-ineff.q
risk.avg<-array(rep(0,length(Mz)*reps*7),
                c(length(Mz),reps,7))

#--- Computation  ------------------------------------
for(i in 1:length(Mz)){
  
  mz<-Mz[i]
  out.Sigma<-gen.spikedcov(n,mz,m.0,K.true,l0,tau,beta)
  sigma.Y<- sqrt(diag(out.Sigma$Sigma.Y))
  
  out.rmtBayes<- rmt.est(K.true,out.Sigma$Sigma,mz,0)
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  result<-foreach(N = 1:reps,.export=c("mvrnorm","EsaBcv","POET","Factmle_cov"))%dopar%{
    
    set.seed(N)
    W<- mvrnorm(mz,rep(0,n),out.Sigma$Sigma)
    S<- cov(W)*(mz/(mz-1))
    
    K.est<- KN.test(S, mz, 0.05)
    out.rmt<- rmt.est(K.est$numOfSpikes,S,mz,1)
    
    # Naive
    S.naive<- S.est.naive.1(S,K.est)
    
    # Bcv
    out.Bcv<- EsaBcv(W)
    K.Bcv<- out.Bcv$best.r
    if(K.Bcv>0){
      
      S.Bcv<- cov(out.Bcv$estS)+out.Bcv$estSigma*diag(n)
      
    }
    if(K.Bcv==0){
      
      S.Bcv<- out.Bcv$estSigma*diag(n)
      
    }
    #POET
    if(i==1){
      out.poet<- POET(t(W),K.est$numOfSpikes,matrix='vad')
      S.poet<- out.poet$SigmaY
      K.poet<- dim(out.poet$factors)[1]
    }
    if(i>1){
      S.poet<- S.naive
      K.poet<-K.est$numOfSpikes
    }
    #FACTMLE
    out.fact<- Factmle_cov(S,K.est$numOfSpikes)
    S.fact<- diag(out.fact$Psi)+out.fact$Lambda%*%t(out.fact$Lambda)
    nfactors<-c(K.true, K.est$numOfSpikes,K.Bcv,K.poet,K.est$numOfSpikes)
    
    simdata<- gen.data(eta,m,out.Sigma$Sigma,out.Sigma$Sigma.theta,N)
    theta<- simdata$theta
    X<- simdata$X
    
    taubeta.q<- taubeta.q.est(taubeta.grid,out.rmt,m,X)
    taubeta.naive<- taubeta.est(taubeta.grid,S.naive,m,X)
    taubeta.Bcv<- taubeta.est(taubeta.grid,S.Bcv,m,X)
    taubeta.poet<- taubeta.est(taubeta.grid,S.poet,m,X)
    taubeta.fact<- taubeta.est(taubeta.grid,S.fact,m,X)
    
    out.q<- q.est(out.rmt,X,a,taubeta.q[1],
                  taubeta.q[2],eta,m,m.0,0)
    qhat<-out.q$q
    f<-out.q$f
    
    qhat.1<- q.est(out.rmt,X,a,taubeta.q[1],
                   taubeta.q[2],eta,m,m.0,1)$q
    qhat.Bayes<- q.Bayes(out.Sigma$Sigma,X,a,tau,beta,eta,m,m.0)
    qhat.naive<- q.naive.1(S.naive,X,a,
                           taubeta.naive[1],taubeta.naive[2],eta,m,m.0)
    qhat.Bcv<- q.Bcv(S.Bcv,X,a,
                     taubeta.Bcv[1],taubeta.Bcv[2],eta,m,m.0)
    qhat.poet<- q.poet(S.poet,X,a,
                       taubeta.poet[1],taubeta.poet[2],eta,m,m.0)
    qhat.fact<- q.factmle(S.fact,X,a,
                          taubeta.fact[1],taubeta.fact[2],eta,m,m.0)
    
    risk.q<- risk.est(qhat,theta,sigma.Y,a,b)
    risk.1<- risk.est(qhat.1,theta,sigma.Y,a,b)
    risk.Bayes<- risk.est(qhat.Bayes,theta,sigma.Y,a,b)
    risk.naive<- risk.est(qhat.naive,theta,sigma.Y,a,b)
    risk.Bcv<- risk.est(qhat.Bcv,theta,sigma.Y,a,b)
    risk.poet<- risk.est(qhat.poet,theta,sigma.Y,a,b)
    risk.fact<- risk.est(qhat.fact,theta,sigma.Y,a,b)
    return(list("q"=risk.q,"q.1"=risk.1,"bayes"=risk.Bayes,
                "naive"=risk.naive,"bcv"=risk.Bcv,
                "poet"=risk.poet,"fact"=risk.fact,"f"=f,
                "nfactors"=nfactors,
                "taubeta"=rbind(taubeta.q,taubeta.naive,taubeta.Bcv,
                                taubeta.poet,taubeta.fact)))
  }
  stopCluster(cl)
  registerDoSEQ()
  temprisk.q<- sapply(1:reps,function(i) result[[i]]$q)
  temprisk.1<- sapply(1:reps,function(i) result[[i]]$q.1)
  temprisk.Bayes<- sapply(1:reps,function(i) result[[i]]$bayes)
  temprisk.naive<- sapply(1:reps,function(i) result[[i]]$naive)
  temprisk.Bcv<- sapply(1:reps,function(i) result[[i]]$bcv)
  temprisk.poet<- sapply(1:reps,function(i) result[[i]]$poet)
  temprisk.fact<- sapply(1:reps,function(i) result[[i]]$fact)
 
  taubeta.estimated[i,,]<- rbind(
    rowMeans(sapply(1:reps,function(i) result[[i]]$taubeta[1,])),
    rowMeans(sapply(1:reps,function(i) result[[i]]$taubeta[2,])),
    rowMeans(sapply(1:reps,function(i) result[[i]]$taubeta[3,])),
    rowMeans(sapply(1:reps,function(i) result[[i]]$taubeta[4,])),
    rowMeans(sapply(1:reps,function(i) result[[i]]$taubeta[5,])))
  
  taubeta.std[i,,]<- rbind(
    apply(sapply(1:reps,function(j) result[[j]]$taubeta[1,]),1,sd),
    apply(sapply(1:reps,function(j) result[[j]]$taubeta[2,]),1,sd),
    apply(sapply(1:reps,function(j) result[[j]]$taubeta[3,]),1,sd),
    apply(sapply(1:reps,function(j) result[[j]]$taubeta[4,]),1,sd),
    apply(sapply(1:reps,function(j) result[[j]]$taubeta[5,]),1,sd))/sqrt(reps)
  
  factors.estimated[i,]<- rowMeans(sapply(1:reps,function(j) result[[j]]$nfactors))
  factors.std[i,]<- apply(sapply(1:reps,function(j) result[[j]]$nfactors),1,sd)/sqrt(reps)
  
  if(i == 1){
    result.i1<-result
    shrink.f<- cbind(sapply(1:reps,function(j) result[[j]]$f))
  }
  
  risk.avg[i,,]<-c(temprisk.Bayes,temprisk.q,temprisk.1,
                   temprisk.naive,temprisk.Bcv,
                   temprisk.poet,temprisk.fact)
  tt<-colMeans(risk.avg[i,,])
  risk.diff<-tt[3]-tt[1]
  ineff.q[i]<- (tt[2]-tt[1])/risk.diff
  ineff.naive[i]<- (tt[4]-tt[1])/risk.diff
  ineff.Bcv[i]<- (tt[5]-tt[1])/risk.diff
  ineff.poet[i]<- (tt[6]-tt[1])/risk.diff
  ineff.fact[i]<- (tt[7]-tt[1])/risk.diff
  print(i)
}

 ineff.mat<- c(ineff.q,ineff.naive,ineff.Bcv,ineff.fact)
 name<- rep(c('CASP','Naive','Bcv','FactMLE'),each=8)

 shrink<- as.data.frame(rowMeans(shrink.f))
 names(shrink)<- "f"
 shrink$lower<- sapply(1:n,function(i) quantile(shrink.f[i,],0.10))
 shrink$upper<- sapply(1:n,function(i) quantile(shrink.f[i,],0.90))
 sorted.f<- order(shrink$f)
 shrink$f<-shrink$f[sorted.f]
 shrink$lower<-shrink$lower[sorted.f]
 shrink$upper<-shrink$upper[sorted.f]
 shrink$dim<- 1:n

 plotdata1<- as.data.frame(ineff.mat)
 plotdata1$name<- factor(name)
 plotdata1$x<- rep(Mz,4)
 plotdata2<- as.data.frame(rep(1,8))
 plotdata2$x<- Mz
 names(plotdata2)<- c("Y","x")
 g1<-ggplot()+geom_line(data=plotdata2,aes(x = x, y = Y),size=1.5)+
   geom_line(data=plotdata1, aes(x = x, y = ineff.mat,color=name),size=1.5)+
   geom_point(data=plotdata1,aes(x = x, y = ineff.mat,color=name,
                                 shape=name),size=1.5)+
   scale_x_continuous(breaks = Mz)+
   scale_y_continuous(breaks = seq(0.95,5,by=0.5))+
   xlab(expression(m))+ylab("REE")+theme_bw()+
   theme(legend.position='top',legend.title=element_blank(),
         legend.background = element_rect(fill="white",
                                          size=1, linetype="solid", colour ="black"),
         legend.text=element_text(size=20),
         axis.text.x = element_text(size=20),
         axis.text.y = element_text(size=20),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size=20))

 g2<- ggplot()+geom_point(data=shrink,aes(x=dim,y=f),color="red")+
   geom_line(data=shrink,aes(x=dim,y=f),color="red",size=1.5)+
   geom_ribbon(data=shrink,aes(x=dim,ymin=lower,ymax=upper),alpha=0.5,fill="gray")+
   ylab("Shrinkage Factors")+xlab("Dimensions")+
   theme_bw()+
   theme(axis.text.x = element_text(size=20),
         axis.text.y = element_text(size=20),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size=20))

 grid.arrange(g1,g2,ncol=2)
 rm('result')
 save.image(paste(getwd(),'/exp12_linexloss.RData',sep=''))
 
 g <- arrangeGrob(g1,g2, ncol=2) #generates g
 ggsave(file="exp12linexloss.pdf", g,width = 35, height=15, units = "cm",dpi=500,scale=1)

