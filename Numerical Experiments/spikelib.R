#' Detect spikes in a covariance model - KN test
#' 
#' details
#' This method calculates the number of spikes in a spiked covariance model
#' by iteratively comparing a standardized eigenvalue to the qunatile of Tracy-Wisdom
#' distribution for given significance level. The procedure is described in 
#' (Kritchman and Nadler, 2009). 
#' 
#' param lambdas eigenvalues to which the distribution is fitted
#' param numObs  number of data observations
#' param alpha significance level for the test
#' 
KN.test <- function(S, numObs, alpha = 0.05) {
  
  lambdas<- eigen(S)$values
  p <- length(lambdas); max.spikes <- min(numObs, p) - 1; sigma.sq <- 0
  s_alpha <- (-3/2 * log(4*sqrt(pi) * alpha/100 ))^(2/3)
  
  for(spikes in 1:max.spikes) {
    
    mu_np <- 1/numObs*(sqrt(numObs-1/2) + sqrt((p - spikes)-1/2))^2
    sigma_np <- sqrt(mu_np/numObs) * (1/sqrt(numObs-1/2) + 1 / sqrt((p-spikes)-1/2) )^(1/3)
    
    sigma.sq.temp <- sigma.sq.est(lambdas,numObs,spikes) 
    if(lambdas[spikes] <= sigma.sq.temp*(mu_np + s_alpha * sigma_np)) break
    sigma.sq <- sigma.sq.temp
    #print(spikes)
  }
  
  est.sigma.sq <- if(spikes > 1) sigma.sq else  sum(lambdas)/p
  return(list(numOfSpikes = spikes - 1, sigma.sq = est.sigma.sq))
}

#' Detect spikes in a covariance model - PY
#' 
PY <- function(lambdas, numObs,C) {
  
  p <- length(lambdas); max.spikes <- min(numObs, p) - 1; 
  spikes<-0
  sigma.sq <- sigma.sq.est.PY(lambdas,spikes)
  dn<- 1/C
  delta<- lambdas[1:(p-1)]-lambdas[2:p]
  
  for(i in 1:max.spikes) {
    
    if (delta[i+1]<dn & delta[i+2]<dn){
      
      spikes<- i
      break
    }
    spikes<- i
    sigma.sq <- sigma.sq.est.PY(lambdas,spikes)
    
  }
  
  
  est.sigma.sq<- sigma.sq.est.PY(lambdas,spikes)
  return(list(numOfSpikes = spikes, sigma.sq = est.sigma.sq))
}



#' Calculate unknown noise estimate
#' 
#' details
#' This method calculates the unknown unbiased estimate of noise for the eigenvalues
#' based on the procedure described in (Kritchman and Nadler, 2009). It involves
#' solving a couple of nonlinear equations.
#' 
#' param lambdas eigenvalues to which the distribution is fitted
#' param numObs  number of data observations
#' param numOfSpikes number of spikes in the spike covariance model
#' 
sigma.sq.est <- function(lambdas , numObs, numOfSpikes) {
  
  p <- length(lambdas)
  sigma.sq <- 1/(p-numOfSpikes) * sum(lambdas[(numOfSpikes+1):p])
  
  while(TRUE) {
    
    tmp.b <- lambdas[1:numOfSpikes] + sigma.sq - sigma.sq*(p-numOfSpikes)/numObs 
    discriminant <- tmp.b^2 - 4 * lambdas[1:numOfSpikes]*sigma.sq
    if (any( discriminant < 0)) break
    
    rho <- 0.5*(tmp.b + sqrt(discriminant)) 
    sigma.sq.new <- 1/(p-numOfSpikes) * (sum(lambdas) - sum(rho))
    if (abs(sigma.sq.new - sigma.sq)/sigma.sq < 1e-8) break
    sigma.sq <- max(sigma.sq.new,1e-6)
  }
  
  return(sigma.sq)
}

#' Calculate unknown noise estimate from 
#' D.Passemier,J.-f.Yao,On determining the number of spikes in a 
#' high-dimensional spiked population model
#' 
sigma.sq.est.PY <- function(lambdas,numOfSpikes) {
  
  p <- length(lambdas)
  sigma.sq.est.PY<- NA
  if(numOfSpikes<p){
    sigma.sq.est.PY <- 1/(p-numOfSpikes) * sum(lambdas[(numOfSpikes+1):p])
  }
 return(sigma.sq.est.PY)
}