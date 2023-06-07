library(rjags)
library(MASS)
source("CLNMA_functions.R")


set.seed(1)
mu1 = c(0.5,1)
mu2 = c(0,-0.5)
tau1 =c(0.25,0.36)
tau2 =c(0.36,0.16)

tau = c(tau1,tau2)
rho = matrix(c(1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1),
             nrow = 4)
tau_BC = tau1+tau2-2*rho[1,2]*sqrt(tau1*tau2)


betweenv = diag(tau)
for(i in 1:3){
  for(j in (i+1):4){
    betweenv[i,j] =betweenv[j,i]= rho[i,j]*sqrt(tau[i]*tau[j])
  }
}
ss1 =1
ss2= 1
rho_w = matrix(c(1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1),
               nrow = 4)
nab = nac=nbc=nabc= 5

start = Sys.time()
# simulate data
dataout = gendata(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2)
# transfer data from long format to wide format
N_total = dim(dataout)[1]
BA_outcome1 = numeric(dim(dataout)[1])
CA_outcome1 = numeric(dim(dataout)[1])
BC_outcome1 = numeric(dim(dataout)[1])
BA_outcome2 = numeric(dim(dataout)[1])
CA_outcome2 = numeric(dim(dataout)[1])
BC_outcome2 = numeric(dim(dataout)[1])
BA_sd1 = numeric(dim(dataout)[1])
CA_sd1 = numeric(dim(dataout)[1])
BC_sd1 = numeric(dim(dataout)[1])
BA_sd2 = numeric(dim(dataout)[1])
CA_sd2 = numeric(dim(dataout)[1])
BC_sd2 = numeric(dim(dataout)[1])
for(i in 1:N_total){
  temp_record = dataout[i,]
  temp_t1 = dataout[i,"t1"]
  temp_t2 = dataout[i,"t1"]
  if(temp_t1 == 1){
    if(temp_t2 == 2){ 
      BC_outcome1[i] = dataout[i,"outcome1"]
      BC_sd1[i] = dataout[i,"sd1"]
      BC_outcome2[i] = dataout[i,"outcome2"]
      BC_sd2[i] = dataout[i,"sd2"]
    } else {
      BA_outcome1[i] = dataout[i,"outcome1"]
      BA_sd1[i] = dataout[i,"sd1"]
      BA_outcome2[i] = dataout[i,"outcome2"]
      BA_sd2[i] = dataout[i,"sd2"]
    } 
  } else {
    CA_outcome1[i] = dataout[i,"outcome1"]
    CA_sd1 = dataout[i,"sd1"]
    CA_outcome2[i] = dataout[i,"outcome2"]
    CA_sd2 = dataout[i,"sd2"]
  }
}

replace_zero_with_mean <- function(x){
  x[x==0] = NA
  x[is.na(x)] =mean(x, na.rm=T)
  return(x)
}

BA_outcome1 = replace_zero_with_mean(BA_outcome1)
CA_outcome1 = replace_zero_with_mean(CA_outcome1)
BC_outcome1 = replace_zero_with_mean(BC_outcome1)
BA_outcome2 = replace_zero_with_mean(BA_outcome2)
CA_outcome2 = replace_zero_with_mean(CA_outcome2)
BC_outcome2 = replace_zero_with_mean(CA_outcome2)
BA_sd1 = replace_zero_with_mean(BA_sd1)
CA_sd1 = replace_zero_with_mean(CA_sd1)
BC_sd1 = replace_zero_with_mean(BC_sd1)
BA_sd2 = replace_zero_with_mean(BA_sd2)
CA_sd2 = replace_zero_with_mean(CA_sd2)
BC_sd2 = replace_zero_with_mean(BC_sd2)


# Reference :
# http://doingbayesiandataanalysis.blogspot.com/2017/06/bayesian-estimation-of-correlations-and.html

y = cbind(BA_outcome1, CA_outcome1, BC_outcome1,BA_outcome2, CA_outcome2, BC_outcome2)
y = rmvnorm(30, mean =rep(0,6), sigma = diag(6))
# Assemble data for sending to JAGS:
datalist = list(
  y = y,
  ss = cbind(BA_sd1^2, CA_sd1^2, BC_sd1^2,
            BA_sd2^2, CA_sd2^2, BC_sd2^2),
  Ntotal = nrow(y),
  Nvar = ncol(y),
  rho_w = diag(x=1,nrow=ncol(y)),
  rho_b = diag(x=1,nrow=ncol(y))
)


# Define the model
modelString = "
model {
  # Likelihood
  for ( i in 1:Ntotal ) {
    
    y[i,1:Nvar] ~ dmnorm(mu[1:Nvar] ,InvCovMat[i,1:Nvar,1:Nvar] )
    
    #InvCovMat[i,1:Nvar,1:Nvar] = Inv_Between_Study_CovMat
    #CovMat = Between_Study_CovMat
    
    #InvCovMat[i,1:Nvar,1:Nvar] = inverse(CovMat[i,1:Nvar,1:Nvar])
    CovMat[i,1:Nvar,1:Nvar] = With_Study_CovMat[i,1:Nvar,1:Nvar] + Between_Study_CovMat[1:Nvar,1:Nvar]
    InvCovMat[i,1:Nvar,1:Nvar] = inverse(CovMat[i,1:Nvar,1:Nvar])
    With_Study_CovMat[i,1:Nvar,1:Nvar] = diag_ss[i, 1:Nvar, 1:Nvar] %*% rho_w %*% diag_ss[i, 1:Nvar, 1:Nvar]
    diag_ss[i,  1:Nvar, 1:Nvar] =  temp_ss[i, 1:Nvar, 1:Nvar]
    #diag_tau[i, 1:Nvar, 1:Nvar] =  temp_tau[i, 1:Nvar, 1:Nvar]
    
    
    for(j in 1:Nvar){
        for(k in 1:Nvar){
           temp_ss[i,j,k]  = ifelse(j == k, ss[i,j],0)
        }
    }
    
   #  for(j in 1:Nvar){
  #      for(k in 1:Nvar){
  #         temp_tau[i,j,k]  = ifelse(j == k,exp(2*tau[j]),0)
  #      }
  #  }
    
  }

  
  # Priors
  # mu
  for ( varIdx in 1:Nvar ) { mu[varIdx] ~ dnorm(0,exp(2*tau[varIdx]) ) }
  
  # within study correlation
  for ( Idx1 in 2:Nvar ) {
    for ( Idx2 in 1:(Idx1-1)) {
    rho_w[Idx1,Idx2] ~ dunif(0,1)
    rho_w[Idx2,Idx1] ~ dunif(0,1)
    }
  }
  
  # between study correlation
  for ( Idx1 in 2:Nvar ) {
    rho_b[Idx1,Idx1]  ~ dunif(0,1)
    for ( Idx2 in 1:(Idx1-1)) {
      rho_b[Idx1,Idx2] ~ dunif(0,1)
      rho_b[Idx2,Idx1] ~ dunif(0,1)
    }
  }
  


  # between study covariance
  #Inv_Between_Study_CovMat ~ dwish(rho_b, Nvar)
  #Between_Study_CovMat = inverse(Inv_Between_Study_CovMat)
  
   for ( Idx1 in 1:Nvar ) {
      tau[Idx1] ~ dunif(0,1)
    }
  
  for(j in 1:Nvar){
      for(k in 1:Nvar){
           temp_tau[j,k]  = ifelse(j == k, exp(-2*tau[j]), 0)
        }
    }
  
  Between_Study_CovMat[1:Nvar,1:Nvar] = temp_tau[1:Nvar, 1:Nvar] %*% rho_b %*% temp_tau[1:Nvar, 1:Nvar]
  
  for ( varIdx in 1:Nvar ) {
    sigma2[varIdx] =Between_Study_CovMat[varIdx, varIdx]
  }  
}


" # close quote for modelString
writeLines( modelString , con="Jags-MultivariateNormal-model.txt" )

# Run the chains:
nChain = 3
nAdapt = 500
nBurnIn = 500
nThin = 10
nStepToSave = 20000

jagsModel = jags.model( file="Jags-MultivariateNormal-model.txt" ,
                        data=datalist , n.chains=nChain , n.adapt=nAdapt )
update( jagsModel , n.iter=nBurnIn )
codaSamples = coda.samples( jagsModel ,
                            variable.names=c("mu","sigma2") ,
                            n.iter=nStepToSave/nChain*nThin , thin=nThin )
time_end = Sys.time()
difference <- difftime(timeEnd, timeStart, units='mins')
print(paste("Running Time is ", round(difference,4)))
