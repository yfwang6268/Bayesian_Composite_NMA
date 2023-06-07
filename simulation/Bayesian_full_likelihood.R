library(rjags)
library(MASS)
source("CLNMA_functions.R")


set.seed(1)
mu1 = c(0.5,1)
mu2 = c(0,-0.5)
tau1 =c(0.25,0.36)
tau2 =c(0.36,0.16)

tau = c(tau1,tau2)
rho = matrix(c(1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,1),
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
rho_w = matrix(c(1,0.2,0.2,0.2,0.2,1,0.2,0.2,0.2,0.2,1,0.2,0.2,0.2,0.2,1),
               nrow = 4)
nab = nac=nbc=nabc= 5

simulation_time = 2

posterior_mu = NULL
posterior_sd = NULL
  

for(sim_t in 1: simulation_time){
  start = Sys.time()
  # simulate data
  dataout = gendata(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2)
  # transfer data from long format to wide format
  N_total = nab+nac+nbc+nabc
  BA_outcome1 = numeric(N_total)
  CA_outcome1 = numeric(N_total)
  BC_outcome1 = numeric(N_total)
  BA_outcome2 = numeric(N_total)
  CA_outcome2 = numeric(N_total)
  BC_outcome2 = numeric(N_total)
  BA_sd1 = numeric(N_total)
  CA_sd1 = numeric(N_total)
  BC_sd1 = numeric(N_total)
  BA_sd2 = numeric(N_total)
  CA_sd2 = numeric(N_total)
  BC_sd2 = numeric(N_total)
  Two_Arm_VarIdx = numeric(N_total) # for two arm study
  
  for(i in 1:N_total){
    temp_record = dataout[i,]
    if(i <= nab+nac+nbc){
      # Two arm study
      temp_t1 = dataout[i,"t1"]
      temp_t2 = dataout[i,"t2"]
      
      if(temp_t1 == 1){
        if(temp_t2 == 2){ 
          BC_outcome1[i] = dataout[i,"outcome1"]
          BC_sd1[i] = dataout[i,"sd1"]
          BC_outcome2[i] = dataout[i,"outcome2"]
          BC_sd2[i] = dataout[i,"sd2"]
          Two_Arm_VarIdx[i] = 3
        } else {
          BA_outcome1[i] = dataout[i,"outcome1"]
          BA_sd1[i] = dataout[i,"sd1"]
          BA_outcome2[i] = dataout[i,"outcome2"]
          BA_sd2[i] = dataout[i,"sd2"]
          Two_Arm_VarIdx[i] = 1
        } 
      } else {
        CA_outcome1[i] = dataout[i,"outcome1"]
        CA_sd1 = dataout[i,"sd1"]
        CA_outcome2[i] = dataout[i,"outcome2"]
        CA_sd2 = dataout[i,"sd2"]
        Two_Arm_VarIdx[i] = 2
      }
    } else {
      # Three arm study
      for(j in 1:3){
        temp_t1 = dataout[i+ (j-1)*nabc,"t1"]
        temp_t2 = dataout[i+ (j-1)*nabc,"t2"]
        if(temp_t1 == 1){
          if(temp_t2 == 2){ 
            BC_outcome1[i] = dataout[i+ (j-1)*nabc,"outcome1"]
            BC_sd1[i] = dataout[i+ (j-1)*nabc,"sd1"]
            BC_outcome2[i] = dataout[i+ (j-1)*nabc,"outcome2"]
            BC_sd2[i] = dataout[i+ (j-1)*nabc,"sd2"]
          } else {
            BA_outcome1[i] = dataout[i+ (j-1)*nabc,"outcome1"]
            BA_sd1[i] = dataout[i+ (j-1)*nabc,"sd1"]
            BA_outcome2[i] = dataout[i+ (j-1)*nabc,"outcome2"]
            BA_sd2[i] = dataout[i+ (j-1)*nabc,"sd2"]
          } 
        } else {
          CA_outcome1[i] = dataout[i+ (j-1)*nabc,"outcome1"]
          CA_sd1 = dataout[i+ (j-1)*nabc,"sd1"]
          CA_outcome2[i] = dataout[i+ (j-1)*nabc,"outcome2"]
          CA_sd2 = dataout[i+ (j-1)*nabc,"sd2"]
        }     
      }
    }
  }
  
  
  
  
  # Below Rjags code is taking reference from :
  # http://doingbayesiandataanalysis.blogspot.com/2017/06/bayesian-estimation-of-correlations-and.html
  
  y = cbind(BA_outcome1, BA_outcome2, CA_outcome1,CA_outcome2,BC_outcome1, BC_outcome2)
  ss = cbind(BA_sd1, BA_sd2, CA_sd1, CA_sd2, BC_sd1,BC_sd2)
  #y = rmvnorm(30, mean =rep(0,6), sigma = diag(6))
  colnames(y) = c("BA1", "BA2", "CA1", "CA2", "BC1", "BC2")
  colnames(ss) = c("BA1", "BA2", "CA1", "CA2", "BC1", "BC2")
  # Assemble data for sending to JAGS:
  datalist = list(
    y = y,
    Two_Arm_VarIdx = Two_Arm_VarIdx,
    ss = ss,
    Ntotal_two_arm = nab+nbc+nbc,
    Ntotal_three_arm = nabc,
    Nvar = ncol(y),
    rho_w = diag(x=1,nrow=ncol(y)),
    rho_b = diag(x=1,nrow=ncol(y))
  )
  
  
  # Define the model
  modelString = "
model {
  # Likelihood for two-arm study
  for ( i in 1: Ntotal_two_arm ) {
    Two_Arm_CovMat[i,1,1] <- ss[i,2*Two_Arm_VarIdx[i]-1]^2 + tau2[2*Two_Arm_VarIdx[i]-1]
    Two_Arm_CovMat[i,2,2] <- ss[i,2*Two_Arm_VarIdx[i]]^2 + tau2[2*Two_Arm_VarIdx[i]]
    Two_Arm_CovMat[i,1,2] <- ss[i,2*Two_Arm_VarIdx[i]-1] * ss[i,2*Two_Arm_VarIdx[i]-1]*rho_w[2*Two_Arm_VarIdx[i], 2*Two_Arm_VarIdx[i]-1] + 
                            sqrt(tau2[2*Two_Arm_VarIdx[i]-1] * tau2[2*Two_Arm_VarIdx[i]])*rho_b[2*Two_Arm_VarIdx[i], 2*Two_Arm_VarIdx[i]-1]
    Two_Arm_CovMat[i,2,1] <- Two_Arm_CovMat[i,1,2]
    Two_Arm_InvCov[i,1:2, 1:2] <- inverse(Two_Arm_CovMat[i,,])
   
    y[i,c(Two_Arm_VarIdx[i],Two_Arm_VarIdx[i]+1)] ~ dmnorm(dist_mu[c(Two_Arm_VarIdx[i], 1+Two_Arm_VarIdx[i])] ,Two_Arm_InvCov[i,,])
   

  }
  
  # Likelihood for three-arm study
  for(i in 1:Ntotal_three_arm){
      for(j in 1:6){
        for(k in 1:6){
            Three_Arm_WithinCov_Diag[i,j,k] <- ifelse(j == k, ss[i+Ntotal_two_arm,j]^2, 0)
        }
      }
      
      Three_Arm_WithinCov[i,1:Nvar, 1:Nvar] <- Three_Arm_WithinCov_Diag[i,,] %*% (rho_w + t(rho_w) + Identity_Matrix) %*% Three_Arm_WithinCov_Diag[i,,]
      Three_Arm_CovMat[i,1:Nvar, 1:Nvar] <- Three_Arm_WithinCov[i,1:Nvar, 1:Nvar]  + Three_Arm_BetweenCov
      Three_Arm_InvCov[i,1:Nvar, 1:Nvar] <- inverse(Three_Arm_CovMat[i,,])
      
      
      
      y[i+Ntotal_two_arm,] ~ dmnorm(dist_mu,Three_Arm_InvCov[i,,])
  }
  
  

  
  # Priors
  # mu
  for ( varIdx in 1:(Nvar-2) ) { 
    mu[varIdx] ~ dnorm(0,tau2[varIdx])
  }
  

  dist_mu[1:6] = c(mu,mu[1] - mu[3],mu[2] - mu[4])

  
  
  for(Idx1 in 1:Nvar){
    for(Idx2 in 1:Nvar){
      Identity_Matrix[Idx1, Idx2] <-  ifelse(Idx1 == Idx2, 1,0)
    }
  }
  
  # within study correlation
  for ( Idx1 in 2:Nvar ) {
    for ( Idx2 in 1:(Idx1-1)) {
      rho_w[Idx1,Idx2] ~ dunif(0,1)
    }
  }
 

  
  # between study correlation
  for ( Idx1 in 2:Nvar ) {
    for ( Idx2 in 1:(Idx1-1)) {
      rho_b[Idx1,Idx2] ~ dunif(0,1)
      
    }
  }


  


  # between study covariance
  #Inv_Between_Study_CovMat ~ dwish(rho_b, Nvar)
  #Between_Study_CovMat = inverse(Inv_Between_Study_CovMat)
  
   for ( Idx1 in 1:Nvar ) {
      log_tau[Idx1] ~ dunif(0,1)
      tau2[Idx1] = exp(-2*log_tau[Idx1])
   }
    
  for(j in 1:6){
        for(k in 1:6){
            Three_Arm_BetweenCov_Diag[j,k] <- ifelse(j == k, tau2[j], 0)
        }
  }
  
  Three_Arm_BetweenCov<- Three_Arm_BetweenCov_Diag %*% (rho_b + t(rho_b) + Identity_Matrix) %*% Three_Arm_BetweenCov_Diag
  
  # for(j in 1:Nvar){
  #     for(k in 1:Nvar){
  #          temp_tau[j,k]  = ifelse(j == k, exp(-2*tau[j]), 0)
  #       }
  #   }
  
 # Between_Study_CovMat[1:Nvar,1:Nvar] = temp_tau[1:Nvar, 1:Nvar] %*% rho_b %*% temp_tau[1:Nvar, 1:Nvar]
  
  #for ( varIdx in 1:Nvar ) {
  #  sigma2[varIdx] =Between_Study_CovMat[varIdx, varIdx]
  #}  
}


" # close quote for modelString
writeLines( modelString , con="Jags-MultivariateNormal-model.txt" )

# Run the chains:
nChain = 1
nAdapt = 500
nBurnIn = 500
nThin = 10
nStepToSave = 10000



jagsModel = jags.model( file="Jags-MultivariateNormal-model.txt" ,
                        data=datalist , n.chains=nChain , n.adapt=nAdapt )
update( jagsModel , n.iter=nBurnIn )
codaSamples = coda.samples( jagsModel ,
                            variable.names=c("dist_mu") ,
                            n.iter=nStepToSave/nChain*nThin , thin=nThin )
time_end = Sys.time()
difference <- difftime(time_end, start, units='mins')


posterior_mu = rbind(posterior_mu, colMeans(codaSamples[[1]]))
posterior_sd = rbind(posterior_sd, apply(codaSamples[[1]],2,sd))

print(paste(sim_t,"th Simulation using " , round(difference,4), " minutes"))

}
colnames(posterior_mu) <- c("BA1", "BA2", "CA1", "CA2", "BC1", "BC2")
colnames(posterior_sd) <- c("BA1", "BA2", "CA1", "CA2", "BC1", "BC2")

filename <- paste("simulation_result_full_likelihood_Bayesian_", Sys.Date(),".RData", sep="")
save(posterior_mu,posterior_sd,file = filename)

