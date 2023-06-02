library(netmeta)
library(MASS)

source("CLNMA_functions.R")
source("functions_gemtc.R")

simulation_times = 500
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
nab = nac=nbc=nabc=5

estimated_mu1 = NULL
estimated_mu2 = NULL


for(t in 1:simulation_times){
  # simulate data
  dataout = gendata(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2)
  dataout1 = subset(dataout, select = c("ID", "outcome1", "sd1", "t1", "t2"))
  colnames(dataout1) = c("author", "TE", "seTE", "treat1", "treat2")
  dataout2 = subset(dataout, select = c("ID", "outcome2", "sd2", "t1", "t2"))
  colnames(dataout2) = c("author", "TE", "seTE", "treat1", "treat2")
  
  
  ## Model fitting
  m.netmeta1<- netmeta(TE = TE,
                       seTE = seTE,
                       treat1 = treat1,
                       treat2 = treat2,
                       studlab = author,
                       data = dataout1,
                       sm = "SMD",
                       fixed = FALSE,
                       random = TRUE,
                       reference.group = 3,
                       details.chkmultiarm = TRUE,
                       sep.trts = " vs ")
  
  m.netmeta2<- netmeta(TE = TE,
                       seTE = seTE,
                       treat1 = treat1,
                       treat2 = treat2,
                       studlab = author,
                       data = dataout2,
                       sm = "SMD",
                       fixed = FALSE,
                       random = TRUE,
                       reference.group = 3,
                       details.chkmultiarm = TRUE,
                       sep.trts = " vs ")
  
  # extract relative effects
  temp_estimated_mu1 = m.netmeta1$TE.random[-nrow(m.netmeta1$TE.random),3]
  estimated_mu1 = rbind(estimated_mu1,temp_estimated_mu1)
  temp_estimated_mu2 = m.netmeta2$TE.random[-nrow(m.netmeta2$TE.random),3]
  estimated_mu2 = rbind(estimated_mu2,temp_estimated_mu2)
  
  #print(paste("Simulation ", t, " is done using ", round((Sys.time() - start)/60, 4), " seconds"))
}

filename <- paste("netmeta_simulation_high_correlation_",  Sys.Date(),".RData", sep="")
save(estimated_mu1, estimated_mu2,file = filename)
