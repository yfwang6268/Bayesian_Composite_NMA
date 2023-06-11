library(gemtc)
library(rjags)
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
rho = matrix(c(1,0.5,0.1,0.1,0.5,1,0.1,0.1,0.1,0.1,1,0.5,0.1,0.1,0.5,1),
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
nab = nac=nbc=nabc=20

# mcmc setting
n.chain = 4
n.adapt = 500
n.iter = 10000
thin = 10

estimated_mu1 = NULL
estimated_mu2 = NULL
estimated_ci1 = NULL
estimated_ci2 = NULL

for(t in 1:simulation_times){
  start <- Sys.time()
  
  # prepare the dataset
  dataout = gendata(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2)
  dataout1 = subset(dataout,  select = c("ID", "t1", "t2", "outcome1", "sd1"))
  dataout1 = prepare_dataset_for_gemtc(dataout1, nab, nac, nbc, nabc)
  dataout2 = subset(dataout,  select = c("ID", "t1", "t2", "outcome2", "sd2"))
  dataout2 = prepare_dataset_for_gemtc(dataout2, nab, nac, nbc, nabc)
  
  # std.err of the reference arm must < std.err of the relative effects
  valid_dataou1 = valid_std_err(dataout1)
  valid_dataou2 = valid_std_err(dataout2)
  
  while(!(valid_dataou1 & valid_dataou2)){
    dataout = gendata(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2)
    dataout1 = subset(dataout,  select = c("ID", "t1", "t2", "outcome1", "sd1"))
    dataout1 = prepare_dataset_for_gemtc(dataout1, nab, nac, nbc, nabc)
    dataout2 = subset(dataout,  select = c("ID", "t1", "t2", "outcome2", "sd2"))
    dataout2 = prepare_dataset_for_gemtc(dataout2, nab, nac, nbc, nabc)
    valid_dataou1 = valid_std_err(dataout1)
    valid_dataou2 = valid_std_err(dataout2)
  }

 # Network graph  
  network1 <- mtc.network(data.re = dataout1)
  network2 <- mtc.network(data.re = dataout2)
  
  # model compilation
  model1 <- mtc.model(network1, 
                      linearModel = "random",
                      n.chain = n.chain)
  
  model2 <- mtc.model(network2, 
                      linearModel = "random",
                      n.chain = n.chain)
  
  
  # mcmc simulation
  mcmc1 <- mtc.run(model1, n.adapt = n.adapt, n.iter = n.iter, thin = thin)
  mcmc2 <- mtc.run(model2, n.adapt = n.adapt, n.iter = n.iter, thin = thin)
  
  
  #  Generating the estimated mu  results
  temp_result1 = relative.effect(mcmc1, t1=3)
  temp_result2 = relative.effect(mcmc2, t1=3)
  temp_samplers1 = NULL
  temp_samplers2 = NULL
  temp_ci1 = NULL
  temp_ci2 = NULL
  
  for(i in 1:n.chain){
    temp_samplers1 = rbind(temp_samplers1, temp_result1$samples[[i]])
    temp_samplers2 = rbind(temp_samplers2, temp_result2$samples[[i]])
  }
  temp_estiamted_mu1 = colMeans(temp_samplers1)
  temp_estiamted_mu1[3] = temp_estiamted_mu1[1] - temp_estiamted_mu1[2]
  estimated_mu1 = rbind(estimated_mu1, temp_estiamted_mu1)
  temp_estiamted_mu2 = colMeans(temp_samplers2)
  temp_estiamted_mu2[3] = temp_estiamted_mu2[1] - temp_estiamted_mu2[2] 
  estimated_mu2 = rbind(estimated_mu2, temp_estiamted_mu2)

  for(i in 1:ncol(estimated_mu1)){
    temp_lb = quantile(temp_samplers1[,i], 0.025)
    temp_ub = quantile(temp_samplers1[,i], 0.975)
    temp_ci1 = c(temp_ci1, temp_lb, temp_ub)
  }
  estimated_ci1 = rbind(estimated_ci1, temp_ci1)

  for(i in 1:ncol(estimated_mu2)){
    temp_lb = quantile(temp_samplers2[,i], 0.025)
    temp_ub = quantile(temp_samplers2[,i], 0.975)
    temp_ci2 = c(temp_ci2, temp_lb, temp_ub)
  }  
  estimated_ci2 = rbind(estimated_ci2, temp_ci2)
  
  print(paste("Simulation ", t, " is done using ", round(Sys.time() - start, 4), " seconds"))
}

filename <- paste("gemtc_sim_unequal_correlation_large_sample_",  Sys.Date(),".RData", sep="")
save(estimated_mu1, estimated_ci1, estimated_mu2,estimated_ci2, file = filename)
