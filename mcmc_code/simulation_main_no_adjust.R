
#source("adjusted_metropolis_hasting_algorithm.R")
source("recurisive_mcmc.R")
#source("log_composite_likelihood.R")
source("sampling_from_proposed_distribution.R")
source("data_simulation.R")
source("CLNMA_functions.R")
library(parallel)
library(MASS)


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
nab = nac=nbc=nabc=5



chain_length <- 5000

burn_in_rate <- 0.5

number_of_simulation <- 1

# Single Core Running

posterior_mu <- NULL
posterior_mu_var_posterior <- NULL
posterior_mu_var_sandwich <- NULL
adjustment.method <- NULL

for(t in 1:number_of_simulation){
  start <- Sys.time()
  temp_simulated_data <- gendata(nab, nac, nbc, nabc, mu1, mu2, betweenv, rho_w, ss1,ss2)
  temp_posterior <- Gibbs_Sampler_Overall(temp_simulated_data , chain_length, burn_in_rate, adjustment.method = adjustment.method)
  estimated_mu = rbind(posterior_mu, temp_posterior[1:6])
  posterior_mu_var_posterior= rbind(posterior_mu_var_posterior, temp_posterior[7:12])
  posterior_mu_var_sandwich = rbind(posterior_mu_var_sandwich, temp_posterior[13:18])
  print(paste("Simulation ", t, " is done using ", round(Sys.time() - start, 4), " seconds"))
}

filename <- paste("simulation_result_no_adjustment", Sys.Date(),".RData", sep="")
save(posterior_mu,estimated_mu_var_posterior,estimated_mu_var_sandwich,file = filename)
