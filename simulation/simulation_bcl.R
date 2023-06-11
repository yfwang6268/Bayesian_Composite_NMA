source("recurisive_mcmc.R")
source("sampling_from_proposed_distribution.R")
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



chain_length <- 1000

burn_in_rate <- 0.5

number_of_simulation <- 1

# Single Core Running

posterior_mu <- NULL
posterior_mu_var_posterior <- NULL
posterior_mu_var_sandwich <- NULL
posterior_mu_ci <- NULL

for(t in 1:number_of_simulation){
  start <- Sys.time()
  temp_simulated_data <- gendata(nab, nac, nbc, nabc, mu1, mu2, betweenv, rho_w, ss1,ss2)
  temp_posterior <- Gibbs_Sampler_Overall(temp_simulated_data , chain_length, burn_in_rate)
 # simulated_data = append(simulated_data, temp_simulated_data)
  estimated_mu = rbind(posterior_mu, temp_posterior[1:6])
  posterior_mu_var_posterior= rbind(posterior_mu_var_posterior, temp_posterior[7:12])
  posterior_mu_var_sandwich = rbind(posterior_mu_var_sandwich, temp_posterior[13:18])
  estimated_ci = rbind(posterior_mu_ci, temp_posterior[(length(temp_posterior)-11):length(temp_posterior)])
  print(paste("Simulation ", t, " is done using ", round((Sys.time() - start)/60, 4), " minutes"))
}


# simulation_result_matrix =  matrix(simulation_result, ncol = 12, byrow = T)
filename <- paste("simu_blc_no_adjustment_", Sys.Date(),".RData", sep="")
save(estimated_mu,posterior_mu_var_posterior,posterior_mu_var_sandwich,estimated_ci,file = filename)

