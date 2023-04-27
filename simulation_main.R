print(utils::getSrcDirectory(function(){}))
setwd(getSrcDirectory(function(){})[1])
source("adjusted_metropolis_hasting_algorithm.R")
source("log_composite_likelihood.R")
source("sampling_from_proposed_distribution.R")
source("data_simulation.R")
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

number_of_simulation <- 500

posterior_mu <- list()
posterior_tau <- list()
simulated_data <- list()

for(t in 1:number_of_simulation){
  start <- Sys.time()
  temp_simulated_data <- simulate_dataset(nab, nac, nbc, nabc, mu1, mu2, betweenv, rho_w, ss1,ss2)
  temp_posterior <- adjusted_metropolis_hasting_algorithm(temp_simulated_data , chain_length, burn_in_rate)
  simulated_data = append(simulated_data, temp_simulated_data)
  posterior_mu = append(posterior_mu, temp_posterior[1])
  posterior_tau = append(posterior_tau, temp_posterior[2])    
  print(paste("Simulation ", t, " is done using ", round(Sys.time() - start, 4), " seconds"))
}

filename <- paste("simulation_result_number_of_simulations_date_", number_of_simulation,"_", Sys.Date(),".RData", sep="")
save(posterior_mu, posterior_tau,simulated_data,file = filename)


true_mu = c(0.5,1,-0.5, 0,-0.5,0.5)
#number_of_simulation = 409
bias = numeric(6)
for(i in 1:number_of_simulation){
  temp_mu = posterior_mu[[i]]
  bias =bias +  abs(true_mu -  colMeans(temp_mu))
}
print(round(bias/number_of_simulation,3))
  

