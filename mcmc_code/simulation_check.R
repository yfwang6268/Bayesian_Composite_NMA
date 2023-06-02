# This script is to check the correctness of the simulation data

# I will simulate the data using my script and the data using Duan's code

# Then I will use the same CLNMA function to estimate the mean and the variance

# If I simulate the data correctly, I will get the same result.

source("CLNMA_functions.R")
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


# Duan's dataset
set.seed(1)
dataout_duan = gendata(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2)
result_duan = tryCatch(CLNMA(dataout_duan),error=function(e) rep(NA,20))
print(result_duan)
# My dataset
set.seed(1)
dataout_mine = simulate_dataset(nab, nac, nbc, nabc, mu1, mu2, betweenv, rho_w, ss1,ss2)
result_mine = tryCatch(CLNMA(dataout_mine),error=function(e) rep(NA,20))
print(result_mine)
