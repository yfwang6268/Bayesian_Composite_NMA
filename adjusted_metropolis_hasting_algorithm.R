sampling_mu_and_tau <- function(dataout, prev_tau, narm = 3){
  proposed_mu = NULL
  proposed_tau = NULL


  for(k in 1:2){
    tau_result = numeric(3)
    mu_result =numeric(3)
    for(t1 in 1:(narm-1)){
      for(t2 in (t1+1):narm){
        temp_data = subset(dataout, dataout$t1 == t1 & dataout$t2 == t2)
        if(k == 1){
          observed_effects  = temp_data$outcome1
          within_study_variance = temp_data$sd1^2
        } else {
          observed_effects = temp_data$outcome2
          within_study_variancence = temp_data$sd2^2
        }
        
        if(t1 == 1){
          if(t2 == 3){
            tau_index = 1 # BA
          } else {
            tau_index = 3 # BC
          }
        } else {
          tau_index = 2 # CA
        }
 
        between_study_variance = prev_tau[(k-1)*3 + tau_index]

        n = length(within_study_variance )

        if(tau_result[tau_index] == 0){
          tau_result[tau_index] = sampling_from_proposed_distribution(between_study_variance, 
                                                                      within_study_variance, 
                                                                      observed_effects)
          temp_variance = 1/sum(1/(within_study_variance + tau_result[tau_index]))         
          mu_result[tau_index] = rnorm(1, mean(observed_effects), sqrt(temp_variance))
          } 
      }
      
    }

    proposed_mu = c(proposed_mu, mu_result)
    proposed_tau = c(proposed_tau, tau_result)
    
  }
  return(list(proposed_mu, proposed_tau))
}


adjusted_metropolis_hasting_algorithm <- function(dataout, chain_length, 
                                                  burn_in_rate, narm = 3){
  
  
  simulation_mu = matrix(nrow = chain_length,ncol = 6)
  simulation_tau = simulation_mu
  
  for(t in 1:chain_length){
    if(t == 1){
      previous_mu = c(0.4,0.8,-0.4,0.1,-0.5,0.6)
      previous_tau = c(0.2,0.3,0.451,0.3,0.1,0.365)
      simulation_mu[t,] = previous_mu
      simulation_tau[t,] = previous_tau

    } else {
      
      # proposed new parameter value
      
      proposed_result = sampling_mu_and_tau(dataout, previous_tau)
      
      simulation_mu[t,] = proposed_result[[1]]
      simulation_tau[t,] = proposed_result[[2]]

      previous_mu = proposed_result[[1]]
      previous_tau = proposed_result[[2]]
    }
  }
  
  store_row = floor(chain_length * burn_in_rate) + 1
  
  simulation_mu = simulation_mu[store_row:chain_length,]
  simulation_tau = simulation_tau[store_row:chain_length,]
  result = list(simulation_mu, simulation_tau)
  return(result)
  
  }
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
