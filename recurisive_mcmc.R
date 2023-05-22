log_composite_likelihood <- function(observed_effects, within_study_variance, tau, mu){
  log_likelihood <- 0
  
  between_study_variance = tau
  estimated_effect_size = mu
  log_likelihood = sum(log(within_study_variance+between_study_variance)) +
                   sum((observed_effects -estimated_effect_size)^2/
                         (within_study_variance+between_study_variance))  
  return(log_likelihood*(-1/2))
}


mh_algorithm <- function(observed_effects, within_study_variance, k, mc_length,t1, t2, mu, burn_in_rate){

  tau2 = numeric(mc_length)

  prev_tau2 = 0.3


  for(t in 1:mc_length){
    between_study_variance = -1
    while(between_study_variance <= 0){
     between_study_variance = rnorm(1, within_study_variance, 0.1)
    }
    
#    between_study_variance = runif(1, prev_tau2 - 0.1, prev_tau2 + 0.1)

    numerator = exp(log_composite_likelihood(observed_effects, within_study_variance, between_study_variance, mu))/between_study_variance/rnorm(between_study_variance, prev_tau2, 0.1)*(1-pnorm(0, prev_tau2, 0.1))
    denominator = exp(log_composite_likelihood(observed_effects, within_study_variance, prev_tau2, mu))/prev_tau2/rnorm(prev_tau2, between_study_variance, 0.1)*(1-pnorm(0, between_study_variance, 0.1))

    a_mh_algorithm = numerator / denominator
    if(runif(1) <= min(a_mh_algorithm ,1)){
      tau2[t] = between_study_variance
      prev_tau2 = between_study_variance
    } else {
      tau2[t] = prev_tau2
    }
  }
  start_index = ceiling(mc_length*burn_in_rate)
  result = tau2[mc_length]
  return(result)
}

  sample_mu <- function(tau2, observed_effects, within_study_variance){
    n = length(observed_effects)
    mu_variance = 1/(sum(1/(within_study_variance + tau2)))
    return(rnorm(1, mean(observed_effects), sqrt(mu_variance/n)))
  }


Gibbs_Sampler_Individual <- function(dataout, chain_length, burn_in_rate, k, t_1, t_2 ,narm = 3){

  temp_data = subset(dataout, t1 == t_1 & t2 == t_2)
  if(k == 1){
    observed_effects  = temp_data$outcome1
    within_study_variance = temp_data$sd1^2
  } else {
    observed_effects = temp_data$outcome2
    within_study_variance = temp_data$sd2^2
  }
  
  simulation_tau2 = numeric(chain_length - ceiling(chain_length * burn_in_rate) + 1)
  simulation_mu =  simulation_tau2
  prev_tau2 = 0.3

  for(t in 1:chain_length){
    prev_mu = sample_mu(prev_tau2, observed_effects, within_study_variance)
    prev_tau2 = mh_algorithm(observed_effects, within_study_variance, k, 500,t1, t2, prev_mu, burn_in_rate)
    simulation_mu[t] = prev_mu
    simulation_tau2[t] = prev_tau2

  }

  store_row = ceiling(chain_length * burn_in_rate)
  simulation_tau = simulation_tau2[store_row:chain_length]
  simulation_mu = simulation_mu[store_row:chain_length]

  result = list(simulation_mu, simulation_tau)
  return(result)
  
  }
  
  
  Gibbs_Sampler_Overall <- function(dataout, chain_length, burn_in_rate, narm = 3){
    simulation_mu = matrix(nrow = 1, ncol = 6)
    simulation_tau = simulation_mu
    for(k in 1:2){
        for(t1 in 1:(narm-1)){
          for(t2 in (t1+1):narm){
            start_time <- Sys.time()
            if(t1 == 1){
              if(t2 == 3){
                tau_index = 1 # BA
              } else {
                tau_index = 3 # BC
              }
            } else {
              tau_index = 2 # CA
            }
            
            if (tau_index != 3){
              sim_result = Gibbs_Sampler_Individual(dataout, chain_length, burn_in_rate, k, t1, t2 ,narm = 3)
              simulation_mu[1, (k-1)*3+tau_index] = mean(sim_result[[1]])
              simulation_tau[1, (k-1)*3+tau_index] = mean(sim_result[[2]])
              
            }
           
          }

        }   
   
    }
   return(c(simulation_mu, simulation_tau))
  }
  
  
  
  
  

  
  
  
  
  
  
  
