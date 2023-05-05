log_composite_likelihood <- function(dataout, tau, mu,t1, t2, k,narm = 3){
  log_likelihood <- 0
  temp_data = subset(dataout, dataout$t1 == t1 & dataout$t2 == t2)
  if(k == 1){
    observed_effects  = temp_data$outcome1
    within_study_variance = temp_data$sd1^2
  } else {
    observed_effects = temp_data$outcome2
    within_study_variance = temp_data$sd2^2
  }

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

    numerator = exp(log_composite_likelihood(dataout, between_study_variance, mu, k, t1, t2, narm = 3))/between_study_variance/rnorm(between_study_variance, prev_tau2, 0.1)
    denominator = exp(log_composite_likelihood(dataout, prev_tau2, mu, k, t1, t2, narm = 3))/prev_tau2/prev_tau2/rnorm(prev_tau2, between_study_variance, 0.1)
    #print(prev_tau2)
    #print(c(pdf_proposed_distribution(between_study_variance, within_study_variance, observed_effects), dgamma(between_study_variance, prev_tau2)))
    a_mh_algorithm = numerator / denominator
    if(runif(1) <= min(a_mh_algorithm ,1)){
      tau2[t] = between_study_variance
      prev_tau2 = between_study_variance
    } else {
      tau2[t] = prev_tau2
    }
  }
  start_index = ceiling(mc_length*burn_in_rate)
  result = mean(tau2[start_index:mc_length])
  return(result)
}

  sample_mu <- function(tau2, observed_effects, within_study_variance){
    mu_variance = 1/(sum(1/(within_study_variance + tau2)))
    return(rnorm(1, mean(observed_effects), sqrt(mu_variance)))
  }


Gibbs_Sampler_Individual <- function(dataout, chain_length, burn_in_rate, k, t1, t2 ,narm = 3){

  temp_data = subset(dataout, dataout$t1 == t1 & dataout$t2 == t2)
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
  
  for(t = 1:chain_length){
    prev_mu = sample_mu(prev_tau2, observed_effects, within_study_variance)
    prev_tau2 = mh_algorithm(observed_effects, within_study_variance, k, mc_length,t1, t2, prev_mu, burn_in_rate)
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
    simulation_mu = matrix(nrow = chain_length, ncol = 6)
    simulation_tau = simulation_mu
    for(k in 1:2){
      for (t in 1:chain_length){
        for(t1 in 1:(narm-1)){
          for(t2 in (t1+1):narm){
            if(t1 == 1){
              if(t2 == 3){
                tau_index = 1 # BA
              } else {
                tau_index = 3 # BC
              }
            } else {
              tau_index = 2 # CA
            }
            sim_result = Gibbs_Sampler_Individual(dataout, chain_length, burn_in_rate, k, t1, t2 ,narm = 3)
            simulation_mu[t, (k-1)*3+tau_index] = mean(sim_result[[1]])
            simulation_tau[t, (k-1)*3+tau_index] = mean(sim_result[[2]])
          }
        }
      }
    }

   return(list(simulation_mu, simulation_tau))
  }
  
  
  
  
  

  
  
  
  
  
  
  
