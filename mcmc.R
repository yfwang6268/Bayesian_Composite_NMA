mh_algorithm <- function(dataout, k, mc_length,t1, t2, burn_in_rate){

  tau2 = numeric(mc_length)

  prev_tau2 = 0.3

  temp_data = subset(dataout, dataout$t1 == t1 & dataout$t2 == t2)
  if(k == 1){
    observed_effects  = temp_data$outcome1
    within_study_variance = temp_data$sd1^2
  } else {
    observed_effects = temp_data$outcome2
    within_study_variance = temp_data$sd2^2
  }

  for(t in 1:mc_length){
#    between_study_variance = -1
#    while(between_study_variance <= 0){
#      between_study_variance = rnorm(1, prev_tau2, 0.1)
#    }
    between_study_variance = runif(1, prev_tau2 - 0.1, prev_tau2 + 0.1)    
    numerator = pdf_proposed_distribution(between_study_variance, within_study_variance, observed_effects)#/dnorm(1,between_study_variance, prev_tau2,0.1)
    denominator = pdf_proposed_distribution(prev_tau2, within_study_variance, observed_effects)#/dnorm(prev_tau2, between_study_variance,0.1)
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
  tau2 = tau2[ceiling(chain_length*burn_in_rate):chain_length]
  return(tau2)
}

sample_mu <- function(tau2, observed_effects, within_study_variance){
  mu_variance = 1/(sum(1/(within_study_variance + tau2)))
  return(rnorm(1, mean(observed_effects), sqrt(mu_variance)))
}

Gibbs_Sampler_Individual <- function(dataout, chain_length, burn_in_rate, k, t1, t2, narm = 3){

  temp_data = subset(dataout, dataout$t1 == t1 & dataout$t2 == t2)
  if(k == 1){
    observed_effects  = temp_data$outcome1
    within_study_variance = temp_data$sd1^2
  } else {
    observed_effects = temp_data$outcome2
    within_study_variance = temp_data$sd2^2
  }
  
  
  simulation_tau = mh_algorithm(dataout,k,chain_length, t1,t2, burn_in_rate)


  simulation_mu = sapply(simulation_tau, sample_mu, observed_effects = observed_effects, within_study_variance = within_study_variance) 
  result = list(simulation_mu, simulation_tau)
  return(result)
  
}
  
  
Gibbs_Sampler_Overall <- function(dataout, chain_length, burn_in_rate, narm = 3){
  start_index = ceiling(chain_length * burn_in_rate)
  simulation_mu = matrix(nrow = chain_length - start_index + 1, ncol = 6)
  simulation_tau = simulation_mu
  for(k in 1:2){

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
          sim_result = Gibbs_Sampler_Individual(dataout, chain_length, burn_in_rate,k, t1, t2 ,narm = 3)
          simulation_mu[, (k-1)*3+tau_index] = mean(sim_result[[1]])
          simulation_tau[, (k-1)*3+tau_index] = mean(sim_result[[2]])
        }
      }

  }

 return(list(simulation_mu, simulation_tau))
}

  
  
  
  

  
  
  
  
  
  
  
