mh_algorithm <- function(dataout, k, mc_length,t1, t2){

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
    between_study_variance = rnorm(1, prev_tau2, 0.1)
    numerator = pdf_proposed_distribution(between_study_variance, within_study_variance, observed_effects)/dnorm(between_study_variance, prev_tau2, 0.1)
    denominator = pdf_proposed_distribution(prev_tau2, within_study_variance, observed_effects)/dnorm(prev_tau2, between_study_variance, 0.1)
    a_mh_algorithm = numerator / denominator
    if(runif(1) <= min(a,1)){
      tau2[t] = between_study_variance
      prev_tau2 = between_study_variance
    } else {
      tau2[t] = prev_tau2
    }
  }
  return(tau2)
}


Gibbs_Sampler <- function(dataout, chain_length, burn_in_rate, k, t1, t2 ,narm = 3){

  temp_data = subset(dataout, dataout$t1 == t1 & dataout$t2 == t2)
  if(k == 1){
    observed_effects  = temp_data$outcome1
    within_study_variance = temp_data$sd1^2
  } else {
    observed_effects = temp_data$outcome2
    within_study_variance = temp_data$sd2^2
  }
  
  simulation_tau2 = numeric(chain_length)
  
  simulation_tau = mh_algorithm(dataout,k,chain_length, t1,t2)
  store_row = floor(chain_length * burn_in_rate) + 1
  simulation_tau = simulation_tau[store_row:chain_length]


  simulation_mu = numeric(chain_length - store_row + 1)
  sample_mu <- function(tau2, observed_effects, within_study_variance){
    mu_variance = 1/(sum(1/(within_study_variance + tau2)))
    return(rnorm(1, mean(observed_effects), sqrt(mu_variance)))
  }
  simulation_mu = sapply(simulation_tau, sample_mu, observed_effects = observed_effects, within_study_variance = within_study_variance) 
  result = list(simulation_mu, simulation_tau)
  return(result)
  
  }
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
