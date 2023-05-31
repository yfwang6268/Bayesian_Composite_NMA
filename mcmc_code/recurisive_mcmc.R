log_composite_likelihood <- function(observed_effects, within_study_variance, tau, mu){
  log_likelihood <- 0
  
  between_study_variance = tau
  estimated_effect_size = mu
  log_likelihood = sum(log(within_study_variance+between_study_variance)) +
                   sum((observed_effects -estimated_effect_size)^2/
                         (within_study_variance+between_study_variance))  
  return(log_likelihood*(-1/2))
}


sample_tau2 <- function(observed_effects, within_study_variance, k, mc_length,t1, t2, mu, burn_in_rate, adjustment = NULL){
  prev_tau2 = 0.3
  for(t in 1:mc_length){
    between_study_variance = -1
    while(between_study_variance <= 0){
     between_study_variance = rnorm(1, within_study_variance, 0.1)
    }
    
#    between_study_variance = runif(1, prev_tau2 - 0.1, prev_tau2 + 0.1)
    
    if(is.null(adjustment)){
      numerator = exp(log_composite_likelihood(observed_effects, within_study_variance, between_study_variance, mu))/between_study_variance/dnorm(between_study_variance, prev_tau2, 0.1)*(1-pnorm(0, prev_tau2, 0.1))
      denominator = exp(log_composite_likelihood(observed_effects, within_study_variance, prev_tau2, mu))/prev_tau2/dnorm(prev_tau2, between_study_variance, 0.1)*(1-pnorm(0, between_study_variance, 0.1))
    } else if (adjustment == "magn"){
      H_matrix = hessian_matrix(observed_effects, prev_tau2, within_study_variance, mu)
      J_matrix = sandwich_J_matrix(observed_effects, prev_tau2, within_study_variance, mu)
      #ev =  eigen(solve(H_matrix) %*% J_matrix)
      Hinverse_times_J = solve(H_matrix) %*% J_matrix
      #k_adjust = length(ev$values)/sum(ev$values)
      k_adjust = dim(Hinverse_times_J)[2] / sum(diag(Hinverse_times_J))
      numerator = exp(k_adjust * log_composite_likelihood(observed_effects, within_study_variance, between_study_variance, mu))/between_study_variance/dnorm(between_study_variance, prev_tau2, 0.1)*(1-pnorm(0,prev_tau2 , 0.1))
      denominator = exp(k_adjust * log_composite_likelihood(observed_effects, within_study_variance, prev_tau2, mu))/prev_tau2/dnorm(prev_tau2, between_study_variance, 0.1)*(1-pnorm(between_study_variance, prev_tau2 , 0.1))
    }
    a_mh_tau = numerator / denominator
    
    if(!is.na(a_mh_tau)){
      if(runif(1) <= min(a_mh_tau, 1)){
        prev_tau2 = between_study_variance
      } 
    }
  }
  return(prev_tau2 )
}

sample_mu <- function(tau2, observed_effects, within_study_variance, chain.length = 500 ,adjustment = NULL){
    n = length(observed_effects)
    mu_variance = 1/(sum(1/(within_study_variance + tau2)))
    if(is.null(adjustment)){
      return(rnorm(1, mean(observed_effects), sqrt(mu_variance/n)))
    } else if (adjustment == "magn"){
      prev_mu = 0.3
      for(t in 1:chain.length){
        H_matrix = hessian_matrix(observed_effects, tau2, within_study_variance, prev_mu)
        J_matrix = sandwich_J_matrix(observed_effects, tau2, within_study_variance, prev_mu)
        Hinverse_times_J = solve(H_matrix) %*% J_matrix
      
        k = dim(Hinverse_times_J)[2]/sum(diag(Hinverse_times_J))
        proposed_mu = rnorm(1, prev_mu, 0.1)
        numerator = dnorm(proposed_mu, mean(observed_effects), sqrt(mu_variance/n))^k /
                    dnorm(proposed_mu, prev_mu, 0.1)
        denominator = dnorm(prev_mu, mean(observed_effects), sqrt(mu_variance/n))^k /
                      dnorm(prev_mu, proposed_mu, 0.1)
        a_mh_mu = numerator/ denominator
        if(!is.na(a_mh_mu)) {
          if(runif(1) <= min(a_mh_mu, 1)){
            prev_mu = proposed_mu
          }
        }
      }
      return(prev_mu)
    }
  }


Gibbs_Sampler_Individual <- function(dataout, chain_length, burn_in_rate, k, t_1, t_2 ,narm = 3, adjustment.method = NULL){

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
    #start = Sys.time()
    prev_mu = sample_mu(prev_tau2, observed_effects, within_study_variance, adjustment = adjustment.method)
    prev_tau2 = sample_tau2(observed_effects, within_study_variance, k, 500,t1, t2, prev_mu, burn_in_rate, adjustment =  adjustment.method)
    simulation_mu[t] = prev_mu
    simulation_tau2[t] = prev_tau2
    #print(paste("Step ", t, " is done using ", round(Sys.time() - start, 4), " seconds"))
  }

  store_row = ceiling(chain_length * burn_in_rate)
  simulation_tau = simulation_tau2[store_row:chain_length]
  simulation_mu = simulation_mu[store_row:chain_length]

  H_matrix = hessian_matrix(observed_effects, mean(simulation_tau), within_study_variance, mean(simulation_mu))
  J_matrix = sandwich_J_matrix(observed_effects, mean(simulation_tau), within_study_variance, mean(simulation_mu))
  sandwich_matirx = solve(H_matrix) %*% J_matrix %*% solve(H_matrix) / length(observed_effects)
  sandwich_var_mu = sandwich_matirx[1,1]
  sandwich_var_tau2 =   sandwich_matirx[2,2]
  
  result = list(simulation_mu, simulation_tau, sandwich_var_mu, sandwich_var_tau2)
  return(result)
  
  }
  
  
  Gibbs_Sampler_Overall <- function(dataout, chain_length, burn_in_rate, adjustment.method = NULL, narm = 3){
    simulation_mu = matrix(nrow = 1, ncol = 6)
    simulation_tau = simulation_mu
    variance_mu_mcmc = simulation_mu
    variance_tau_mcmc = simulation_mu
    variance_mu_sandwich = simulation_mu
    variance_tau_sandwich = simulation_tau
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
              sim_result = Gibbs_Sampler_Individual(dataout, chain_length, burn_in_rate, k, t1, t2 ,narm = 3, adjustment = adjustment.method)
              simulation_mu[1, (k-1)*3+tau_index] = mean(sim_result[[1]])
              variance_mu_mcmc[1, (k-1)*3+tau_index] = var(sim_result[[1]])
              variance_mu_sandwich[1, (k-1)*3+tau_index] = sim_result[[3]]
              simulation_tau[1, (k-1)*3+tau_index] = mean(sim_result[[2]])
              variance_tau_mcmc[1, (k-1)*3+tau_index] = var(sim_result[[2]])
              variance_tau_sandwich[1, (k-1)*3+tau_index] = sim_result[[4]]
              

            }
           
          }

        }   
   
    }
   return(c(simulation_mu, simulation_tau, 
            variance_mu_mcmc, variance_tau_mcmc, 
            variance_mu_sandwich, variance_tau_sandwich))
  }
  
  
  
  
  

  
  
  
  
  
  
  
