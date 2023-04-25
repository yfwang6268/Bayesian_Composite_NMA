sampling_mu_and_tau <- function(dataout, prev_tau=NULL, narm = 3){
  proposed_mu = NULL
  proposed_tau = NULL


  if(is.null(prev_tau)){
      temp_prev_tau = c(0.2,0.3,0.5,0.3,0.1,0.4)
  } else {
      temp_prev_tau = prev_tau
  }

  
  for(k in 1:2){
    w = numeric(3)
    yw = numeric(3)
    tau_result = numeric(3)

    for(t1 in 1:(narm-1)){
      for(t2 in (t1+1):narm){
        
        temp_data = subset(dataout, dataout$t1 == t1 & dataout$t2 == t2)
        if(k == 1){
          observed_effects  = temp_data$outcome1
          within_study_variance = temp_data$sd1^2
        } else {
          observed_effects = temp_data$outcome2
          within_study_variance = temp_data$sd2^2
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
 
        between_study_variance = temp_prev_tau[(k-1)*3 + tau_index]
        

        if(is.null(temp_prev_tau)){
            tau_result = temp_prev_tau
          } else {
          if(tau_result[tau_index] == 0){
            tau_result[tau_index] = sampling_from_proposed_distribution(between_study_variance, 
                                                                      within_study_variance, 
                                                                      observed_effects)         
            } 
          }
 

       
        w[tau_index] = sum((within_study_variance + tau_result[tau_index])^(-1))
        yw[tau_index] = sum(observed_effects/((within_study_variance + tau_result[tau_index])^(-1)))
        
      }
      
    }
    
    H = matrix(c(w[1]+w[3], -w[3], -w[3], w[2]+w[3]), nrow = 2)
    v = matrix(c(yw[1]-yw[3],yw[2]-yw[3]), nrow = 2)
    mu_k = solve(H, v) 
    mu_k = c(mu_k, mu_k[2] - mu_k[1])
    proposed_mu = c(proposed_mu, mu_k)
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
      proposed_result = sampling_mu_and_tau(dataout)
      previous_mu = proposed_result[[1]]
      previous_tau = proposed_result[[2]]
      simulation_mu[t,] = previous_mu
      simulation_tau[t,] = previous_tau

    } else {
      
      # proposed new parameter value
      
      proposed_result = sampling_mu_and_tau(dataout, previous_tau)
      
      proposed_mu = proposed_result[[1]]
      proposed_tau = proposed_result[[2]]

      

      # calcuate a
      likelihood_ratio  = log_composite_likelihood(dataout, proposed_tau, proposed_mu) /
        log_composite_likelihood(dataout, previous_tau, previous_mu)
      proposed_dist_ratio = product_propose_distribution(dataout, previous_tau, previous_mu) /
        product_propose_distribution(dataout, proposed_tau, proposed_mu)
      prior_dist_ratio = product_prior_distribution(dataout, proposed_tau) /
        product_prior_distribution(dataout, previous_tau)
      
      a_mh_algo = likelihood_ratio * proposed_dist_ratio * prior_dist_ratio  
 
      # accpet/reject
      if(runif(1) <= min(a_mh_algo, 1)){
        simulation_mu[t,] = proposed_mu 
        simulation_tau[t,] = proposed_tau
        previous_mu = proposed_mu
        previous_tau = proposed_tau
      } else {
        simulation_mu[t,] = previous_mu
        simulation_tau[t,] = previous_tau
      }
    }
  }
  
  store_row = floor(chain_length * burn_in_rate) + 1
  
  simulation_mu = simulation_mu[store_row:chain_length,]
  simulation_tau = simulation_tau[store_row:chain_length,]
  result = list(simulation_mu, simulation_tau)
  return(result)
  
  }
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
