hessian_matrix <- function(observed_effects, between_study_variance, within_study_variance, estimated_effect){
  second_partial_log_likelihood_over_partial_tau2 = 0
  second_partial_log_likelihood_over_partial_mu = 0
  second_partial_log_likelihood_over_partial_mu_and_tau = 0
  
  
  for(y in observed_effects){
    second_partial_log_likelihood_over_partial_tau2 = second_partial_log_likelihood_over_partial_tau2 + 
                                                      1/(between_study_variance + within_study_variance)/2 +
                                                      (y - estimated_effect)^2/(between_study_variance + within_study_variance)^3
    
    second_partial_log_likelihood_over_partial_mu = second_partial_log_likelihood_over_partial_mu -
                                                    1/(between_study_variance + within_study_variance)
                                                      
    second_partial_log_likelihood_over_partial_mu_and_tau2 = second_partial_log_likelihood_over_partial_mu_and_tau -
                                                            (y - estiamted_effect)/(between_study_variance + within_study_variance)^2
    
  }
  
  result = matrix(c(second_partial_log_likelihood_over_partial_mu,
                    second_partial_log_likelihood_over_partial_mu_and_tau2,
                    second_partial_log_likelihood_over_partial_mu_and_tau2, 
                    second_partial_log_likelihood_over_partial_tau2),
                  nrow = 2)
  
  return(result)
}



J_matrix <- function(observed_effects, between_study_variance, within_study_variance, estimated_effect){
  first_partial_log_likelihood_over_partial_tau2 = 0
  first_partial_log_likelihood_over_partial_mu = 0
  
  for(y in observed_effects){
    first_partial_log_likelihood_over_partial_tau2 = first_partial_log_likelihood_over_partial_tau2 - 
                                                     1/(between_study_variance + within_study_variance)/2 +
                                                     (y - estimated_effect)^2/(between_study_variance + within_study_variance)^2/2
    
    first_partial_log_likelihood_over_partial_mu = first_partial_log_likelihood_over_partial_mu +
                                                   (y - estimated_effect)/(between_study_variance + within_study_variance)
    
  }
  
  result = matrix(c(first_partial_log_likelihood_over_partial_mu^2,
                    first_partial_log_likelihood_over_partial_mu * first_partial_log_likelihood_over_partial_tau2,
                    first_partial_log_likelihood_over_partial_mu * first_partial_log_likelihood_over_partial_tau2,
                    first_partial_log_likelihood_over_partial_tau2^2),
                  nrow =  2)
  
  return(result)
}