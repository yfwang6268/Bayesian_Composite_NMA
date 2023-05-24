hessian_matrix <- function(observed_effects, between_study_variance, within_study_variance, estimated_effect){

  n = length(observed_effects)
  
  
  second_partial_log_likelihood_over_partial_tau2 = sum(1/(between_study_variance + within_study_variance)/2 +
                                                            (observed_effects - estimated_effect)^2/(between_study_variance + within_study_variance)^3) 
                                                      
    
  second_partial_log_likelihood_over_partial_mu = sum( -1/(between_study_variance + within_study_variance))
                                                    
                                                      
  second_partial_log_likelihood_over_partial_mu_and_tau2 = sum(-(observed_effects - estimated_effect)/(between_study_variance + within_study_variance)^2)
                                                            
    
  
  
  result = matrix(c(second_partial_log_likelihood_over_partial_mu,
                    second_partial_log_likelihood_over_partial_mu_and_tau2,
                    second_partial_log_likelihood_over_partial_mu_and_tau2, 
                    second_partial_log_likelihood_over_partial_tau2),
                  nrow = 2)
  
  return(result)
}



sandwich_J_matrix <- function(observed_effects, between_study_variance, within_study_variance, estimated_effect){
  # first_partial_log_likelihood_over_partial_tau2 = 0
  # first_partial_log_likelihood_over_partial_mu = 0
  # n = length(observed_effects)
  # for(i in 1:n){
    first_partial_log_likelihood_over_partial_tau2 = sum(- 
                                                     1/(between_study_variance + within_study_variance)/2 +
                                                     (observed_effects - estimated_effect)^2/(between_study_variance + within_study_variance)^2/2)
    
    first_partial_log_likelihood_over_partial_mu = sum(
                                                   (observed_effects - estimated_effect)/(between_study_variance + within_study_variance))
    
  # }
  
  result = matrix(c(first_partial_log_likelihood_over_partial_mu^2,
                    first_partial_log_likelihood_over_partial_mu * first_partial_log_likelihood_over_partial_tau2,
                    first_partial_log_likelihood_over_partial_mu * first_partial_log_likelihood_over_partial_tau2,
                    first_partial_log_likelihood_over_partial_tau2^2),
                  nrow =  2)
  
  return(result)
}