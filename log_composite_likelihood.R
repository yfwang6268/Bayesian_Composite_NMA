log_composite_likelihood <- function(dataout, tau, mu ,narm = 3){
  log_likelihood <- 0
  for(k in 1:2){
    if(k == 1){
      temp_tau = tau[1:narm]
      temp_mu = mu[1:narm]
    } else {
      temp_tau = tau[narm+1:narm]
      temp_mu = mu[narm+1:narm]
    }
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
            comparision_index = 1 # BA
          } else {
            comparision_index = 3 # BC
          }
        } else {
          comparision_index = 2 # CA
        }
        between_study_variance = temp_tau[comparision_index]
        estimated_effect_size = temp_mu[comparision_index]
        log_likelihood = log_likelihood + 
                         sum(log(within_study_variance+between_study_variance)) +
                         sum((observed_effects -estimated_effect_size)^2/
                               (within_study_variance+between_study_variance))
        
      }
    }
  }
  return(log_likelihood)
}