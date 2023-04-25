library(invgamma)

pdf_propose_distribution <- function(between_study_variance, within_study_variance, observed_effects){
  
  part1 = 1
  n_divided_by_sigma2 = 0
  for(s in within_study_variance){
    part1 = part1 * 1/sqrt(s+between_study_variance)
    n_divided_by_sigma2 =  n_divided_by_sigma2 + 1/(s+between_study_variance)
   
  }
  part2 = 0
  ybar = mean(observed_effects)
  n = length(observed_effects)
  for(i in 1:n){
    part2 = part2 - (observed_effects[i] - ybar)^2/(2*(between_study_variance+within_study_variance[i]))
  }
  part2 = exp(part2)
  part3 = sqrt(2*pi/n_divided_by_sigma2)
  return(part1*part2*part3)
}

sampling_from_proposed_distribution <- function(between_study_variance, within_study_variance, observed_effects){
  n = length(observed_effects)
  proposed_value = rinvchisq(1, n, between_study_variance)
 
  numerator = pdf_propose_distribution(proposed_value, within_study_variance, observed_effects)*dinvchisq(between_study_variance, 
                                                                                                          n, proposed_value)
  denominator = pdf_propose_distribution(between_study_variance, within_study_variance, observed_effects)*dinvchisq(proposed_value, n, 
                                                                                                                    between_study_variance)
  
  
  
  a_proposed_distribution = numerator/denominator
  
  if(runif(1) <= min(a_proposed_distribution, 1)){
    return(proposed_value)
  } else{
    return(between_study_variance)
  }
}

product_propose_distribution <- function(dataout, tau, mu, narm = 3){
  product_probability = 1
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
        product_probability = product_probability  * pdf_propose_distribution(between_study_variance, 
                                                                              within_study_variance, 
                                                                              observed_effects)
      }
    }
  }
  return(product_probability)
}



product_prior_distribution <- function(dataout, tau, narm = 3){
  
  product_probability = 1
  for(k in 1:2){
    if(k == 1){
      temp_tau = tau[1:narm]
    } else {
      temp_tau = tau[narm+1:narm]
    }
    for(t1 in 1:(narm-1)){
      for(t2 in (t1+1):narm){
        temp_data = subset(dataout, dataout$t1 == t1 & dataout$t2 == t2)
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
        n = dim(temp_data)[1]
        product_probability = product_probability  * dinvchisq(1, n, between_study_variance)
      }
    }
  }  
  return(product_probability)
}
