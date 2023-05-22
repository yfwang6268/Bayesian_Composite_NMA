library(invgamma)

pdf_proposed_distribution <- function(between_study_variance, within_study_variance, observed_effects){
  part1 = prod(sqrt(within_study_variance+between_study_variance)^(-1))
  ybar = mean(observed_effects)
  n = length(observed_effects)
  part2 = exp(-sum((observed_effects - ybar)^2/(2*(within_study_variance+between_study_variance))))
  part3 = sqrt(1/sum(1/(within_study_variance+between_study_variance))*2*pi)
  return(part1*part2*part3*(between_study_variance)^(-1))
}



sampling_from_proposed_distribution <- function(within_study_variance, observed_effects, simu_length=1000){
  n = length(observed_effects)
  simulated_value = numeric(simu_length + 1)
  count = 0
  for(t in 1:(simu_length+ 1)){
      #proposed_value = rinvgamma(1, 1)
      phi = runif(1)
      proposed_value = exp(phi)
      if(t > 1){
        numerator = pdf_proposed_distribution(proposed_value, within_study_variance, observed_effects)
        denominator = pdf_proposed_distribution(previous_value , within_study_variance, observed_effects)
        a_proposed_distribution = numerator/denominator
        if(runif(1) <= min(a_proposed_distribution, 1)){
          simulated_value[t] = proposed_value
          previous_value = proposed_value
          } else {
          simulated_value[t] = previous_value
          }
      } else {
        simulated_value[t] = proposed_value
        previous_value = proposed_value
      }
  
  } 
  start_index = floor(simu_length * 0.5) + 1
  return(mean(simulated_value[start_index:(simu_length + 1)]))
}

product_propose_distribution <- function(dataout, tau, mu, k, narm = 3){
    product_probability = 1

    temp_tau = tau
    temp_mu = mu

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
        effect_size = temp_mu[comparision_index]
        temp_mu_variance = 1/sum(1/(within_study_variance + between_study_variance))


        product_probability = product_probability  * 
                              pdf_proposed_distribution(between_study_variance, within_study_variance, observed_effects) *
                              dnorm(effect_size, mean(observed_effects), sqrt(temp_mu_variance))
      }
    }
  
  return(product_probability)
}



product_prior_distribution <- function(dataout, tau, mu, narm = 3){
  
  product_probability = 1
  for(k in 1:2){
    temp_tau = tau[(k-1)*3 + 1:narm]
    temp_mu = mu[(k-1)*3 + 1:narm]
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
        effect_size = temp_mu[comparision_index]
        n = dim(temp_data)[1]  
        product_probability = product_probability  * dgamma(1/between_study_variance,1) * dnorm(effect_size)
      }
    }
  }  
  return(product_probability)
}
