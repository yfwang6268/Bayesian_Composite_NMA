
tau2_location <- function(t1, t2, narm){
  if(t2 == narm){
    tau_index = t1
  } else{
    if(t1 == 1){
      tau_index = narm - 1 +t2 - t1
    } else{
      tau_index = sum(narm - 1 - seq(1:(t1-1))) + t2-t1 + narm - 1
    }
  }
  return(tau_index)
}


initial_value_setup <- function(dataout, narm){
  # mu_vector = rep(0, narm-1)
  # for(i in 1:(narm-1)){
  #   temp_data = dataout[dataout$t1 == i & dataout$t2 == narm,]
  #   if(nrow(temp_data) >0){
  #     mu_vector[i] = mean(temp_data$outcome)
  #   }
  # }
  mu_vector = runif(narm-1)
  
  # 
  # tau2_value = ifelse(var(dataout$outcome) - mean(dataout$sd^2) > 0, 
  #                     var(dataout$outcome) - mean(dataout$sd^2),
  #                     runif(1))
  
  tau2_value = runif(1)
  return(list(mu_vector, tau2_value))
}

log_composite_likelihood <- function(y, sd, tau2, mu){

  log_likelihood =  (-1/2) * (log(sd^2+tau2) +
                             (y -mu)^2 /(sd^2+tau2)) 
  return(log_likelihood)
}


total_log_composite_likelihood <- function(dataout, mu_vector, tau2, narm){

  temp_mu = ifelse(dataout[, "t2"] == narm, mu_vector[dataout[, "t1"]], mu_vector[dataout[, "t1"]] - mu_vector[dataout[, "t2"]] )
  temp_sd = dataout[, "sd"]
  temp_y = dataout[, "outcome"]
  n = length(temp_sd)

  temp_tau2 = rep(tau2, n)
  total_lcl = sum(mapply(log_composite_likelihood, temp_y, temp_sd, temp_tau2, temp_mu))
  return(total_lcl)
}

square_root_matrix <- function(x){
  s <- svd(x)
  root <- s$u %*% sqrt(diag(s$d)) %*% t(s$v)
  return(root)
}

sample_tau2 <- function(dataout, prev_mu, prev_tau2, narm){

    number_of_tau2 <- 1
    proposed_tau2 <- rtnorm(1, prev_tau2, 1, lb = 0, ub=Inf)
    
    total_lcl_propose_tau2 <- total_log_composite_likelihood(dataout, 
                                                             prev_mu,
                                                             proposed_tau2,
                                                             narm)
    
    total_lcl_previous_tau2 <- total_log_composite_likelihood(dataout, 
                                                              prev_mu,
                                                              prev_tau2, 
                                                              narm)

    prob_proposed_tau2_on_prev_tau2 <- dnorm(proposed_tau2,  prev_tau2, 1, log=TRUE) - log(1 - pnorm(0, prev_tau2, 1))
    prob_prev_tau2_on_proposed_tau2 <- dnorm(prev_tau2, proposed_tau2, 1, log=TRUE) - log(1 - pnorm(0, proposed_tau2, 1))
    prior_proposed_tau2 <- dinvgamma(proposed_tau2, shape = 0.001, rate = 0.001, log=TRUE)
    prior_prev_tau2 <- dinvgamma(prev_tau2, shape = 0.001,  rate = 0.001, log=TRUE)
    numerator = total_lcl_propose_tau2 + prior_proposed_tau2 - prob_proposed_tau2_on_prev_tau2
    denominator = total_lcl_previous_tau2 + prior_prev_tau2 - prob_prev_tau2_on_proposed_tau2
    
    
    
    if(log(runif(1)) < min(0, numerator - denominator)){
      prev_tau2 = proposed_tau2
    } 
  return(prev_tau2)
}


sample_mu <- function(dataout, prev_mu, prev_tau2, narm){
  number_of_mu <- length(prev_mu)
  for(i in 1:number_of_mu){
    
    # studies comparing trt i versus the placebo
    set_A <- dataout[(dataout$t1 == i) & (dataout$t2 == narm),]
    if(nrow(set_A) > 0){
      G1 = sum(1/(2*(set_A$sd^2+prev_tau2)))
      H1 = sum((set_A$outcome)/(2*(set_A$sd^2+prev_tau2))) 

    } else {
      G1 = 0
      H1 = 0
    }
    
    # studies comparing trt i versus other trts
    set_B <- dataout[(dataout$t1 == i) & (dataout$t2 <  narm),]
    if(nrow(set_B) > 0){
      G2 = sum(1/(2*(set_B$sd^2+prev_tau2)))
      H2 = sum((set_B$outcome + prev_mu[set_B$t2])/(2*(set_B$sd^2+prev_tau2)))
    } else {
      G2 = 0
      H2 = 0
    }
    
    
    # studies comparing other trts versus trt i
    set_C <- dataout[(dataout$t1 < i) & (dataout$t2 == i),]
    if(nrow(set_C) > 0){
      G3 = sum(1/(2*(set_C$sd^2+prev_tau2)))
      H3 = -sum((set_C$outcome - prev_mu[set_C$t1])/(2*(set_C$sd^2+prev_tau2)))
    } else {
      G3 = 0
      H3 = 0
    }
    
    G = G1+G2+G3
    H = H1+H2+H3
    
    prev_mu[i] = rnorm(1, H/G, sqrt(1/(2*G)))
    
  }
  return(prev_mu)
}


Gibbs_Sampler_Overall <- function(dataout, chain_length, burn_in_rate, narm){

  initial_values = initial_value_setup(dataout, narm)
  prev_mu = initial_values[[1]]
 
  prev_tau2 = initial_values[[2]]

  estimated_mu = matrix(nrow = chain_length, ncol = (narm - 1))
  estimated_tau2 = matrix(nrow = chain_length, ncol = length(prev_tau2))
  for(t in 1:chain_length){
    prev_tau2 = sample_tau2(dataout, prev_mu, prev_tau2, narm)  
    prev_mu = sample_mu(dataout, prev_mu, prev_tau2, narm)   
    estimated_mu[t,] = prev_mu
    estimated_tau2[t,] = prev_tau2
  }
  start_index = ceiling(chain_length * burn_in_rate)
  estimated_mu = estimated_mu[start_index: chain_length,]
  estimated_tau2 = estimated_tau2[start_index:chain_length,]
 return(list(estimated_mu,estimated_tau2))
}

wrap_up_total_log_composite_likelihood <- function(dataout, theta, narm){
  mu = theta[1:(narm-1)]
  tau = theta[narm]
  result <- total_log_composite_likelihood(dataout, mu, tau, narm)
  return(result)
}


wrap_up_posterior <- function(dataout, theta, narm){

  result <-wrap_up_total_log_composite_likelihood(dataout, theta, narm)

  # prior for tau
  tau2 = theta[narm]
  result = result + dinvgamma(tau2, shape = 0.001, rate = 0.001, log=TRUE)
  # prior for mu
  mu_vector = theta[1:(narm-1)]
  result = result + sum(dnorm(mu_vector,0, 100, log=TRUE))
  return(result)
} 

calculate_hessian_matrix <- function(dataout, theta, narm = 3,d=1e-10){
  number_of_parameter <- length(theta)
  indicator_matrix <- diag(number_of_parameter)
  result <- matrix(0, nrow = number_of_parameter, ncol = number_of_parameter)
  log_posterior0 <- wrap_up_posterior(dataout, theta, narm)
  
  for(i in 1:number_of_parameter){
    for(j in 1:number_of_parameter){
      theta1 <- theta + indicator_matrix[i,]*d + indicator_matrix[j,]*d
      log_posterior1 <- wrap_up_posterior (dataout, theta1, narm)
      theta2 <- theta + indicator_matrix[i,]*d
      log_posterior2 <- wrap_up_posterior (dataout, theta2, narm)
      theta3 <- theta + indicator_matrix[j,]*d
      log_posterior3 <- wrap_up_posterior(dataout, theta3, narm)
      result[i,j] = ((log_posterior1 - log_posterior2) -
                     (log_posterior3 - log_posterior0))/d^2
    }
  }
  result = (result + t(result))/2
  return(result)
}

calculate_score_vector <- function(dataout, theta, narm = 3,d=1e-10){
  number_of_parameter <- length(theta)
  indicator_matrix <- diag(number_of_parameter)
  result <- matrix(0, nrow = number_of_parameter, ncol = 1)
  log_posterior0 <- wrap_up_posterior(dataout, theta, narm)
  for(i in 1:number_of_parameter){
    log_posterior1 <-  wrap_up_posterior(dataout, theta + d*indicator_matrix[i,], narm)
    result[i,1] <- (log_posterior1 - log_posterior0)/d
  }
  return(result)  
}

calculate_score_squares <- function(data_list, narm_vector, theta,
                                    d=1e-10){
  
  number_of_parameters = sum(narm_vector)
  number_of_outcomes = length(data_list)
  SS = matrix(0, nrow = number_of_parameters, ncol = number_of_parameters)
  cum_start_index = 0
  for(i in 1:number_of_outcomes){
    temp_data = data_list[[i]]
    temp_theta = theta[cum_start_index + 1:narm_vector[i]]
    for(j in 1:nrow(temp_data)){
      temp_S = matrix(0,nrow = number_of_parameters,ncol = 1)

      temp_S[cum_start_index + 1:narm_vector[i],1] = calculate_score_vector(temp_data[j,], 
                                                                            temp_theta, 
                                                                            narm_vector[i],d)
      SS = SS + c(temp_S) %*% t(c(temp_S))
    }
    cum_start_index = cum_start_index + narm_vector[i]
  }
  return(SS)
}

  
  
  
  
  
  
  
