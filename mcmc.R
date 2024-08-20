rm(list = ls())
library(Matrix)
source("gibbs_within_mh_equal_tau.R")
library(invgamma)
library(expm)
library(mvtnorm)
library(dplyr)
library(plyr)
library(MASS)
library(TruncatedNormal)

run_gibbs_within_mh <- function(new_data, chain_length, burn_in_rate, narm, thin){
  

  # Run mh within gibbs
  sampler <- Gibbs_Sampler_Overall(new_data, chain_length, burn_in_rate, narm)
  
  # add thin rate
  posterior_mu = sampler[[1]]
  n = dim(posterior_mu)[1]
  posterior_mu = posterior_mu[c(seq(1,n,thin),n),]
  posterior_tau = sampler[[2]]
  posterior_tau = posterior_tau[c(seq(1,n,thin),n)]
  posterior = cbind(posterior_mu, posterior_tau)
  return(posterior)
  
}

mcmc <- function(all_data, chain_length = 50000, burn_in_rate = 0.2, thin = 5, number_of_outcomes = 4){
  
  posterior = NULL
  narm_vector = numeric(number_of_outcomes)
  H = list()
  data_by_list = list()
  theta = NULL
  
 for(i in 1:number_of_outcomes){

    temp_col_index = 4:5 + 2 *(i-1)
    input_data = all_data[,c(1:3,temp_col_index)]

    colnames(input_data) = c("StudyId","t1","t2","outcome","sd")
    input_data = input_data[!is.na(input_data$outcome),]
    # print(input_data)
    temp_trt_no <- sort(unique(c(input_data$t1, input_data$t2)))
    temp_narm <- length(temp_trt_no)
    
    input_data$t1 <- mapvalues(input_data$t1, from = temp_trt_no, to = 1:temp_narm)
    input_data$t2 <- mapvalues(input_data$t2, from = temp_trt_no, to = 1:temp_narm)
    
    n = dim(input_data)[1]
    for(j in 1:n){
      if(input_data$t1[j] > input_data$t2[j]){
        temp_trt = input_data$t1[j]
        input_data$t1[j] = input_data$t2[j]
        input_data$t2[j] = temp_trt
        input_data$outcome[j] = -input_data$outcome[j]
      }
    }
    
    # rearrange the order to make sure t1 < t2
    
    
    temp_t1 <- ifelse(input_data$t1 > input_data$t2, input_data$t2, input_data$t1) 
    temp_t2 <- ifelse(input_data$t1 > input_data$t2, input_data$t1, input_data$t2) 
    temp_outcome <- ifelse(input_data$t1 > input_data$t2, -input_data$outcome, input_data$outcome) 
    input_data$t1 = temp_t1
    input_data$t2 = temp_t2
    input_data$outcome = temp_outcome
    
    
    
    temp_posterior <- run_gibbs_within_mh(input_data, chain_length, burn_in_rate, temp_narm, thin)
 
    temp_theta <- colMeans(temp_posterior)
    
    narm_vector[i] <-  temp_narm
  
    data_by_list[[i]] <- input_data
    
    
    theta <- c(theta, temp_theta)
    data_by_list[[i]] <- input_data
    H[[i]] <- calculate_hessian_matrix(input_data, temp_theta ,temp_narm, 1e-5)

    posterior <- cbind(posterior, temp_posterior)
    print(i)
  }
  
  H <- as.matrix(bdiag(H))
  SS <- calculate_score_squares(data_by_list,  narm_vector, theta)
  
  BS_mean_vector = colMeans(posterior)
  BS_mean_matrix = matrix(rep(BS_mean_vector,nrow(posterior)), nrow = nrow(posterior),byrow = T)
  temp = t(solve(-H) %*% square_root_matrix(SS) %*% square_root_matrix(-H) %*% t(posterior - BS_mean_matrix) + t(BS_mean_matrix))
  adjust_posterior = tryCatch(t(solve(-H) %*% square_root_matrix(SS) %*% square_root_matrix(-H) %*% t(posterior - BS_mean_matrix) + t(BS_mean_matrix)),
                              error = function(e){posterior})
  
  
  return(list(posterior,adjust_posterior, narm_vector))
}

seed = 1
set.seed(seed)
ori_data = read.csv("clean_data.csv")
#ori_data = ori_data[-13,]

result <- mcmc(ori_data)
filename <- paste("sample_seed", seed, ".RData",sep="")
save.image(filename)

