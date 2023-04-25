library(MASS)

simulate_observed_effects <- function(mu1, mu2,between_study_covariance, within_study_correlation, within_study_sd1,within_study_sd2 ){
  
  observed_within_study_sd1 = abs(rnorm(2,within_study_sd1, 0.1))
  observed_within_study_sd2 = abs(rnorm(2,within_study_sd2, 0.1))
  
  observed_within_study_variance = c(observed_within_study_sd1, observed_within_study_sd2)^2
  observed_within_study_covariance = diag(observed_within_study_variance)
  
  
  for(i in 1:3){
    for(j in (i+1):4){
      observed_within_study_covariance[i,j] = observed_within_study_covariance[j,i] = within_study_correlation[i,j]*
                                                                                  sqrt(observed_within_study_variance[i] *
                                                                                       observed_within_study_variance[j])
    }
  }
  
  y = mvrnorm(1, c(mu1, mu2), observed_within_study_covariance + between_study_covariance)
  
  
  return(c(y, observed_within_study_sd1, observed_within_study_sd2))
}


subset_dataset <- function(data, arm,start_index, end_index,t1, t2){
  n = end_index - start_index + 1
  result = data.frame(ID = start_index:end_index,t1 = rep(t1, n), t2 = t2)
  outcome1 = paste("Y_",arm,"_1", sep="")
  sd1 = paste("S_",arm,"_1", sep="")
  outcome2 = paste("Y_",arm,"_2", sep="")
  sd2 = paste("S_",arm,"_2", sep="")  
  result$outcome1 = data[start_index:end_index, outcome1]
  result$outcome2 = data[start_index:end_index, outcome2]
  result$sd1 = data[start_index:end_index, sd1]
  result$sd2 = data[start_index:end_index, sd2]
  return(result)
}
  


simulate_dataset <- function(nab, nac, nbc, nabc, mu1, mu2, 
                             between_study_covariance, within_study_correlation, 
                             within_study_sd1,within_study_sd2) {
n = nab + nac + nabc + nabc

dataset = replicate(n, simulate_observed_effects(mu1, mu2,
                                                 between_study_covariance, within_study_correlation, 
                                                 within_study_sd1,within_study_sd2 ))

dataset = data.frame(t(dataset))

colnames(dataset) = c("Y_BA_1","Y_CA_1","Y_BA_2","Y_CA_2","S_BA_1","S_CA_1","S_BA_2","S_CA_2")

dataset$Y_BC_1 = dataset$Y_BA_1 - dataset$Y_CA_1
dataset$Y_BC_2 = dataset$Y_BA_2 - dataset$Y_CA_2
dataset$S_BC_1 = sqrt(dataset$S_BA_1^2+dataset$S_CA_1^2 - within_study_correlation[1,2]*dataset$S_BA_1*dataset$S_CA_1)
dataset$S_BC_2 = sqrt(dataset$S_BA_2^2+dataset$S_CA_2^2 - within_study_correlation[3,4]*dataset$S_BA_2*dataset$S_CA_2)


#BA Arm
dataset_BA = subset_dataset(dataset, "BA",1,nab,1, 3)
dataset_CA = subset_dataset(dataset, "CA",nab+1,nab+nac,2, 3)
dataset_BC = subset_dataset(dataset, "BC",nab+nac+1,nab+nac+nbc,1, 2)
dataset_ABC_BA = subset_dataset(dataset,"BA" ,nab+nac+nbc + 1,nab+nac+nbc + nabc,1, 3)
dataset_ABC_CA = subset_dataset(dataset, "CA",nab+nac+nbc + 1,nab+nac+nbc + nabc,2, 3)
dataset_ABC_BC = subset_dataset(dataset, "BC",nab+nac+nbc + 1,nab+nac+nbc + nabc,1, 2)

dataout = rbind(dataset_BA, dataset_CA, dataset_BC, dataset_ABC_BA, dataset_ABC_CA, dataset_ABC_BC)
dataout$ID = seq(1, nab+nac+nbc+3*nabc)

return(dataout)
  
}