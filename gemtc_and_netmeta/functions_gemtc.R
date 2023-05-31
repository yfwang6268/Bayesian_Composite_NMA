calculate_reference_std <- function(dataset){
  var_BA = dataset[16:20, "sd"]^2
  var_CA = dataset[21:25, "sd"]^2
  var_BC = dataset[21:25, "sd"]^2
  sd_result = sqrt((var_BA + var_CA - var_BC)/2)
  study_id =  dataset[16:20, "ID"]
  result = cbind(study_id, sd_result)
  colnames(result) = c("ID", "sd") 
  return(result)
}

prepare_dataset_for_gemtc <- function(dataset){
  colnames(dataset)[4:5] = c("outcome","sd")
  reference_sd = calculate_reference_std(dataset)
  dataset = subset(dataset,  t2 != 2)
  result = NULL
  for(i in 1:nrow(dataset)){
    temp_Id = dataset[i, "ID"]
    temp_treatment=  dataset[i,"t1"]
    temp_reference = dataset[i,"t2"]
    temp_diff = dataset[i, "outcome"]
    temp_sd = dataset[i,"sd"]
    result = rbind(result, c(temp_Id, temp_treatment, temp_diff, temp_sd))
    # gemtc only allows one reference level per study
    if(sum(is.na(result[result[,1] == temp_Id,3])) == 0){
      
      
      if(temp_Id %in% reference_sd[,"ID"]){
        
        temp_sd = reference_sd[reference_sd[,"ID"] == temp_Id,"sd"]
        
        
      } else {
        temp_sd = NA
      }
      result = rbind(result, c(temp_Id, temp_reference, NA,temp_sd))
    }  
  }
  colnames(result) = c("study", "treatment", "diff", "std.err")
  return(result)
}

valid_std_err <- function(dataset){
  std.err.ref <- subset(dataset, is.na(dataset[,"diff"]) & !is.na(dataset[,"std.err"]))  
  std.err.rel <- subset(dataset, 
                        dataset[,"study"] %in% std.err.ref[,"study"] & !is.na(dataset[,"diff"]) &  !is.na(dataset[,"std.err"]))

  for (id in unique(std.err.ref[,"study"])){
     std.err.ref.value <- subset(std.err.ref, std.err.ref[,"study"] == id, select = std.err)
     std.err.rel.value.array <- subset(std.err.rel, std.err.rel[,"study"] == id, select = std.err)
     for (s in std.err.rel.value.array){
       if (s <= std.err.ref.value){
         return(FALSE)
       }
     }
  }
  return(TRUE)
}


