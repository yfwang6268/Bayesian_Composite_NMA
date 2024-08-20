rm(list=ls())
library(pcnetmeta)
library(wesanderson)

trt_names_tbl <- read.csv("trt_mapping.csv", header=TRUE)
trt_names <- trt_names_tbl$trt_name

unique_outcomes <- c("CRC", "RecRed", "CALgain", "KTgain")

trt_names_mapping <- function(origin_index, narm, trt_names){
  trt_index <- mapvalues(1:narm, from = 1:narm, to = origin_index)
  result <- trt_names[trt_index]
  result <- result[-length(result)] # remove the placebo
  return(result)
}


study_id_with_trt_id <- function(real_data){
  Unique_StudyId <- unique(as.character(real_data$StudyId))
  s.id <- t.id <- NULL
  for(i in 1:length(Unique_StudyId)){
    temp_subset_data <- real_data[real_data$StudyId == Unique_StudyId[i],]
    temp.trt <- sort(unique(c(temp_subset_data$t1, temp_subset_data$t2)))
    s.id <- c(s.id, rep(i, length(temp.trt)))
    t.id <- c(t.id, temp.trt)
  }
  return(list(s.id, t.id))
}

network_plot_function<- function(real_data, outcome, clr){
  temp_study_and_trt_ids <- study_id_with_trt_id(real_data)
  nma.networkplot(s.id=temp_study_and_trt_ids[[1]], t.id=temp_study_and_trt_ids[[2]], n=NULL, 
                  text.cex = 2,weight.node = FALSE,
                  node.col = clr, title = outcome)
}

number_of_outcomes = length(unique_outcomes)
pal <- wes_palette(name = "Zissou1",n = number_of_outcomes, type = c("continuous"))
pal <- pal[sample(1:number_of_outcomes,number_of_outcomes ,replace = FALSE)]
#par(mar = c(2, 5, 2, 5))
m <- matrix(c(1:4,rep(5,2)),nrow = 3,ncol = 2, byrow = T)

layout(mat = m,widths = rep(0.5,2), heights = c(0.4,0.4,0.2))
#par(mar=c(2,1,1,1))
for(i in 1:number_of_outcomes){
  temp_outcome = unique_outcomes[i]
  temp_data = read.csv(paste(temp_outcome,".csv", sep =""))
  colnames(temp_data) <- c("StudyId", "t1", "t2", "outcome", "sd")
  par(mar=c(2,2,2,2))
  network_plot_function(temp_data, temp_outcome, pal[i])
}
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
length_name <- NULL
for(i in 1:length(trt_names)){
  temp_trt_name = trt_names[i]
  length_name <- c(length_name, paste("(",i, ") ", temp_trt_name, sep=""))
}

legend("center",inset = 0 ,legend = length_name, ncol = 4,cex =1.5,x.intersp = 0, title = "Treatment")