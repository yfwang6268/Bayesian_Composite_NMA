rm(list=ls())
library(dplyr)
library(plyr)
library(kableExtra)
library(microViz)
library(cochrane)

my_colors = distinct_palette(14,pal = "kelly")
unique_outcomes <- c("CRC", "RecRed", "CALgain", "KTgain")
trt.mapping <- read.csv("trt_mapping.csv")
posterior_file_name = "sample_seed"

calculate_probability<- function(order_vector, narm){
  rank_proportion = sapply(seq(1:(narm-1)), 
                           function(v, vector)mean(v == vector), 
                           vector= order_vector)
  return(rank_proportion)
}

calculate_sucra <- function(dataout, narm){
  ranking_result = t(apply(dataout, 1, rank))
  narm=dim(ranking_result)[2] + 1
  ranking_probability = apply(ranking_result, 2, calculate_probability, narm=narm)
  sucra_result = colSums(apply(ranking_probability[1:(narm-2),],2,cumsum))/(narm-2)
  return(sucra_result)
}

temp_filename <- paste(posterior_file_name,  1, ".RData", sep="")
load(temp_filename)
adjust_posterior <- result[[2]]

number_of_outcomes <- length(unique_outcomes)
temp_location <- 0

sucra_result <- matrix(ncol = number_of_outcomes,
                       nrow = nrow(trt.mapping) - 1)
unique_trts = trt.mapping$trt_name[1:(nrow(trt.mapping) - 1)]
colnames(sucra_result) = unique_outcomes
rownames(sucra_result) = unique_trts 
sucra_result <- data.frame(sucra_result)

m <- matrix(c(1:4,5,5),nrow = 3,ncol = 2, byrow = T)

layout(mat = m,widths = rep(0.5,2), heights = c(0.4,0.4,0.2))

for(i in 1:number_of_outcomes){
  temp_outcome = unique_outcomes[i]
  temp_file_name <- paste(temp_outcome, ".csv", sep="")
  temp_data <- read.csv(temp_file_name)    

  temp_trt_no <- sort(unique(c(temp_data$t1, temp_data$t2)))
  temp_narm <- length(temp_trt_no)
  
  temp_origin_trt_no <- mapvalues(1:temp_narm, from = 1:temp_narm, to = temp_trt_no)
  temp_trt_names_result <- trt.mapping$trt_name[temp_origin_trt_no]

  temp_adjust_posterior <- adjust_posterior[,temp_location + 1:(temp_narm - 1)]
  temp_location = temp_location + temp_narm
  colnames(temp_adjust_posterior) <- temp_trt_names_result[-temp_narm]
  
  temp_sucra_result <- calculate_sucra(-temp_adjust_posterior,temp_narm)
  temp_col_index = which(unique_outcomes == temp_outcome)
  for(trt in names(temp_sucra_result)){
    temp_row_index = which(unique_trts == trt)
    sucra_result[temp_row_index, temp_col_index] = temp_sucra_result[trt]
  }
  
  par(mar=c(5,5,5,5))
  
  temp_sucra_result = sucra_result[, temp_col_index]
  temp_sucra_result[is.na(temp_sucra_result)] = 0
  temp_trt_no = order(temp_sucra_result, decreasing = T)
  temp_plot_valaue = temp_sucra_result[temp_trt_no]
  
  temp_bar <- barplot(names = temp_trt_no,
                      height = temp_plot_valaue,
                      main = temp_outcome,
                      ylim = c(0,1.19),
                      col = my_colors[i],
                      cex.main=1.5)
  
  text(temp_bar, temp_plot_valaue+0.1 , round(temp_plot_valaue,2),
       cex=1.75) 
  
}

mtext("Rank",side=1,line=-15,outer=TRUE,cex=1.5,
      font = 2)
mtext("Probability",side=2,line=-2,outer=TRUE,cex=1.5,
      font=2)

par(mar=c(0,4,1,4))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
length_name <- NULL

for(i in c(1,5,2,6,3,7,4)){
  temp_trt_name =  unique_trts[i]
  length_name <- c(length_name, paste("(",i, ") ",temp_trt_name, sep=""))
}

legend("top",x.intersp=-0.5, y.intersp = 1, bty = "y", cex=2.25, 
       legend = length_name, ncol = 4, title = "Treatment")

write.csv(sucra_result, file="SUCRA_result.csv")




