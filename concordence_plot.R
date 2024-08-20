rm(list=ls())
library(microViz)
my_colors = distinct_palette(14,pal = "kelly")

paper_result <- read.csv("paper_estimate.csv")
ori_data <- read.csv("clean_data.csv")
trt_name <- c("CAF+EMD", "CAF+ADM", "CAF+PRP","CAF+HF-DDS",
               "CAF+CM", "CAF","CAF+CTG","CAF+BM")
poster <- NULL
adj_poster <- NULL


temp_filename = paste("sample_seed", 1, ".RData", sep = "")
load(temp_filename)
poster <- result[[1]]
adj_poster <- result[[2]]


narm_vector <- result[[3]]
cum_location_col = 0
m <- matrix(c(1:4),nrow = 2,ncol = 2, byrow = T)
outcome_name <- c("CRC", "RecRed", "CAL gain", "KT gain")
layout(mat = m,widths = rep(0.2,5), heights = rep(0.33,3))


zvalue_plot <- function(paper_result, before_adjust_result, after_adjust_result,
                        outcome, outcome_number,color_palette = my_colors){
  x <- after_adjust_result
  y <- paper_result
  
  plot(x, y,
       xlim = c(min(x) - 5,max(x) + 5), xlab="",
       ylim = c(min(y) - 5,max(y) + 5), ylab="",
       main = outcome,
       col = color_palette[outcome_number],
       pch=19)
  points(before_adjust_result, paper_result, pch = 1, 
         col = color_palette[outcome_number])
  
  abline(h=c(qnorm(0.025), qnorm(0.975)), lty = "dashed", col = "red")
  abline(v=c(qnorm(0.025), qnorm(0.975)), lty = "dashed", col = "red")
}

m <- matrix(c(1:4,rep(5,2)),nrow = 3,ncol = 2, byrow = T)

layout(mat = m,
       widths = c(rep(0.5,2)), 
       heights = c(rep(0.45,2),0.1)) 



 for(i in 1:4){

  temp_narm = narm_vector[i]
  temp_poster = poster[,1:temp_narm + cum_location_col]
  temp_poster = temp_poster[,1:(temp_narm-1)]
  temp_adj_poster = adj_poster[,1:temp_narm + cum_location_col]
  temp_adj_poster = temp_adj_poster[,1:(temp_narm-1)]
  
  temp_zvalue =  colMeans(temp_poster)/apply(temp_poster,2,sd)
  temp_adj_zvalue = colMeans(temp_adj_poster)/apply(temp_adj_poster,2,sd)
  cum_location_col = temp_narm + cum_location_col
  
  temp_paper_col_inds <- c(1,2,3:5+(i-1)*3)
  temp_ori_col_inds <- c(2,3,4:5+(i-1)*2)
  
  temp_paper_result <- paper_result[,temp_paper_col_inds]
  colnames(temp_paper_result) = c("t1", "t2", "outcome", "lb", "up")
  temp_paper_result = temp_paper_result[!is.na(temp_paper_result$outcome),]
  
  temp_ori_data <- ori_data[,temp_ori_col_inds]
  colnames(temp_ori_data) = c("t1", "t2", "outcome", "sd")
  temp_ori_data = temp_ori_data[!is.na(temp_ori_data$outcome),]
  temp_trt_no = sort(unique(c(temp_ori_data$t1, temp_ori_data$t2)))
  
  if(i == 1){
    temp_paper_sd <- (log(temp_paper_result[,5]) - log(temp_paper_result[,4]))/2/qnorm(0.95)
    temp_paper_z_value = log(temp_paper_result$outcome)/temp_paper_sd
  } else {
    temp_paper_sd <- (temp_paper_result[,5] - temp_paper_result[,4])/2/qnorm(0.95)
    temp_paper_z_value = temp_paper_result$outcome / temp_paper_sd
  }
  
  temp_comp_trt1 <- temp_paper_result$t1
  temp_comp_trt2 <- temp_paper_result$t2
  number_of_comp <- length(temp_comp_trt1)
  temp_zvalue = numeric(number_of_comp)
  temp_adj_zvalue = numeric(number_of_comp)
  temp_comp_names = character(number_of_comp)

  for(j in 1:number_of_comp){
    temp_comp_names[j] = paste(trt_name[temp_comp_trt1[j]],
                              "vs",
                               trt_name[temp_comp_trt2[j]])
    
    if(temp_comp_trt1[j] == 8){
      temp_map_col <- which(temp_trt_no == temp_comp_trt2[j])
      temp_est = temp_poster[,temp_map_col]
      temp_adj_est = temp_adj_poster[,temp_map_col]
    } else if (temp_comp_trt2[j] == 8){
      temp_map_col <- which(temp_trt_no == temp_comp_trt1[j])
      temp_est = -temp_poster[,temp_map_col]
      temp_adj_est = -temp_adj_poster[,temp_map_col]
    } else {
      temp_map_col1 <- which(temp_trt_no == temp_comp_trt1[j])
      temp_map_col2 <- which(temp_trt_no == temp_comp_trt2[j])
      temp_est = temp_poster[,temp_map_col2] - temp_poster[,temp_map_col1]
      temp_adj_est = temp_poster[,temp_map_col2] - temp_adj_poster[,temp_map_col1]
    }
    temp_zvalue[j] = mean(temp_est)/sd(temp_est)
    temp_adj_zvalue[j] = mean(temp_adj_est)/sd(temp_adj_est)
  }
  
  zvalue_result = cbind(temp_paper_z_value, 
                        temp_zvalue,
                        temp_adj_zvalue)
  par(mar=c(5, 5, 5, 5))
  zvalue_plot(temp_paper_z_value, temp_zvalue, temp_adj_zvalue,outcome_name[i],i)
  
  if(i == 4){
    mtext("Z value obtained by BCL/BCL-OFS methods",side=1,line=-8,outer=TRUE,cex=1, 
          font = 2)
    mtext("Z value obtained by Lu and Ades' approach",side=2,line=-1.5,outer=TRUE,cex=1,
          font=2)
    # plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    # legend("bottom",legend = c("No Adjustment", "OFS Adjustment"), 
    #        pch = c(1, 19),
    #        ncol = 2,
    #        cex = 1.5)
  }
  
  
  rownames(zvalue_result) = temp_comp_names
  write.csv(zvalue_result, paste(outcome_name[i],"_zvalue.csv",sep=""))
 }
par(mar=c(1,1,0.25,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",legend = c("BCL", "BCL-OFS"),
       ncol=2,
       pch = c(1, 19),
       cex = 1.75)
