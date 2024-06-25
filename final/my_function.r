my_zsore <- function(data_column) {
  data_column = (data_column-mean(data_column))/sd(data_column)
  return(data_column)
}



## preselect gene expression explaining 50% variations over all samples
low_expressed_gene_filter = function(log_data, percentage){
  var.x1 <- apply(t(log_data), 2, var)
  var.sort <- sort(var.x1, decreasing=TRUE)
  sum.a <- cumsum(var.sort)
  half.a <- var.sort[which(sum.a>(sum.a[length(var.x1)]/percentage))[1]]
  x1.cov <- t(log_data)[,var.x1>=half.a]
  x1.cov = log2(x1.cov+1)
  return(x1.cov)
}

initialization <- function() {
  library("glmnet")
  #library("pcr")
  library("RColorBrewer")
  library("ComplexHeatmap")
  library("circlize")
  library("tidyverse")
  #library("hrbrthemes")
  library("viridis")
  library("stats")
  library("caret")
  library("pROC")
  #library("DMwR")
  #library("PerfMeas")
  #library("qpgraph")
  library("PRROC")
  library("e1071")
  library("ROCR")
  library("unikn")
  library("animation")
  library("MASS")
  library("ggplot2")
  library("reshape2")
  library(caret)
  library(ggfortify)
  library(gridExtra)
  library("dplyr")
  library("randomForest")
  library("edgeR")
  library("readxl")
  library("ggpubr")
}