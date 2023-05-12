data = read.csv2("/Users/xuhaifeng/Documents/PhD_project_4/results/NEW_SAMPLES_2023/results_for_plot.csv")

data = data.frame(t(data))
data = data[2:nrow(data),]
data = apply(data, 2, as.numeric)
rownames(data) = c("Leave_one_out_CV", "3_fold_CV", "5_fold_CV", "5_fold_CV_and_SMOTE")
colnames(data) = c("AUROC", "AUPRC","Accuracy","Balanced_accuracy", "Sensitivity", "Specificity", "F1_score")

library(data.table)
array2dataTable <- function(x) {
  
  # if is matrix, add third dimension (as.data.table does not melt matrices)
  x.is.matrix <- FALSE
  if (length(dim(x))==2) {
    x.is.matrix <- TRUE
    cat("\nNote: x is a matrix, converting it to array with 3rd dim==3 ..")
    dim(x) <- c(dim(x), 1L)
  }
  # add dimnames
  if (is.null(dimnames(x))) {
    cat("\nNote: Array has no dimnames, using seq of integers ..\n")
    dimnames(x) <- lapply(dim(x), function(X) as.character(seq.int(1, X)))
  }
  DT <- as.data.table(x, na.rm = TRUE)
  if (x.is.matrix==TRUE) DT[,V3:=NULL] # remove third column if converting from 2D matrix
  print(str(DT))
  return(DT)
}

#set baseline value for comparison
baseline = c(0.5,0.232, 0.5, 0.5, 0.5, 0.5, 0.5)
data = rbind(data, baseline)

data = as.matrix(data)
temp = array2dataTable(data)
temp$V1[temp$V1 == 1] = "Leave_one_out_CV"
temp$V1[temp$V1 == 2] = "3_fold_CV"
temp$V1[temp$V1 == 3] = "5_fold_CV"
temp$V1[temp$V1 == 4] = "5_fold_CV_and_SMOTE"
temp$V1[temp$V1 == 5] = "baseline"
temp$V1 = factor(temp$V1, levels = c("3_fold_CV", "5_fold_CV", "5_fold_CV_and_SMOTE", "Leave_one_out_CV", "baseline"))

temp$V2[temp$V2 == 1] = "AUROC"
temp$V2[temp$V2 == 2] = "AUPRC"
temp$V2[temp$V2 == 3] = "Accuracy"
temp$V2[temp$V2 == 4] = "Balanced_accuracy"
temp$V2[temp$V2 == 5] = "Sensitivity"
temp$V2[temp$V2 == 6] = "Specificity"
temp$V2[temp$V2 == 7] = "F1_score"
temp$V2 = factor(temp$V2, levels = c("AUROC", "AUPRC", "Accuracy", "Balanced_accuracy",
                                     "Sensitivity", "Specificity", "F1_score"))

colnames(temp) = c("Method", "Metric", "Value")


ggplot(temp, aes(fill=Method, y=Value, x=Metric)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(labels=c("3-fold CV", "5-fold CV", "5-fold with SMOTE", "Leave-one-out CV", "baseline"), values = c("#F8766D", "#00BA38", "#619CFF", "#FFF173", "grey")) + 
  scale_x_discrete(labels=c("AUROC", "AUPRC", "Accuracy", "Balanced Accuracy",
                            "Sensitivity", "Specificity", "F1 score")) +
  theme(panel.background = element_blank(), axis.title = element_blank()) +  ylim(0, 1)
  
  
setwd("/Users/xuhaifeng/Documents/PhD_project_4/results/NEW_SAMPLES_2023/")
data2 = read.csv2("results_for_plot cutoff_comparison.csv")
data2[,2:ncol(data2)] = apply(data2[,2:ncol(data2)], 2, as.numeric)
rownames(data2) = data2[,1]
data2 = data2[,2:ncol(data2)]
data2 = as.matrix(data2)
colnames(data2) = c("3_fold_Proportion", "3_fold_0_5", "5_fold_Proportion", "5_fold_0_5", "5_fold_SMOTE_Proportion", "5_fold_SMOTE_0_5", "LOOCV_Proportion", "LOOCV_0_5")

temp = array2dataTable(data2)
temp$V1[temp$V1 == 1] = "Accuracy"
temp$V1[temp$V1 == 2] = "Balanced_accuracy"
temp$V1[temp$V1 == 3] = "Sensitivity"
temp$V1[temp$V1 == 4] = "Specificity"
temp$V1 = factor(temp$V1, levels = c("Accuracy", "Balanced_accuracy", "Sensitivity", "Specificity"))
temp$V2[temp$V2 == 1] = "3_fold_Proportion"
temp$V2[temp$V2 == 2] = "3_fold_0_5"
temp$V2[temp$V2 == 3] = "5_fold_Proportion"
temp$V2[temp$V2 == 4] = "5_fold_0_5"
temp$V2[temp$V2 == 5] = "5_fold_SMOTE_Proportion"
temp$V2[temp$V2 == 6] = "5_fold_SMOTE_0_5"
temp$V2[temp$V2 == 7] = "LOOCV_Proportion"
temp$V2[temp$V2 == 8] = "LOOCV_0_5"
temp$V2 = factor(temp$V2, levels = c("3_fold_Proportion", "3_fold_0_5", "5_fold_Proportion", "5_fold_0_5",
                                     "5_fold_SMOTE_Proportion", "5_fold_SMOTE_0_5", 
                                     "LOOCV_Proportion", "LOOCV_0_5"))
colnames(temp) = c("metric", "Method", "value")

ggplot(temp, aes(fill=Method, y=value, x=metric)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(labels=c("3-fold (Pos. proportion)", "3-fold (0.5)", 
                             "5-fold (Pos. proportion)", "5-fold (0.5)", 
                             "5-fold SMOTE (Pos. proportion)", "5-fold SMOTE (0.5)",
                             "LOOCV (Pos. proportion)", "LOOCV (0.5)"), 
                    values = c("#F8766D","darkred", "#00BA38","darkgreen", "#619CFF","darkblue", "#FFF173", "orange")) + 
  scale_x_discrete(labels=c("Accuracy", "Balanced Accuracy", "Sensitivity", "Specificity")) +
  theme(panel.background = element_blank(), axis.title = element_blank()) +  ylim(0, 1)

  
five_fold = data2[,4:5]
five_fold = as.matrix(five_fold)
colnames(five_fold) = c("Proportion", "0_5")
temp = array2dataTable(five_fold)
temp$V1[temp$V1 == 1] = "Accuracy"
temp$V1[temp$V1 == 2] = "Balanced_accuracy"
temp$V1[temp$V1 == 3] = "Sensitivity"
temp$V1[temp$V1 == 4] = "Specificity"
temp$V1 = factor(temp$V1, levels = c("Accuracy", "Balanced_accuracy", "Sensitivity", "Specificity"))
temp$V2[temp$V2 == 1] = "proportion"
temp$V2[temp$V2 == 2] = "0_5"
temp$V2 = factor(temp$V2, levels = c("proportion", "0_5"))
colnames(temp) = c("metric", "cutoff", "value")
five_fold_plot = ggplot(temp, aes(fill=cutoff, y=value, x=metric)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(labels=c("cutoff = Pos. proportion", "cutoff = 0.5"), values = c("#F8766D", "#619CFF")) + 
  scale_x_discrete(labels=c("Accuracy", "Balanced Accuracy", "Sensitivity", "Specificity")) +
  theme(panel.background = element_blank(), axis.title = element_blank()) +  ylim(0, 1)

three_fold_plot + five_fold_plot


dat1 <- data.frame(animal = c("bear", "lion", "tiger"),
                   percentage = c(25, 87, 14))
dat2 <- data.frame(shape = c("circle", "square"),
                   percentage = c(17, 67))
dat3 <- data.frame(color = c("red", "blue", "green"),
                   percentage = c(48, 5, 11))

# Combining all the data
all <- list(dat1, dat2, dat3)
all <- lapply(all, function(dat) {
  dat$type <- colnames(dat)[1]
  colnames(dat)[1] <- "variable"
  dat
})
all <- do.call(rbind, all)
