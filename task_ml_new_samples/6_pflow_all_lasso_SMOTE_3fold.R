##### Run this script from here #####
library(smotefamily)
setwd("/Users/xuhaifeng/Documents/PhD_project_1/scripts/")
source("my_function.r")
initialization()
setwd("/Users/xuhaifeng/Documents/PhD_project_4/data/processed_data_HX/")
old_pflow = readRDS("pflow_for_ml.rds")
setwd("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2023")
pflow = readRDS("baseline_pflow_all_cleaned.rds")
pflow[,3:ncol(pflow)] = apply(pflow[,3:ncol(pflow)], 2, as.numeric)

## classify CR into PR together as responder
pflow$labels[pflow$labels == "responder (CR)"] = "responder"
## remove CR samples
#pflow = pflow[pflow$labels != "responder (CR)",]

###LOOCV
# initialization of the nested cross-validation

#setwd("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2023")
#saveRDS(test_index, "best_test_index.rds")
#test_index = readRDS("best_test_index.rds")

auc_list = list()
prauc_list = list()
ba_list = list()
accuracy_list = list()
GeneList = list()
cutoff_list = list()
f1_list = list()
sens_list = list()
spec_list = list()
set.seed(1)

smote_or_not = FALSE
nfolds_diy = 3
test_index = createFolds(y = pflow$labels, k = nfolds_diy, list = TRUE, returnTrain = FALSE)
par(mfrow=c(nfolds_diy,2))
par(mfrow=c(1,1))
plot_list = list()

for (i in 1:nfolds_diy) {
  test_data <- pflow[as.numeric(test_index[[i]]), ]
  train_data <- pflow[as.numeric(-test_index[[i]]), ]
  #test_data <- randomized_data[as.numeric(test_index[[i]]), ]
  #train_data <- randomized_data[-as.numeric(test_index[[i]]), ]
  x.test = test_data[, 3:ncol(test_data)]
  y.test = test_data$labels
  x.train = train_data[, 3:ncol(train_data)]
  y.train = train_data$labels
  
  y.train = as.character(y.train)
  y.train[y.train == "non_responder"] = 0
  y.train[y.train == "responder"] = 1
  y.train = as.numeric(y.train)
  y.train = factor(y.train, levels = c(0,1))
  y.test = as.character(y.test)
  y.test[y.test == "non_responder"] = 0
  y.test[y.test == "responder"] = 1
  y.test = as.numeric(y.test)
  
  if(smote_or_not == TRUE){
    #SMOTE oversampling
    temp = SMOTE(x.train, target = y.train, K = 7, dup_size = 0)
    temp = temp$data
    y.train = temp$class
    y.train = factor(as.numeric(y.train), levels = c(0,1))
    x.train = temp[,1:ncol(temp)-1]
  }
  
  
  
  pr_label = vector(mode = "logical", length = length(y.test))
  for (pr_i in 1:length(y.test)) {
    if(y.test[pr_i] == 1){
      pr_label[pr_i] = TRUE
    }
  }
  
  if(smote_or_not == TRUE){
    measure_diy = "auc"
  }else{
    measure_diy = "mse"
  }
  
  y.test = as.factor(y.test)
  cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), nfolds = nfolds_diy,
                         alpha = 1, family = "binomial", type.measure=measure_diy)
  #plot(cv_output)
  model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
                alpha = 1, family = "binomial", type.measure=measure_diy)
  
  #get the common genes
  betaGene = model$beta[,model$lambda==cv_output$lambda.min]
  betaGene = betaGene[betaGene != 0]
  betaGene = as.data.frame(betaGene)
  GeneList[[i]] = rownames(betaGene)
  
  
  result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
  label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
  result = as.numeric(result)
  
  #results[i,1] = result
  
  roc = roc(response = y.test, predictor = result, direction = "<", ci=TRUE, plot=FALSE)
  #plot(roc, yaxt = "n",cex.axis=1.5, cex.lab = 1.5, main = "")
  #text(0.5, 0.07, paste("AUC = ",round(auc(roc),3)), cex = 1.5)
  # Draw the y-axis.
  #axis(side = 2, labels = FALSE)
  #axis(side = 2,
  #    ## Rotate the labels.
  #     las = 2,
  #     ## Adjust the label position.
  #     mgp = c(3, 0.75, 0), cex.axis=1.5, cex.lab = 1.5)
  auc_list[[i]] = pROC::auc(roc)
  
  ciobj <- ci.se(roc, specificities=seq(0, 1, l=25))
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  this_auc = paste("AUC =", round(pROC::auc(roc),3))
  
  if(smote_or_not == FALSE){
    the_color = "pink"
  }else{
    the_color = "steelblue"
  }
  
  temp = ggroc(roc) + theme_minimal() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal() + 
    geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = the_color, alpha= 0.2) + ggtitle(capture.output(roc$ci)) +
    annotate("text", x=0.25, y=0.25, label= this_auc) 
  plot_list[[i]] = temp
  
  # ROCR_pred_test <- prediction(result,y.test)
  # ROCR_perf_test <- performance(ROCR_pred_test,'tpr','fpr')
  # #plot(ROCR_perf_test,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
  # cost_perf = performance(ROCR_pred_test, "cost", cost.fp = 1, cost.fn = 2)  
  # the_cutoff = ROCR_pred_test@cutoffs[[1]][which.min(cost_perf@y.values[[1]])]
  # if(the_cutoff == Inf){
  #   cost_perf = performance(ROCR_pred_test, "cost")  
  #   the_cutoff = ROCR_pred_test@cutoffs[[1]][which.min(cost_perf@y.values[[1]])]
  # }
  
  cutoff_label = as.numeric(as.character(y.train))
  the_cutoff = length(cutoff_label[cutoff_label==1])/length(cutoff_label)
  #the_cutoff = 0.5
  
  class_results = result
  class_results[class_results>=the_cutoff] = 1
  class_results[class_results<the_cutoff] = 0
  class_results = as.factor(as.character(class_results))
  labels_all = as.factor(as.character(y.test))
  CM = confusionMatrix(data = class_results, reference = labels_all)
  
  precision =  CM$byClass[6] #this function has some problem so the recall is actually precision
  recall =  CM$byClass[5] #this function has some problem so the recall is actually precision
  f1 = 2*precision*recall/(precision+recall)
  
  f1_list[[i]] = f1
  sens_list[[i]] = CM$byClass[2] #this function has some problem so the sens is actually spec
  spec_list[[i]] = CM$byClass[1] 
  
  accuracy_list[[i]] = CM$overall[1]
  ba_list[[i]] = CM$byClass[[11]]
  cutoff_list[[i]] = the_cutoff
  #CM
  
  #wpr<-pr.curve(result, weights.class0 = pr_label, curve = TRUE,
  #              max.compute = T, min.compute = T, rand.compute = T)
  #prauc_list[[i]] = wpr$auc.integral
  #plot(wpr, add = FALSE, color = "red", yaxt = "n",cex.axis=1.5, cex.lab = 1.5, main = "")
  
  #no_skill = length(y.test[y.test==1]) / length(y.test)
  #lines(x= seq(from=0, to=1, by=0.05), y = rep_len(no_skill, 21),  type = "c", col = "red", lwd = 3)
  #text(0.5, 0.07, paste("baseline = ", round(no_skill,3)), cex = 1.5)
  
  
  # Draw the y-axis.
  #axis(side = 2, labels = FALSE)
  #axis(side = 2,
       ## Rotate the labels.
  #     las = 2,
       ## Adjust the label position.
  #     mgp = c(3, 0.75, 0), cex.axis=1.5, cex.lab = 1.5)
  
  
  #calculate average mse for this fold
  # mse = cbind(as.data.frame(result), as.data.frame(as.numeric(as.character(y.test))))
  # temp_c = as.data.frame(matrix(0,nrow(mse),1))
  # mse = cbind(mse, temp_c)
  # for (j in 1:nrow(mse)) {
  #   mse[j,3] = (mse[j,1] - mse[j,2])^2
  # }
  #mse_results[1,i] = mean(mse[j,3])
}
mean(unlist(auc_list))
mean(unlist(prauc_list))
mean(unlist(accuracy_list))
mean(unlist(ba_list))
mean(unlist(sens_list))
mean(unlist(spec_list))
mean(unlist(f1_list))

GeneList = unlist(GeneList)
temp = as.data.frame(table(GeneList))

#cutoff_list
ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
          #labels = c("Fold 1", "Fold 2", "Fold 3"),
          ncol = 3, nrow = 1)

ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],
          #labels = c("Fold 1", "Fold 2", "Fold 3"),
          ncol = 3, nrow = 2)
