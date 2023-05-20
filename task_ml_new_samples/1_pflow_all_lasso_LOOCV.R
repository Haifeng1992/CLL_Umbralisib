##### Run this script from here #####
library(smotefamily)
setwd("/Users/xuhaifeng/Documents/PhD_project_1/scripts/")
source("my_function.r")
initialization()
setwd("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2023")
pflow = readRDS("baseline_pflow_all_cleaned.rds")
pflow[,3:ncol(pflow)] = apply(pflow[,3:ncol(pflow)], 2, as.numeric)

## classify CR into PR together as responder
pflow$labels[pflow$labels == "responder (CR)"] = "responder"

###LOOCV
# initialization of the nested cross-validation
results = matrix(0,nrow(pflow), 1)
GeneList = list()

set.seed(1)
par(mfrow=c(1,1))

# run this line if you want to randomlize the labels
#pflow$labels = sample(pflow$labels)

for (i in 1:nrow(pflow)) {
  test_data <- pflow[i, ]
  train_data <- pflow[-i, ]
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
  
  pr_label = vector(mode = "logical", length = length(y.test))
  for (pr_i in 1:length(y.test)) {
    if(y.test[pr_i] == 1){
      pr_label[pr_i] = TRUE
    }
  }
  
  y.test = as.factor(y.test)
  cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), nfolds = 5,
                         alpha = 1, family = "binomial", type.measure="mse")
  #plot(cv_output)
  #coef(cv_output)
  model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
                alpha = 1, family = "binomial", type.measure="mse")
  
  #get the important genes
  betaGene = model$beta[,model$lambda==cv_output$lambda.min]
  betaGene = betaGene[betaGene != 0]
  betaGene = as.data.frame(betaGene)
  GeneList[[i]] = betaGene
  GeneList[[i]] = rownames(betaGene)
  
  result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
  label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
  result = as.numeric(result)
  
  results[i,1] = result
}

labels_all = pflow$labels
labels_all[labels_all == "non_responder"] = 0
labels_all[labels_all == "responder"] = 1
labels_all = as.numeric(labels_all)
results = as.numeric(results)

roc = roc(response = labels_all, predictor = results, direction = "<",ci=TRUE)
pROC::auc(roc)
plot(roc, yaxt = "n",cex.axis=1.5, cex.lab = 1.5, main = "")
text(0.5, 0.07, paste("AUC = ",round(pROC::auc(roc),3)), cex = 1.5)
# Draw the y-axis.
axis(side = 2, labels = FALSE)
axis(side = 2,
     ## Rotate the labels.
     las = 2,
     ## Adjust the label position.
     mgp = c(3, 0.75, 0), cex.axis=1.5, cex.lab = 1.5)

roc_random = roc(response = labels_all, predictor = results, direction = "<")
roc_real = roc(response = labels_all, predictor = results, direction = "<")
roc.test(roc_real,roc_random,method=c("delong"))

ciobj <- ci.se(roc, specificities=seq(0, 1, l=25))
dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])
this_auc = paste("AUC =", round(pROC::auc(roc),3))
ggroc(roc) + theme_minimal() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal() + 
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "lightgreen", alpha= 0.2) + ggtitle(capture.output(roc$ci))+
  annotate("text", x=0.25, y=0.25, label= this_auc) 

class_results = results
class_results[class_results>=0.232] = 1
class_results[class_results<0.232] = 0
class_results = as.factor(as.character(class_results))
labels_all = as.factor(as.character(labels_all))
CM = confusionMatrix(data = class_results, reference = labels_all, positive = "1")
CM

class_results = as.numeric(as.character(class_results))

pr_label = vector(mode = "logical", length = length(labels_all))
for (pr_i in 1:length(labels_all)) {
  if(labels_all[pr_i] == 1){
    pr_label[pr_i] = TRUE
  }
}
wpr<-pr.curve(results, weights.class0 = pr_label, curve = TRUE,
              max.compute = T, min.compute = T, rand.compute = T)
#plot(wpr)
plot(wpr, add = FALSE, color = "lightgreen", yaxt = "n",cex.axis=1.5, cex.lab = 1.5, main = "")

au_prc = wpr$auc.integral
no_skill = length(labels_all[labels_all==1]) / length(labels_all)
lines(x= seq(from=0, to=1, by=0.05), y = rep_len(no_skill, 21),  type = "c", col = "red", lwd = 3)
text(0.5, 0.07, paste("baseline = ", round(no_skill,3)))

# Draw the y-axis.
axis(side = 2, labels = FALSE)
axis(side = 2,
     ## Rotate the labels.
     las = 2,
     ## Adjust the label position.
     mgp = c(3, 0.75, 0), cex.axis=1.5, cex.lab = 1.5)

GeneList = unlist(GeneList)
GeneList = as.data.frame(table(GeneList))
genes = as.character(GeneList[GeneList$Freq>10,1])
