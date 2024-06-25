##### Run this script from here #####
library(smotefamily)
#setwd("/Users/xuhaifeng/Documents/PhD_project_1/scripts/")
#source("my_function.r")
setwd("final/1_1_pflow_lasso_LOOCV.R")
source(".../CLL_Umbralisib/final/my_function.r")
initialization() #load required packages
setwd("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2023")
pflow = readRDS("baseline_pflow_all_cleaned.rds")
pflow[,3:ncol(pflow)] = apply(pflow[,3:ncol(pflow)], 2, as.numeric)

## classify CR into PR together as responder
pflow$labels[pflow$labels == "responder (CR)"] = "responder"

# add new responder based on literature
new_res = c("PHA0210", "NYA0206")
res = pflow[pflow$labels == "responder",]
old_res = res$ids
res = c(old_res, new_res)
#res = old_res
res = pflow[is.element(pflow$ids, res),]
res$labels = "responder"

new_non_responder = pflow$ids[is.element(pflow$ids, res$ids) == FALSE]
non_res = pflow[is.element(pflow$ids, new_non_responder),]
non_res$labels = "non_responder"
non_res = non_res[non_res$ids!="NYA0204",] #remove the sample with no clear info
pflow = rbind(res, non_res)

###LOOCV
# initialization of the nested cross-validation
results = matrix(0,nrow(pflow), 1)
ProteinList = list()
rownames(pflow) = c(1:55)

set.seed(100)
par(mfrow=c(1,1))

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
  y.train = factor(y.train, levels = c(0,1), ordered = TRUE)
  y.test = as.character(y.test)
  y.test[y.test == "non_responder"] = 0
  y.test[y.test == "responder"] = 1
  y.test = as.numeric(y.test)
  
  #temp = SMOTE(x.train, target = y.train, K = 5, dup_size = 0)
  #temp = temp$data
  
  #y.train = temp$class
  #y.train = factor(as.numeric(y.train), levels = c(0,1))
  #x.train = temp[,1:ncol(temp)-1]
  
  pr_label = vector(mode = "logical", length = length(y.test))
  for (pr_i in 1:length(y.test)) {
    if(y.test[pr_i] == 1){
      pr_label[pr_i] = TRUE
    }
  }
  
  y.test = as.factor(y.test)
  cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), nfolds = 3,
                         alpha = 1, family = "binomial", type.measure="mse")
  #plot(cv_output)
  #coef(cv_output)
  model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
                alpha = 1, family = "binomial", type.measure="mse")
                #, lambda=cv_output$lambda.min)
  
  #get the common genes
  betaGene = model$beta[,model$lambda==cv_output$lambda.min]
  betaGene = betaGene[betaGene != 0]
  betaGene = as.data.frame(betaGene)
  ProteinList[[i]] = betaGene
  ProteinList[[i]] = rownames(betaGene)
  
  result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
  label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
  result = as.numeric(result)
  
  model = svm(x.train, y.train, probability = TRUE)
  sam_result = predict(model,x.test, probability = TRUE)
  sam_result = attr(sam_result,"probabilities")
  sam_result = sam_result[1]
  
  results[i,1] = result

}

ProteinList = unlist(ProteinList)
ProteinList = as.data.frame(table(ProteinList))
rownames(ProteinList) = ProteinList$ProteinList
ProteinList = ProteinList[order(ProteinList$Freq, decreasing= TRUE), ]
#barplot(ProteinList$Freq, names.arg = ProteinList$ProteinList, las=2)
proteins = as.character(ProteinList[ProteinList$Freq>10,1])
proteins

labels_all = pflow$labels
labels_all[labels_all == "non_responder"] = 0
labels_all[labels_all == "responder"] = 1
labels_all = as.numeric(labels_all)
results = as.numeric(results)

roc = roc(response = labels_all, predictor = results, direction = "<",ci=TRUE)
library("verification")
roc.area(labels_all, results)
ciobj <- ci.se(roc, specificities=seq(0, 1, l=25))
dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])
this_auc = paste("AU-ROC =", format(round(pROC::auc(roc),2),nsmall = 2))
#this_auc = paste("AU-ROC =", round(pROC::auc(roc),2))
p = ggroc(roc, colour = "black") + 
  scale_x_continuous(breaks = seq(0, 1, 0.2), trans = "reverse", name = "1 - Specificity") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_minimal() + 
  geom_abline(slope=1, intercept = 1, linetype = "longdash", alpha=0.7, color = "red") + 
  coord_equal() + 
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "lightgreen", alpha= 0.2) + 
  ggtitle(paste("",capture.output(roc$ci)))+
  #annotate("text", x=0.25, y=0.25, label= this_auc, size = 8) +
  #xlab("Specificity") + 
  ylab("Sensitivity")

p 
p = p+ theme(axis.text=element_text(size=23),
       axis.title=element_text(size=21),
       plot.title = element_text(size = 21),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(),
       panel.background = element_blank(),
       axis.text.x = element_text(color="black"),
       axis.title.x = element_text(colour = "black"),
       axis.ticks = element_line(color = "black"),
       axis.text.y = element_text(color="black"),
       axis.title.y = element_text(colour = "black"),
       #axis.ticks = 
       #panel.background = element_rect(fill='transparent'),
       #plot.background = element_rect(fill='transparent', color=NA),
       axis.line = element_line(colour = "black"))
       #panel.border = element_rect(colour = "black", fill=NA, size=1))
p


## cutoff based metrics
class_results = results
frac = nrow(res)/nrow(pflow)
class_results[class_results>=frac] = 1
class_results[class_results<frac] = 0
class_results = factor(as.character(class_results), levels = c(0,1), ordered = TRUE)
labels_all = factor(as.character(labels_all), levels = c(0,1), ordered = TRUE)
CM = confusionMatrix(data = class_results, reference = labels_all, positive = "1")
CM

### PR curves
library("yogiroc")
results = results

pr_label = vector(mode = "logical", length = length(labels_all))
for (pr_i in 1:length(labels_all)) {
  if(labels_all[pr_i] == 1){
    pr_label[pr_i] = TRUE
  }
}

TrueData = cbind.data.frame(
  pr_label, results, results
)

#create yogiroc2 object


source("~/Documents/PhD_project_4/scripts/temp_yogi.r")
predictions_tg = as.matrix(results)
colnames(predictions_tg) = "prediction"
yrobj <- yr2(truth=TrueData$pr_label, scores=predictions_tg,names = "")
auprc.CI <- function(yr2,monotonized=FALSE) {
  do.call(rbind,lapply(1:length(yr2),function(i) {
    precCI <- prcCI(yr2[[i]][,"tp"],yr2[[i]][,"tp"]+yr2[[i]][,"fp"])
    auprcs <- apply(precCI,2,function(ppv) {
      if (monotonized) {
        ppv <- monotonize(ppv)
      }
      calc.auc(yr2[[i]][,"tpr.sens"],ppv)
    })
  }))
}
auprc.CI(yrobj)
auprc.pvrandom(yrobj)

# edit the original function for similar plot styles to ggplot
draw.prc.CI_edited <- function(yr2,col=seq_along(yr2),lty=1,
    monotonized=TRUE,balanced=FALSE,legend="bottomleft",
    sampling=c("accurate","quickDirty"),nsamples=1000L,monotonizedSampling=FALSE,
    ...) {

  stopifnot(inherits(yr2,"yr2"))
  sr <- switch(match.arg(sampling,c("accurate","quickDirty")),
    quickDirty=sampleRatesQD,accurate=sampleRates
  )
  if (length(lty) < length(yr2)) {
    lty <- rep(lty,length(yr2))
  }
  # mon <- function(xs) if (monotonized) monotonize(xs) else xs
  ppv <- function(i) configure.prec(yr2[[i]],monotonized,balanced)
  plabel <- ifelse(balanced,"Balanced precision","Precision")
  plot(
    yr2[[1]][,"tpr.sens"],ppv(1),
    type="l",
    xlab="Recall", ylab=plabel,
    xlim=c(0,1),ylim=c(0,1),col="black", lty=lty[[1]], ...
  )
  if(length(yr2) > 1) {
    for (i in 2:length(yr2)) {
      lines(
        yr2[[i]][,"tpr.sens"],ppv(i),
        col=col[[i]], lty=lty[[i]], ...
      )
    }
  } 
  for (i in 1:length(yr2)) {
    # x <- yr2[[i]][,"tpr.sens"]
    prior <- yr2[[i]][1,"tp"]/(yr2[[i]][1,"tp"]+yr2[[i]][1,"fp"])
    # precCI <- prcCI(yr2[[i]][,"tp"],yr2[[i]][,"tp"]+yr2[[i]][,"fp"])
    precCI <- inferPRCCI(samplePRCs(yr2[[i]],N=nsamples,monotonized=monotonizedSampling,sr=sr))
    precCI[,-1] <- apply(precCI[,-1],2,function(column) {
      if (balanced) {
        column <- balance.prec(column,prior)
      } 
      # if (monotonized) {
      #   column <- monotonize(column)
      # }
      column
    })
    polygon(c(precCI[,1],rev(precCI[,1])),
            c(precCI[,"0.025"],rev(precCI[,"0.975"])),
            col=yogitools::colAlpha(col[[i]],0.1),border=NA
    )
  }
  if (!is.na(legend)) {
    legend(legend,sprintf("%s (AUPRC=%.02f;R90P=%.02f)",
        names(yr2),auprc(yr2,monotonized,balanced),recall.at.prec(yr2,0.9,monotonized,balanced)
    ),col=col,lty=lty)
  }
}

wpr<-pr.curve(results, weights.class0 = pr_label, curve = TRUE,
              max.compute = T, min.compute = T, rand.compute = T)
#autoplot(wpr)
#plot(wpr)

au_prc = wpr$auc.integral
#draw PRC curves with confidence itervals,
draw.prc.CI_edited(yrobj,col=c("green"), 
                   legend = NA, 
                   axes=F, balanced=FALSE, 
                   monotonized = FALSE, cex.lab = 1.5)

box(bty="l")
axis(side = 1.8, labels = FALSE, tick = FALSE)
axis(side = 2,
     tick = FALSE,
     ## Rotate the labels.
     las = 2,
     ## Adjust the label position.
     #yaxt = "n",
     mgp = c(3, 0.75, 0), cex.axis=1.8, cex.lab = 1.8)
axis(1, cex.axis = 1.8, tick = FALSE)
no_skill = 0.273 # this is the proportion of positive class in the training set
a = seq(from=0, to=1, by=0.025)
b = rep_len(no_skill, 41)
lines(x=a , y = b,  
      type = "l", col = "red", lty=2, lwd=1)

text(0.6, 0.16, paste("Random AU-PRC = ", round(no_skill,2)), cex = 1.8)
text(0.3, 0.45, paste("AU-PRC = ", round(au_prc,2)), cex = 1.8)



