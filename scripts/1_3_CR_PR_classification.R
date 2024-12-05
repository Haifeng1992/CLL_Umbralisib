##### Run this script from here #####
setwd("/data folder path") # insert your data path here
source("my_function.r")
initialization()
pflow = readRDS("baseline_pflow_all_cleaned.rds")
pflow[,3:ncol(pflow)] = apply(pflow[,3:ncol(pflow)], 2, as.numeric)

## classify CR into PR together as responder
## these CR were originally labelled as CR but the 
## final clinical report showed they are PR
pflow$labels[pflow$labels == "responder (CR)"] = "responder"

## "PHA0210" and "NYA0206" are found as responder in literature
new_res = c("PHA0210", "NYA0206")
res = pflow[pflow$labels == "responder",]
old_res = res$ids
res = c(old_res, new_res)
res = pflow[is.element(pflow$ids, res),]
res$labels = "responder"

new_non_responder = pflow$ids[is.element(pflow$ids, res$ids) == FALSE]
non_res = pflow[is.element(pflow$ids, new_non_responder),]
non_res$labels = "non_responder"
non_res = non_res[non_res$ids!="NYA0204",]
pflow = rbind(res, non_res)

vali_set = read.csv2("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2023/validation/pflow_jenifer_brown.csv")
cols = c(2,4,6,8,10,12,14,16,18,20,22,24)
rownames(vali_set) = vali_set$X
vali_set = vali_set[, cols]
vali_set = vali_set[1:32, ] # remove the unrelated protein according to Yanping's answer

vali_set = t(vali_set)
vali_set = vali_set[,2:ncol(vali_set)]
is.element(colnames(vali_set), colnames(pflow))
vali_set = vali_set[,c(1,3:31)]
is.element(colnames(vali_set), colnames(pflow))


proteins = c("MEK1 (pS298)",         
             "p38 MAPK (pT180/182)",
             "STAT1 (pS727)",
             "p90RSK (pS380)",      
             "Bcl-2 (pS70)")

PR_train = pflow[pflow$labels == "responder", ]
PR_train_sub = cbind(PR_train$labels, PR_train[,proteins])
rownames(PR_train_sub) = PR_train$ids
colnames(PR_train_sub)[1] = "label"

vali_set_sub = as.data.frame(vali_set)
vali_sub_labels = c("PR", "CR", "CR", "CR", "CR", "CR", 
                    "PR", "CR", "CR", "PR", "PR", "CR")
vali_all_sub = cbind(vali_sub_labels, vali_set_sub)
vali_all_sub_PR = vali_all_sub[vali_all_sub$vali_sub_labels == "PR",]
vali_all_sub_CR = vali_all_sub[vali_all_sub$vali_sub_labels == "CR",]

colnames(vali_all_sub_PR)[1] = "label"
vali_all_sub_PR = cbind(vali_all_sub_PR$label, vali_all_sub_PR[,proteins])
colnames(vali_all_sub_PR)[1] = "label"
new_train = rbind(PR_train_sub, vali_all_sub_PR)
new_train$label[new_train$label == "responder"] = "PR"

colnames(vali_all_sub_CR)[1] = "label"
vali_all_sub_CR = cbind(vali_all_sub_CR$label, vali_all_sub_CR[,proteins])
colnames(vali_all_sub_CR)[1] = "label"
new_train = rbind(new_train, vali_all_sub_CR)

# Mode 1: combine all responder samples into one data set, then use LOOCV
# Mode 2: LOOCV on Cohort 2, but each sub-training set is 
# the rest 11 samples + all responders in cohort 1
# change this to 2 if you want to only test on Cohort 2 (default: 1)
mode = 1 
if(mode == 2){
  cohort_1 = new_train[1:15, ]
  new_train = new_train[16:27, ]
}


results = matrix(0,nrow(new_train), 1)
rownames(new_train) = c(1:nrow(new_train))
set.seed(1234)
par(mfrow=c(1,1))
for (i in 1:nrow(new_train)) {
  test_data <- new_train[i, ]
  train_data <- new_train[-i, ]
  if(mode == 2){
    train_data = rbind(cohort_1, train_data)
  }
  
  x.test = test_data[, 2:ncol(test_data)]
  y.test = test_data$label
  x.train = train_data[, 2:ncol(train_data)]
  y.train = train_data$label
  
  y.train = as.character(y.train)
  y.train[y.train == "PR"] = 0
  y.train[y.train == "CR"] = 1
  y.train = as.numeric(y.train)
  y.train = factor(y.train, levels = c(0,1), ordered = TRUE)
  y.test = as.character(y.test)
  y.test[y.test == "PR"] = 0
  y.test[y.test == "CR"] = 1
  y.test = as.numeric(y.test)
  
  pr_label = vector(mode = "logical", length = length(y.test))
  for (pr_i in 1:length(y.test)) {
    if(y.test[pr_i] == 1){
      pr_label[pr_i] = TRUE
    }
  }
  
  y.test = as.factor(y.test)
  
  # keep the original Lasso model as potential comparison
  # Note the lasso results are recorded with the object "result",
  # and not used in the final evaluation
  cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), nfolds = 3,
                         alpha = 1, family = "binomial", type.measure="mse")

  model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
                alpha = 1, family = "binomial", type.measure="mse")
  
  result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
  label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
  result = as.numeric(result)
  
  
  # create SVM model here
  # SVM results are recorded with the object "sam_result", 
  # and stored in the entire "results" list for final evaluation
  y.train = factor(y.train, levels = c(0,1))
  model = svm(x.train, y.train, probability = TRUE)
  sam_result = predict(model,x.test, probability = TRUE)
  sam_result = attr(sam_result,"probabilities")
  sam_result = sam_result[1]
  
  results[i] = sam_result
  
}
labels_all = new_train$label
labels_all[labels_all == "PR"] = 0
labels_all[labels_all == "CR"] = 1
labels_all = as.numeric(labels_all)
results = as.numeric(results)
roc = roc(response = labels_all, predictor = results, direction = "<",ci=TRUE)
plot(roc)
auc(roc)

# create a random roc and use Delong's test to compare with the actual ROC 
random = runif(n=27, min=1e-12, max=.9999999999)
roc2 = roc(response = labels_all, predictor = random, direction = "<",ci=TRUE)
plot(roc2)
auc(roc2)
roc.test(roc, roc2)

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


library("yogiroc")


class_results = results
frac = 1-(8/27)
class_results[class_results>=frac] = 1
class_results[class_results<frac] = 0
class_results = factor(as.character(class_results), levels = c(0,1), ordered = TRUE)
labels_all = factor(as.character(labels_all), levels = c(0,1), ordered = TRUE)
CM = confusionMatrix(data = class_results, reference = labels_all, positive = "1")
CM
#
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
no_skill = 8/27
a = seq(from=0, to=1, by=0.025)
b = rep_len(no_skill, 41)
lines(x=a , y = b,  
      type = "l", col = "red", lty=2, lwd=1)

text(0.6, 0.16, paste("Random AU-PRC = ", round(no_skill,2)), cex = 1.8)



