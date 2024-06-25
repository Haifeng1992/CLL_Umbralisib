##### Run this script from here #####
source("scripts/my_function.r")
initialization() #load required packages

# Due to the ethical policy, the training data is not shared here,
# but it can be required from the corresponding author Sigrid S Skåland: sigrid.skanland@ous-research.no
# The data contains 56 samples as rows, including 15 responders, 40 non-responders, and 1 sample with no clear information; 
# and it has 32 columns, including the IDs, labels (treatment response), and 30 (phospho) proteins

setwd("/data folder path") # insert your data path here
pflow = readRDS("baseline_pflow_all_cleaned.rds") # read training data
pflow[,3:ncol(pflow)] = apply(pflow[,3:ncol(pflow)], 2, as.numeric)

## Preprocessing. This is kept for potential detail-checking
# classify CR into PR together as responder
pflow$labels[pflow$labels == "responder (CR)"] = "responder"
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



## Read independent test set
# Due to the ethical policy, the test data is not shared here,
# but it can be required from the corresponding author Sigrid S Skåland: sigrid.skanland@ous-research.no
# The data contains 12 samples as rows and all of them are responders;
# and it has 30 (phospho) proteins as columns (the same as the training set)

vali_set = readRDS("pflow_jenifer_brown.rds") 

## Read the training labels and data
# Read the training labels
y.train = pflow$labels
y.train = as.character(y.train)
y.train[y.train == "non_responder"] = 0
y.train[y.train == "responder"] = 1
y.train = as.numeric(y.train)
y.train = factor(y.train, levels = c(0,1))
# Read the training data
x.train = pflow[, 3:ncol(pflow)]

## Modelling
# Tuning parameters for Lasso
cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), nfolds = 5,
                       alpha = 1, family = "binomial", type.measure="mse")
plot(cv_output)
#coef(cv_output)
model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
              alpha = 1, family = "binomial", type.measure="mse")

# Get the important proteins
betaGene = model$beta[,model$lambda==cv_output$lambda.min]
betaGene = betaGene[betaGene != 0]
betaGene = as.data.frame(betaGene)
betaGene # Here shows that only the five proteins got from LOOCV were kept in the model 

result = predict(model, vali_set, type = "response", s=cv_output$lambda.min)
#label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
result = as.numeric(result)
result

## read training data prediction from LOOCV
tg_results = readRDS("final_prediction_results.rds")
results_tg_res = tg_results[1:15]
results_tg_nres = tg_results[16:55]

## Perform Wilcoxon test between each group
wilcox.test(result, results_tg_res)
wilcox.test(result, results_tg_nres)
wilcox.test(results_tg_res, results_tg_nres)

r_list = list()
r_list[[3]] = result
r_list[[2]] = results_tg_res
r_list[[1]] = results_tg_nres
names(r_list) = c("Non Responder (training)", "Responder (training)", "Responder (validation)")

## Draw boxplots
# Note the final version is not created here but with GraphPad Prism 9 for better illustrations 
temp = plyr::ldply(r_list, cbind)
colnames(temp) = c("cohort", "value")
ggplot(temp, aes(x=cohort, y=value, color = cohort)) + 
  geom_violin() + 
  geom_jitter(size = 2, aes(colour = cohort)) + 
  scale_colour_manual(values = c("#228B22" ,"#007dff", "purple"), labels=c("Non-responder (training)","Responder (training)","Responder (test)")) +
  #labs(x = NULL) +
  theme(text=element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.title = element_blank(), 
        #legend.text = element_text(size=18), 
        legend.position = "none",
        axis.text = element_text(size = 20),
        #axis.text.x=element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = "", y="predicted probabilities") +
  scale_x_discrete(labels=c("NR\n(training)","R\n(training)","R\n(test)")) +
  geom_hline(yintercept = 0.273, lty="dashed") 

