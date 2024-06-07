#### correlation
#install.packages("clValid")
# initialization 
library("clValid")
setwd("/Users/xuhaifeng/Documents/PhD_project_1/scripts/")
source("my_function.r")
initialization()
library(ggfortify)
setwd("/Users/xuhaifeng/Documents/PhD_project_4/projects_2024/dss_2024/data")
data = readRDS("dss_cleaned_data.rds")
df = cbind(as.character(data[[2]]), data[[1]])
library(FeatureImpCluster)
library(flexclust)
df2 = cbind(df$`BGB-11417+Copanlisib`, 
            df$`Copanlisib+Venetoclax`,
            df$`BGB-11417+Duvelisib`,
            df$`Venetoclax+ZSTK474`)
colnames(df2) = c("BGB-11417+Copanlisib", "Copanlisib+Venetoclax",
                  "BGB-11417+Duvelisib", "Venetoclax+ZSTK474")
rownames(df2) = rownames(df)

setwd("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2023")
pflow = readRDS("new_pflow_data_55.rds")
rownames(pflow) = pflow$ids
length(intersect(rownames(pflow), rownames(df2)))
common_ids = intersect(rownames(pflow), rownames(df2))
pflow = pflow[common_ids,]


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

dss_vali = read.csv2("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2024/DSS_JBrown_common_samples_new.csv")
dss_vali = t(dss_vali)
colnames(dss_vali) = dss_vali[1,]
dss_vali = dss_vali[-1,]
dss_vali = as.data.frame(dss_vali)
dss_vali = dss_vali[dss_vali$`screening/responding/progression` == "screening", ]
rownames(dss_vali) = dss_vali$id
dss_vali = dss_vali[,3:ncol(dss_vali)]
dss_vali = apply(dss_vali, 2, as.numeric)

vali_set = vali_set[-3,]
vali_set = vali_set[-7,]
dss_vali = as.data.frame(dss_vali)
rownames(dss_vali) = rownames(vali_set)



dss = df2[common_ids,]
dss = as.data.frame(dss)
dss_vali = dss_vali[, colnames(dss)]
dss = rbind(dss, dss_vali)

pflow_pure = pflow[, colnames(vali_set)]
pflow_pure = rbind(pflow_pure, vali_set)

#pair
k = 1
cor_list = list()
for (i in 1:4) {
  for (j in 1:30) {
    if(cor(dss[,i], pflow_pure[,j], method = "pearson")> 0.5 | cor(dss[,i], pflow_pure[,j], method = "spearman")< -0.5){
      this_pair = c(colnames(dss)[i], colnames(pflow_pure)[j], cor(dss[,i], pflow_pure[,j]))
      cor_list[[k]] = this_pair
      k = k +1
    }
  }
  
}

plot_list = list()

new_labels = rep("Res_vali", 10)
labels = c(pflow$labels, new_labels)
labels[labels == "non_responder"] = "NR (Co. 1)"
labels[labels == "responder"] = "R (Co. 1)"
labels[labels == "Res_vali"] = "R (Co. 2)"


cor_list = list()
p_list = list()
p_list_adjusted = list()
dss_list = list()
pflow_list = list()
k = 1
for (i in 1:4) {
  for (j in 1:30) {
    temp = cor.test(dss[,i], pflow_pure[,j], method = "pearson")
    cor_list[[k]] = temp$estimate
    p_list[[k]] = temp$p.value
    p_list_adjusted[[k]] = temp$p.value
    dss_list[[k]] = colnames(dss)[i]
    pflow_list[[k]] = colnames(pflow_pure)[j]
    #this_pair = c(colnames(dss)[i], colnames(pflow)[j], cor(dss[,i], pflow[,j]))
    #cor_list[[k]] = this_pair
    k = k +1
  }
}

p_list = unlist(p_list)
cor_list = unlist(cor_list)
dss_list = unlist(dss_list)
pflow_list = unlist(pflow_list)
p_list_adjusted = unlist(p_list_adjusted)
p_list_adjusted <- p.adjust(p_list_adjusted, method = "fdr")

stats_these = cbind.data.frame(dss_list, pflow_list, cor_list, p_list, p_list_adjusted)

stats_these = stats_these[stats_these$cor_list>0.5 | stats_these$cor_list< (-0.5), ]

my_round <- function(x, n)  {
  x1 <- round(x, n)
  if(round(x, n) == 0) {
    n1 <- stringr::str_locate(x, "[^-0.]")[1] -  str_locate(x, fixed("."))[1]
    print(n1)
    x1 <- round(x, n1)
    
  }
  return(x1)
  
  
}
options(scipen=999)
plot_list = list()

for (i in 1:6) {
  this_cor = stats_these[i,3]
  this_p = stats_these[i,4]
  this_adjust_p = stats_these[i,5]
  this_dss = stats_these[i,1]
  this_dss_values = dss[, this_dss]
  this_pflow = stats_these[i,2]
  this_pflow_values = pflow_pure[, this_pflow]
  temp = cbind.data.frame(this_dss_values, this_pflow_values, labels)
  
  colnames(temp) = c("DSS", "pFLOW", "Response")
  if(this_pflow == "S6-ribosomal protein (pS235/236)"){
    this_pflow = "S6-rib. prot. (pS235/236)"
  }
  plot_list[[i]] = 
    ggscatter(temp, x= "DSS", y= "pFLOW",
              #facet.by  = "Response", 
              fill  = "Response",
              shape = 21, size = 4,
              add = "reg.line", 
              conf.int = FALSE, 
              cor.coef = FALSE, #set this to true if you want to output R and p
              cor.method = "pearson", 
              title=paste("R","=",round(stats_these$cor_list[i],2),", ",
                          "p=",my_round(stats_these$p_list[i],4), "\n", 
                          "adjusted p=", my_round(stats_these$p_list_adjusted[i],3), sep = ""),
              #title="",
              show.legend.text = FALSE,
              cor.coef.size = 7,
              xlab= paste(this_dss, "(DSS)"),
              ylab = paste(this_pflow,"\n" ,"(Arcsinh ratio)"),
              palette = c("#91C0F9", "#097F98", "#F2B77C")) + 
    theme(axis.text=element_text(size=18),
          legend.position = "none",
          axis.title=element_text(size=18)) +
    font("legend.title", size = 15, face = "bold") +
    font("legend.text",  size = 15) 
}


ggarrange(plotlist = plot_list,
          
          #labels = c("A", "B", "C"),
          ncol = 3, nrow = 2)

# draw non-significant pairs
stats_these = cbind.data.frame(dss_list, pflow_list, cor_list, p_list, p_list_adjusted)
stats_these = stats_these[stats_these$cor_list>0.5 | stats_these$cor_list< (-0.5), ]
these_proteins = unique(stats_these$pflow_list)
stats_these = cbind.data.frame(dss_list, pflow_list, cor_list, p_list, p_list_adjusted)
stats_these = stats_these[is.element(stats_these$pflow_list, these_proteins),]
stats_these$dss_list = gsub("BGB-11417", "Sonrotoclax", stats_these$dss_list)
stats_these$dss_list[stats_these$dss_list == "Copanlisib+Venetoclax"] = "Venetoclax+Copanlisib"
stats_these$pflow_list[stats_these$pflow_list == "S6-ribosomal protein (pS235/236)"] = "S6-rib. prot. (pS235/236)"
colnames(dss) = gsub("BGB-11417", "Sonrotoclax", colnames(dss))
colnames(dss)[colnames(dss) == "Copanlisib+Venetoclax"] = "Venetoclax+Copanlisib"
colnames(pflow_pure)[colnames(pflow_pure) == "S6-ribosomal protein (pS235/236)"] = "S6-rib. prot. (pS235/236)"

library("grDevices")
for (i in 1:12) {
  this_cor = stats_these[i,3]
  this_p = stats_these[i,4]
  this_adjust_p = stats_these[i,5]
  this_dss = stats_these[i,1]
  this_dss_values = dss[, this_dss]
  this_pflow = stats_these[i,2]
  this_pflow_values = pflow_pure[, this_pflow]
  temp = cbind.data.frame(this_dss_values, this_pflow_values, labels)
  
  colnames(temp) = c("DSS", "pFLOW", "Response")
  #if(this_pflow == "S6-ribosomal protein (pS235/236)"){
  #  this_pflow = "S6-rib. prot. (pS235/236)"
  #}
  plot_list[[i]] = 
    ggscatter(temp, x= "DSS", y= "pFLOW",
              #facet.by  = "Response", 
              fill  = "Response",
              shape = 21, size = 4,
              add = "reg.line", 
              conf.int = FALSE, 
              cor.coef = FALSE, #set this to true if you want to output R and p
              cor.coef.size = 7,
              cor.method = "pearson", 
              title=paste("R=",round(stats_these$cor_list[i],2),", ",
                          "p=",my_round(stats_these$p_list[i],4), "\n", 
                          "adjusted p=", my_round(stats_these$p_list_adjusted[i],3), sep = ""),
              #title="",
              show.legend.text = FALSE,
              xlab= paste(this_dss, "(DSS)"),
              ylab = paste(this_pflow,"\n" ,"(Arcsinh ratio)"),
              palette = c("#91C0F9", "#097F98", "#F2B77C")) + 
    theme(axis.text=element_text(size=18),
          legend.position = "none",
          title=element_text(size=15),
          axis.title=element_text(size=18)) +
    font("legend.title", size = 15, face = "bold") +
    font("legend.text",  size = 15) 
  
}
plot_list[1] 
plot_list[2]
plot_list[3] 
plot_list[4]
plot_list[5] 
plot_list[6]
plot_list[7] 
plot_list[8] 
plot_list[9] 
plot_list[10] 
plot_list[11] 
plot_list[12] 
