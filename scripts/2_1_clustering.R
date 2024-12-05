# initialization 
library("clValid")
library("ggfortify")
library("FeatureImpCluster")
library("flexclust")
library("ComplexHeatmap")
source("/scripts/my_function.r")
initialization() #load required packages

# Due to the ethical policy, the training data is not shared here,
# but it can be required from the corresponding author Sigrid S Skåland: sigrid.skanland@ous-research.no
# The data is a list consisting 2 elements, while the first element is the DSS values, and the second it the treatment response
# The "df" object is data frame with 31 samples as rows, including 10 responders, 20 non-responders, and 1 sample with no clear information; 
# and it has 182 colums, including the treatment response as the 1st column, and the DSS of 181 drug/drug combinations as the rest columns.

setwd("/data folder path") # insert your data path here
data = readRDS("dss_cleaned_data.rds")
df = cbind(as.character(data[[2]]), data[[1]])

# select only non-responder data
set.seed(10) # Random seed is VERY important to k-means clustering, do not change here
non_res = df[df$`as.character(data[[2]])`== "SD" | df$`as.character(data[[2]])` == "PD",]
non_res = non_res[, 2:ncol(non_res)]
non_res[1, 1]
res_2 <- flexclust::kcca(non_res,k=2, family=kccaFamily("kmeans"))

set.seed(10) # Random seed is VERY important to k-means clustering, do not change here
FeatureImp_res <- FeatureImpCluster(res_2,as.data.table(non_res), biter = 5)
#View(FeatureImp_res$featureImp)
temp = FeatureImp_res$misClassRate
temp = as.matrix(temp)
temp2 <- temp[,colSums(temp) != 0]
boxplot(temp2)
temp3 = colSums(temp2)/0.05
barplot(temp3)
temp3 = as.matrix(temp3)
temp3 = t(temp3)
colnames(temp3)
colnames(temp3) = c( # add "\n" for better illustrtaion
  "Sonrotoclax+\nCopanlisib",  
  "Copanlisib+\nVenetoclax", 
  "Sonrotoclax+\nDuvelisib",   
  "Venetoclax+\nZSTK474"   
)
barplot(temp3, ylim = c(0,50), col = "grey")

## keep only important features
these_features = FeatureImp_res$featureImp
these_features = cbind(names(these_features), these_features)
these_features = these_features[these_features[,2]!=0,]
these_features = as.character(these_features[,1])
selected_matrix = non_res[,these_features]
colnames(selected_matrix) = c("","","","") #only for better drawing
res_with_selected_features <- kcca(selected_matrix,k=2)
barplot(res_with_selected_features)
coul <- c("red", "blue") 
barplot(res_with_selected_features, col=coul)
selected_matrix = non_res[,these_features] # get original colnames back

clusters_all = res_2@cluster
clusters_select = res_with_selected_features@cluster
clusters_all == clusters_select


# compare between all features, 4 important features, and random 4 features
library("factoextra")
res.km_all <- kmeans(scale(non_res), 2)
# Note that cluster 1 and cluster 2 are just names so they can be opposite in different tests
# the important thing is that the same patients are clustered into the same groups.
# So one can switch them for better illustration, but here it is not switched
#res.km_all$cluster[res.km_all$cluster == 1] = 99 
#res.km_all$cluster[res.km_all$cluster == 2] = 1
#res.km_all$cluster[res.km_all$cluster == 99] = 2

res.km <- kmeans(scale(selected_matrix), 2)
res.km$cluster
#saveRDS(res.km, "clustering_results.rds")

#res.km_random <- kmeans(scale(non_res[,5:8]), 2)
#colnames(non_res[,5:8])
#res.km_random$cluster
#res.km_random$cluster == res.km_all$cluster

random_id = sample(1:181, 4)
res.km_random <- kmeans(scale(non_res[,random_id]), 2)
colnames(non_res[,random_id])
res.km_random$cluster
res.km_random$cluster == res.km_all$cluster

####
cluserting_results = res.km$cluster
cluserting_results = as.character(cluserting_results)
cluserting_results

cluster_1 = cluserting_results
cluster_1[cluster_1 == "1"] = TRUE
cluster_1[cluster_1 == "2"] = FALSE
cluster_1 = as.logical(cluster_1)
cluster_1_n = rownames(selected_matrix)
cluster_1 = cluster_1_n[cluster_1]

cluster_2 = cluserting_results
cluster_2[cluster_2 == "1"] = FALSE
cluster_2[cluster_2 == "2"] = TRUE
cluster_2 = as.logical(cluster_2)
cluster_2_n = rownames(selected_matrix)
cluster_2 = cluster_2_n[cluster_2]

col_fun = colorRamp2(c(0, 49, 98), c("blue", "white", "red"))

selected_matrix = non_res[,these_features] # get original colnames back
selected_matrix = as.matrix(selected_matrix)
selected_matrix = selected_matrix[,1:4]
selected_matrix = apply(selected_matrix, 2, as.numeric)
rownames(selected_matrix) = rownames(non_res)
colnames(selected_matrix) = gsub("BGB-11417", "Sonrotoclax", colnames(selected_matrix))
colnames(selected_matrix)[colnames(selected_matrix) == "Copanlisib+Venetoclax"] = "Venetoclax+Copanlisib"

Heatmap(selected_matrix,
        name = "DSS", 
        #right_annotation = ha, 
        #row_km = 10,
        rect_gp = gpar(col = "white", lwd = 2), # add white borders around cells
        #column_title = "Clustering based on Euclidean distance.",
        column_title_gp = gpar(fontsize = 14),
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        #clustering_method_rows = "kmeans",
        #row_km = 2,
        show_row_names = TRUE,
        row_split = factor(cluserting_results, levels = c("2", "1")),
        row_title = "Cluster %s",
        row_title_gp = gpar(fontsize = 14),
        cluster_row_slices = FALSE,
        #cluster_rows = cluster_within_group(selected_matrix, group),
        # euclidean, maximum, manhattan, canberra, minkowski, pearson, spearman, kendall
        row_names_side = "left",
        row_title_side = "right",
        row_dend_side = "right", 
        #row_dend_reorder = TRUE,
        #row_split = target_matrix,
        #row_title = NULL,
        #column_dend_side = "top", 
        show_row_dend = TRUE,
        show_column_dend = FALSE,
        #show_row_dend = FALSE,
        column_names_gp = gpar(fontsize = 17), row_names_gp = gpar(fontsize = 17), 
        #column_title = "euclidean",
        col = col_fun
        #cell_fun = function(j, i, x, y, width, height, fill) {
        #  if(p_value_matrix[i, j] < 0.05)
        #    #grid.text(sprintf("%.2f", result_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
        #    grid.text(sprintf("X"), x, y, gp = gpar(fontsize = 12))
        #}
)

selected_matrix = t(selected_matrix)
rownames(selected_matrix) = c("Sonrotoclax + Copanlisib",
                              "Venetoclax + Copanlisib", 
                              "Sonrotoclax + Duvelisib",  
                              "Venetoclax + ZSTK474")
Heatmap(selected_matrix,
        name = "DSS", 
        #right_annotation = ha, 
        #row_km = 10,
        rect_gp = gpar(col = "white", lwd = 2), # add white borders around cells
        #column_title = "Clustering based on Euclidean distance.",
        row_title_gp = gpar(fontsize = 14),
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        #clustering_method_rows = "kmeans",
        #row_km = 2,
        show_column_names = TRUE,
        column_split = factor(cluserting_results, levels = c("2", "1")),
        column_title = "Cluster %s",
        column_title_gp = gpar(fontsize = 14),
        cluster_column_slices = FALSE,
        #cluster_rows = cluster_within_group(selected_matrix, group),
        # euclidean, maximum, manhattan, canberra, minkowski, pearson, spearman, kendall
        row_names_side = "left",
        row_title_side = "right",
        row_dend_side = "right", 
        #row_dend_reorder = TRUE,
        #row_split = target_matrix,
        #row_title = NULL,
        #column_dend_side = "top", 
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        #show_row_dend = FALSE,
        column_names_gp = gpar(fontsize = 15), row_names_gp = gpar(fontsize = 15), 
        #column_title = "euclidean",
        col = col_fun
        #cell_fun = function(j, i, x, y, width, height, fill) {
        #  if(p_value_matrix[i, j] < 0.05)
        #    #grid.text(sprintf("%.2f", result_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
        #    grid.text(sprintf("X"), x, y, gp = gpar(fontsize = 12))
        #}
)


selected_matrix = df[,these_features]
res = df[df$`as.character(data[[2]])`!= "SD" & df$`as.character(data[[2]])` != "PD",]
res = res[res$`as.character(data[[2]])`!="GREY", ]
res = res[, these_features]

cluster_1 = c("PHA0203", "FMB0204", "SEA0201", "SEA0202", "HPA0205",
              "FMB0203","PHA0207", "ROA0204", "HPA0206", "SEA0203", "FMB0202")
cluster_2 = c("PHA0212", "HPA0202", "PHA0215", "KNB0203", "PHA0213", 
              "HPA0201", "KNB0201", "PHA0202", "PHA0204")
cluster_2 = selected_matrix[cluster_2,]
cluster_2 = as.data.frame(cluster_2)
cluster_1 = selected_matrix[cluster_1,]
cluster_1 = as.data.frame(cluster_1)
cluster_2 = cbind.data.frame(cluster_2, rep("cluster 2", 9))
cluster_1 = cbind.data.frame(cluster_1, rep("cluster 1", 11))
colnames(cluster_2)[5] = "cluster"
colnames(cluster_1)[5] = "cluster"
res$cluster = "responder"
dat_f = rbind(cluster_1, cluster_2, res)
colnames(dat_f) = c("BGB-11417\n+Copanlisib",  
                    "Copanlisib\n+Venetoclax",
                    "BGB-11417\n+Duvelisib", 
                    "Venetoclax\n+ZSTK474",
                    "cluster")
#write.csv(dat_f, "dss_cluster_boxplot.csv")
dat_m = melt(dat_f)

colnames(dat_m) = c("Cluster", "Drug combination", "DSS")

ggplot(dat_m, aes(x = `Drug combination`, y = DSS, fill = Cluster)) +
  geom_boxplot(outlier.colour="black", outlier.shape=1,
               outlier.size=4)+
  scale_fill_manual("Cluster", values = c("cluster 2" = "blue", "cluster 1" = "red"))+
  # Implement a grouped bar chart
  #geom_bar(position = "dodge", stat = "identity") + 
  theme(
    axis.text=element_text(size=19),
    axis.title=element_text(size=19),
    plot.title = element_text(size = 16),
    legend.title=element_text(size=16), 
    legend.text=element_text(size=16), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    axis.line = element_line(colour = "black"))
#panel.border = element_rect(colour = "black", fill=NA, size=1))


### draw heatmap for the entire data
df = cbind(as.character(data[[2]]), data[[1]])
df = df[df$`as.character(data[[2]])`!="GREY",]
cluster_1 = c("PHA0203", "FMB0204", "SEA0201", "SEA0202", "HPA0205",
              "FMB0203","PHA0207", "ROA0204", "HPA0206", "SEA0203", "FMB0202")
cluster_2 = c("PHA0212", "HPA0202", "PHA0215", "KNB0203", "PHA0213", 
              "HPA0201", "KNB0201", "PHA0202", "PHA0204")
responders = rownames(df)[is.element(rownames(df), c(cluster_1, cluster_2)) == FALSE]

df = df[,2:ncol(df)]
cluster_1 = df[cluster_1, ]
cluster_2 = df[cluster_2, ]
responders = df[responders, ]
df = rbind(cluster_1, cluster_2, responders)
colnames(df) = gsub("BGB-11417", "Sonrotoclax", colnames(df))
colnames(df)[colnames(df) == "Copanlisib+Venetoclax"] = "Venetoclax+Copanlisib"


cluster_labesls = c(rep("Cluster 1", 11), rep("Cluster 2", 9), rep("Responders", 10))
chr_level = c("Cluster 1", "Cluster 2", "Responders")
ha = columnAnnotation(Labels = cluster_labesls, 
                      #labels_gp = gpar(col = "black", fontsize = 22),
                      #at = chr_level,
                      show_annotation_name = FALSE,
                      col = list(Labels = c("Cluster 1" = "red", 
                                            "Cluster 2" = "blue", 
                                            "Responders" = "darkgreen")),
                      annotation_legend_param = list(Labels = list(
                        ncol = 1, 
                        title = "", 
                        #title_position = "topcenter",
                        at = chr_level))
)


max(df)
min(df)
col_fun = colorRamp2(c(0, 49, 98), c("blue", "white", "red"))
df = t(df)

Heatmap(df,
        name = "DSS", 
        #right_annotation = ha, 
        #row_km = 10,
        rect_gp = gpar(col = "white", lwd = 2), # add white borders around cells
        #column_title = "Clustering based on Euclidean distance.",
        row_title_gp = gpar(fontsize = 11),
        column_split = factor(cluster_labesls, levels = c("Cluster 1", "Cluster 2", "Responders")),
        top_annotation = ha, 
        clustering_distance_rows = "euclidean",
        #clustering_method_rows = "kmeans",
        #row_km = 2,
        #column_split = factor(cluserting_results, levels = c("1", "2")),
        #column_title = "Cluster %s",
        column_title_gp = gpar(fontsize = 11),
        cluster_column_slices = FALSE,
        #cluster_rows = cluster_within_group(selected_matrix, group),
        # euclidean, maximum, manhattan, canberra, minkowski, pearson, spearman, kendall
        row_names_side = "left",
        row_title_side = "right",
        row_dend_side = "right", 
        #row_dend_reorder = TRUE,
        #row_split = target_matrix,
        #row_title = NULL,
        #column_dend_side = "top", 
        show_row_dend = TRUE,
        show_column_names = TRUE,
        show_column_dend = TRUE,
        #show_row_dend = FALSE,
        column_names_gp = gpar(fontsize = 11), row_names_gp = gpar(fontsize = 11), 
        #column_title = "euclidean",
        col = col_fun
        #cell_fun = function(j, i, x, y, width, height, fill) {
        #  if(p_value_matrix[i, j] < 0.05)
        #    #grid.text(sprintf("%.2f", result_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
        #    grid.text(sprintf("X"), x, y, gp = gpar(fontsize = 12))
        #}
)
