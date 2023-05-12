pflow_new = readxl::read_xlsx("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2023/baseline_pflow_all.xlsx")
labels = pflow_new[3,2:57]
labels[is.na(labels)] = "non_responder"
labels = as.character(labels)
ids = pflow_new[2,2:57]
ids = stringr::str_replace_all(ids, "\\-", "")
ids = substr(ids, 1,7)
ids = as.character(ids)

data = pflow_new[8:38,]
data = as.matrix(data)
data = as.data.frame(data)
rownames(data) = data[,1]
data = data[,2:ncol(data)]

#removed bcl-2 since there's no experimental results from Yanping
data = data[-2,]
features = rownames(data)
data = apply(data,2,as.numeric)
data = as.data.frame(data)
rownames(data) = features
data = t(data)
data = cbind(ids, labels, data)
rownames(data) = 1:nrow(data)
data = as.data.frame(data)
setwd("/Users/xuhaifeng/Documents/PhD_project_4/data/data_2023")
saveRDS(data, "baseline_pflow_all_cleaned.rds")

