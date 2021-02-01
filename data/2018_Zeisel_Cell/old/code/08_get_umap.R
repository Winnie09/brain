setwd("/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Zeisel_Cell/")
pr = readRDS('./result/pca.rds')
pr = pr$x[,1:20]
library(umap)
set.seed(12345)
umap = umap(pr)$layout
saveRDS(umap,'./result/umap.rds')

