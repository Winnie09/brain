data = as.character(commandArgs(trailingOnly = T)[[1]])
num.PC = as.numeric(commandArgs(trailingOnly = T)[[2]])
# data = 'atlasOnly_homolog_sub'
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
pca = readRDS(paste0('./brain/atlas/result/',data,'/pca/pca.rds'))
pca = pca[,1:num.PC]
library(umap)
set.seed(12345)
u = umap(pca)$layout
dir.create(paste0('./brain/atlas/result/',data,'/umap',num.PC,'/'),showWarnings = F,recursive = T)
saveRDS(u,paste0('./brain/atlas/result/',data,'/umap',num.PC,'/umap.rds'))
