source('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/resource/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/')
cnt = readRDS('./processed/log2/mat.rds')
pr <- pca_scranFilter(cnt)
saveRDS(pr[c('sdev','x')],'./result/pca.rds')


