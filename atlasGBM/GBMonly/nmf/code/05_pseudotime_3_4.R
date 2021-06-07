factorscore <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pathwayscore.rds')

library(TSCAN)
pd <- cbind(factorscore[,3],factorscore[,4])
clu <- kmeans(pd,2)$cluster
ord <- TSCANorder(exprmclust(t(pd),cluster=clu,reduce = F),orderonly = T)
saveRDS(ord,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pseudotime_3_4.rds')


