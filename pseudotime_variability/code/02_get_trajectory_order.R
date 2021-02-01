setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain')
# setwd('/Users/wenpinhou/Dropbox/brain/')
pos = readRDS('./pseudotime/old/20200223/plot/plot/pseudotime_order.rds')
u <- readRDS('./atlas/old/20200212/result/sub/umap/umap.rds')
library(TSCAN)
set.seed(12345)
clu <- kmeans(u,12,iter.max = 1000)$cluster
em <- exprmclust(t(u),reduce=F,cluster = clu)
order <- TSCANorder(em,listbranch = T,orderonly = F,MSTorder = c(6,10,3))
saveRDS(order, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/order.rds')

