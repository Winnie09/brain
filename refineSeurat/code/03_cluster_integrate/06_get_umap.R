source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
mat <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/mat_c_transform.rds')
pr = PCA(mat, save.pca = TRUE, plot.statistics=F, result.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/',plot.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/plot/03_cluster_integrate/',findVariableGenes=F, numPC=10)
str(pr)
library(umap)
u = umap(pr)$layout
str(u)
rownames(u) = colnames(mat)
saveRDS(u, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/umap.rds')
