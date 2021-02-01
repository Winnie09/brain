library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')

mat1 <-  log2(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')+1)
row.names(mat1) <- toupper(row.names(mat1)) ## mouse
mat1 <- mat1[!duplicated(row.names(mat1)),]


mat2 <- log2TPM_from_fludigm_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds'))
row.names(mat2) <- toupper(row.names(mat2))  ## mouse
mat2 <- mat2[!duplicated(row.names(mat2)),]


gn <- intersect(rownames(mat1), rownames(mat2))
mat = cbind(mat1[gn,],mat2[gn,])

pr = PCA(mat, save.pca = TRUE, plot.statistics=TRUE, result.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/02_pca/',plot.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/02_pca/')

library(umap)
u = umap(pr)$layout
rownames(u) = rownames(mat)
saveRDS(u, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/02_pca/umap.rds')




