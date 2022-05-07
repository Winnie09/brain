library(Seurat)
library(Matrix)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
clu <- Idents(d)
n <- names(clu)
clu <- as.numeric(as.character(clu)) + 1
names(clu) <- n
saveRDS(clu,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')

