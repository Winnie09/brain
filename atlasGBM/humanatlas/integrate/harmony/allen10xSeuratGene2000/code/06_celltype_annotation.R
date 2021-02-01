gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/res/active.ident.rds')
tb <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/diff/celltype_annotation.csv', as.is = TRUE)
anno <- tb[,2]
names(anno) <- tb[,1]
gbm <- gbm
gbm.anno <- anno[gbm]
names(gbm.anno) <- names(gbm)
saveRDS(gbm.anno, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/res/celltype_annotation.rds')
