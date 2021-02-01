gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/active.ident.rds')
tb <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/integrate/harmony/seuratGene/diff/celltype_annotation.csv', as.is = TRUE)
anno <- tb[,2]
names(anno) <- tb[,1]
gbm <- gbm
gbm.anno <- anno[gbm]
names(gbm.anno) <- names(gbm)
saveRDS(gbm.anno, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/celltype_annotation.rds')


