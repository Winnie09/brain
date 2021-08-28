library(ggplot2)
library(Seurat)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/project.rds')

## get GBM projected predicted celltype and the annotation
gbm.pred.clu = as.numeric(d@meta.data[, 'predicted.celltype'])
tb = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/humancortexatlas_celltype_annotation.csv', as.is = T)
gbm.pred.ct = tb[match(gbm.pred.clu, tb[,1]),3]
sort(table(gbm.pred.clu))
sort(table(gbm.pred.ct))
names(gbm.pred.ct) = rownames(d@meta.data)
d@meta.data[,'predicted.celltype_anno'] <- gbm.pred.ct[rownames(d@meta.data)]
saveRDS(gbm.pred.ct, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/gbm_predict_celltpye_anno.rds')

r <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/humanAtlas_harmony.rds')

## get atlas cell type annotation
atlas.ct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/celltype.rds')
r@meta.data[,'celltype.anno'] = atlas.ct[rownames(r@meta.data)]

## save seurat objects with cell type anno in meta.data
saveRDS(d, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/project_with_celltype_anno.rds')
saveRDS(r, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/humanAtlas_harmony_celltype_anno.rds')

