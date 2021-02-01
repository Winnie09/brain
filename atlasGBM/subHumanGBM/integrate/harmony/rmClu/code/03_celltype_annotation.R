library(Seurat)
seu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/subHumanGBM/integrate/harmony/rmClu/result/sub_atlasGBM_harmony.rds')
meta <- seu@meta.data

gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/celltype_annotation.rds')
human <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/res/celltype_annotation.rds')

anno <- c(human, gbm)

meta$celltype_anno <- anno[meta$cell]
saveRDS(anno[meta$cell], '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/subHumanGBM/integrate/harmony/rmClu/res/celltype_annotation.rds')

seu@meta.data <- meta

p <- DimPlot(seu, group.by = 'celltype_anno', split.by = 'celltype_anno', ncol = 4)
library(ggplot2)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/subHumanGBM/integrate/harmony/rmClu/plot/umap_celltype_anno_facet.png', p, width = 15, height = 15, dpi = 200)

p <- DimPlot(seu, group.by = 'celltype_anno')
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/subHumanGBM/integrate/harmony/rmClu/plot/umap_celltype_anno.png', p, width = 9, height = 3.5, dpi = 200)

