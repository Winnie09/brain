library(Seurat)
library(ggplot2)
library(patchwork)
library(here)
setwd(here())
atlas <- readRDS('atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/humanAtlas_harmony.rds')
gbm <- readRDS('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')

atlas$RNA@var.features <- rownames(atlas$RNA@scale.data)
gbm$RNA@var.features <- rownames(gbm$RNA@scale.data)

anchors <- FindTransferAnchors(
  reference = atlas,
  query = gbm,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)
saveRDS(anchors, 'atlasGBM/mouseGBM/integrate/seuratv4/36nonNormal_seuratGene/res/anchors2.rds')

ct <- readRDS('atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/celltype_annotation.rds')

res <- MapQuery(
  anchorset = anchors,
  query = gbm,
  reference = atlas,
  refdata = ct,
  reference.reduction = "pca", 
  reduction.model = "umap"
)

saveRDS(res, 'atlasGBM/mouseGBM/integrate/seuratv4/36nonNormal_seuratGene/res/seurat_obj2.rds')

umap <- res@reductions$umap@cell.embeddings
saveRDS(umap, 'atlasGBM/mouseGBM/integrate/seuratv4/36nonNormal_seuratGene/res/umap_embeddings2.rds')











