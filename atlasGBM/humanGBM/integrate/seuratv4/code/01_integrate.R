library(Seurat)
library(ggplot2)
library(patchwork)
library(here)
setwd(here())
atlas <- readRDS('atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/res/humanAtlas_harmony.rds')

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
saveRDS(anchors, 'atlasGBM/humanGBM/integrate/seuratv4/res/anchors.rds')

ct <- readRDS('atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/res/celltype_annotation.rds')

identical(names(ct), colnames(atlas@assays$RNA@counts))
res <- MapQuery(
  anchorset = anchors,
  query = gbm,
  reference = atlas,
  refdata = ct,
  reference.reduction = "pca", 
  reduction.model = "umap"
)

saveRDS(res, 'atlasGBM/humanGBM/integrate/seuratv4/res/seurat_obj.rds')
