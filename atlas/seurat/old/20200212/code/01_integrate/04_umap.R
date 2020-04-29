setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
library(Seurat)
d.integrated <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/res/sub/integrated_GBM_reference.rds')
meta = d.integrated@meta.data  ## 282449
fullmeta = readRDS('./brain/atlas/res/full/meta.rds')
length(intersect(rownames(meta),fullmeta$cell))
meta2 = cbind(meta, fullmeta[match(rownames(meta), fullmeta$cell), ])
saveRDS(meta2,'./brain/atlas/res/sub/meta.rds')

d.integrated@meta.data = meta2
library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(d.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
d.integrated <- ScaleData(d.integrated, verbose = FALSE)
d.integrated <- RunPCA(d.integrated, npcs = 30, verbose = FALSE)
d.integrated <- RunUMAP(d.integrated, reduction = "pca", dims = 1:30)
u = d.integrated@reductions$umap@cell.embeddings
dir.create('./brain/atlas/res/umap/',showWarnings = F,recursive = T)
saveRDS(u,'./brain/atlas/res/umap/umap.rds')

pca = d.integrated@reductions$pca@cell.embeddings
dir.create('./brain/atlas/res/pca/',showWarnings = F,recursive = T)
saveRDS(pca,'./brain/atlas/res/pca/pca.rds')
