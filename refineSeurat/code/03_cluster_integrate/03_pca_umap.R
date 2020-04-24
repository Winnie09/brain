library(Seurat)
library(data.table)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
source('./resource/function.R')
d.integrated = readRDS('./brain/refineSeurat/result/03_cluster_integrate/integrated.rds')
mat = as.matrix(d.integrated@assays$RNA@counts)

meta = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/meta.rds')

d.integrated@meta.data = meta
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
dir.create('./brain/refineSeurat/result/03_cluster_integrate/umap/', showWarnings = F, recursive = T)
saveRDS(u,'./brain/refineSeurat/result/03_cluster_integrate/umap/umap.rds')

pca = d.integrated@reductions$pca@cell.embeddings
dir.create('./brain/refineSeurat/result/03_cluster_integrate/pca/',showWarnings = F,recursive = T)
saveRDS(pca,'./brain/refineSeurat/result/03_cluster_integrate/pca/pca.rds')
