data = as.character(commandArgs(trailingOnly = T)[[1]])
# data = 'atlasOnly_homolog_sub'
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
library(Seurat)
d.integrated <- readRDS(paste0('./brain/atlas/result/',data,'/integrated.rds'))
library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(d.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
d.integrated <- ScaleData(d.integrated, verbose = FALSE)
d.integrated <- RunPCA(d.integrated, npcs = 30, verbose = FALSE)
pca = d.integrated@reductions$pca@cell.embeddings
dir.create(paste0('./brain/atlas/result/',data,'/pca/'),showWarnings = F,recursive = T)
saveRDS(pca,paste0('./brain/atlas/result/',data,'/pca/pca.rds'))

d.integrated <- RunUMAP(d.integrated, reduction = "pca", dims = 1:15)
u = d.integrated@reductions$umap@cell.embeddings
dir.create(paste0('./brain/atlas/result/',data,'/umap/'),showWarnings = F,recursive = T)
saveRDS(u,paste0('./brain/atlas/result/',data,'/umap/umap.rds'))

