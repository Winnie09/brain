key = as.character(commandArgs(trailingOnly = T)[[1]])
# key = 'referenceVec_allenDEG'
library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas')
d.integrated <- readRDS(paste0('./result/',key,'/integrated.rds'))
ref.integrated <- readRDS('./result/integrated.rds')
refmeta <- ref.integrated@meta.data
# cn = colnames((d.integrated@assays)$RNA@counts)
# scn = setdiff(cn,refmeta$cell)
# p = sub('_.*','',scn)
meta = d.integrated@meta.data
tmp = cbind(meta,  refmeta[match(rownames(meta),refmeta$cell),]) 
v = tmp[grepl('GBM',tmp$orig.ident),'orig.ident']
tab = table(v)
v = paste0(v,'(',tab[match(v,names(tab))],')')
tmp[grepl('GBM',tmp$orig.ident),'study'] <- v
tmp[grepl('GBM',tmp$orig.ident),'cell'] <- rownames(tmp)[grepl('GBM',tmp$orig.ident)]
tmp = tmp[,-(4:6)]
saveRDS(tmp,paste0('./result/',key,'/meta.rds'))
d.integrated@meta.data = tmp
library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(d.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
d.integrated <- ScaleData(d.integrated, verbose = FALSE)
d.integrated <- RunPCA(d.integrated, npcs = 30, verbose = FALSE)
d.integrated <- RunUMAP(d.integrated, reduction = "pca", dims = 1:30)
saveRDS(d.integrated,paste0('./result/',key,'/GBM_integratedUMAP.rds'))

u = d.integrated@reductions$umap@cell.embeddings
dir.create(paste0('./result/',key,'/umap/'),showWarnings = F,recursive = T)
saveRDS(u,paste0('./result/',key,'/umap/umap.rds'))

pca = d.integrated@reductions$pca@cell.embeddings
dir.create(paste0('./result/',key,'/pca/'),showWarnings = F,recursive = T)
saveRDS(pca,paste0('./result/',key,'/pca/pca.rds'))
