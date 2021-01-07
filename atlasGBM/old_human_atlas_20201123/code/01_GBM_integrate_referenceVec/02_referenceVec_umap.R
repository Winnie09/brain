library(Seurat)
d.integrated <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/GBM_integrated.rds')
ref.integrated <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/integrated.rds')
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
saveRDS(tmp,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/meta.rds')
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
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/GBM_integratedUMAP.rds')

u = d.integrated@reductions$umap@cell.embeddings
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/umap/',showWarnings = F,recursive = T)
saveRDS(u,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/umap/umap.rds')

pca = d.integrated@reductions$pca@cell.embeddings
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/pca/',showWarnings = F,recursive = T)
saveRDS(pca,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/pca/pca.rds')

