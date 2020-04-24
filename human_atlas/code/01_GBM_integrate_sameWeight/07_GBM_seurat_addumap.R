library(Seurat)
d.integrated <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/GBM_integrated.rds')
ref.integrated <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/integrated.rds')
refmeta <- ref.integrated@meta.data

cn = colnames((d.integrated@assays)$RNA@counts)
scn = setdiff(cn,refmeta$cell)
p = sub('_.*','',scn)

meta = d.integrated@meta.data
tmp = cbind(meta,  refmeta[match(rownames(meta),rownames(refmeta)),])
v = tmp[grepl('GBM',tmp$orig.ident),'orig.ident']
tab = table(v)
v = paste0(v,'(',tab[match(v,names(tab))],')')
tmp[grepl('GBM',tmp$orig.ident),'study'] <- v
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

saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/GBM_integratedUMAP.rds')
