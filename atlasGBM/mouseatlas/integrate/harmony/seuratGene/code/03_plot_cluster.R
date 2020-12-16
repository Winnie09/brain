library(Seurat)
library(here)
setwd(here())
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(scattermore)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
rdir = 'atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/'  ##### 
pdir = 'atlasGBM/mouseatlas/integrate/harmony/seuratGene/plot/'  ##### 
brain <- readRDS(paste0(rdir, 'humanAtlas_harmony.rds'))

p <- DimPlot(brain, label = TRUE, repel = TRUE) + NoLegend()
ggsave(paste0(pdir, max(as.numeric(brain@active.ident)),'clusters.png'),p,
       width = 7,
       height = 5,
       dpi = 200)

for (i in colnames(brain@meta.data)[12: ncol(brain@meta.data)]){
  if (length(unique(brain@meta.data[,i])) < 100 & length(unique(brain@meta.data[,i])) > 0){
    print(i)  
    p <- DimPlot(brain, reduction = "umap", group.by = i,label = TRUE, repel = TRUE) + NoLegend()
    ggsave(paste0(pdir, 'seurat_plot_', i, '.png'),
           width = 7,
           height = 5,
           dpi = 200)}
}
  
