library(Seurat)
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(scattermore)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
rdir <- 'atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/res/'  ####
pdir <- 'atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/plot/'  ####
brain <- readRDS(paste0(rdir, 'humanAtlas_harmony.rds'))
p <- DimPlot(brain, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
ggsave(paste0(pdir, max(as.numeric(brain@active.ident)),'clusters.png'),
       width = 7,
       height = 5,
       dpi = 200)


for (i in colnames(brain@meta.data)){
  if (length(unique(brain@meta.data[,i]))< 100 & length(unique(brain@meta.data[,i])) > 1){
    print(i)
    p <- DimPlot(brain, reduction = "umap", group.by = i,label = TRUE, repel = TRUE) + NoLegend()
    ggsave(paste0(pdir, 'seurat_plot_', i, '.png'),
           width = 7,
           height = 5,
           dpi = 200)
  }
}
  








