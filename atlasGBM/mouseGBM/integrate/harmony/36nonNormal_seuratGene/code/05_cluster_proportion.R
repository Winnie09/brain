library(Seurat)
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
library(RColorBrewer)
library(pheatmap)
rdir = 'atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/res/'  ##### 
pdir = 'atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/plot/'  ##### 
brain <- readRDS(paste0(rdir, 'atlasGBM_harmony.rds'))
clu <- as.character(brain@active.ident)
names(clu) <- colnames(brain@assays$RNA@counts)
saveRDS(clu, paste0(rdir, max(clu), 'clusters.rds'))

p <- as.character(brain@meta.data[,'Sample.ID'])
p <- sapply(p, function(i) ifelse(grepl('GBM',i), i, 'humanatlas'))
tab <- cumsum(table(clu))/length(clu)
prop <- sapply(unique(p),function(i) {
  tab <- table(clu[p==i])/sum(p==i)
  tmp <- rep(0,length(unique(clu)))
  names(tmp) <- sort(unique(clu))
  tmp[names(tab)] <- as.vector(tab)
  tmp
})
rownames(prop) <- paste0('cluster', rownames(prop))
saveRDS(prop, paste0(rdir, max(clu), 'cluster_proportion.rds'))

meta <- read.csv('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/doc/metas_20200329.csv',as.is=T)
meta <- meta[match(colnames(prop)[grepl('GBM',colnames(prop))], meta$Sample.ID), -1]
colann <- rbind(NA, meta)
rownames(colann) <- colnames(prop)

color.ls <- lapply(colnames(colann), function(i){
  v <- colorRampPalette(rev(brewer.pal(n = 7, name =
  "Spectral")))(length(unique(meta[, i])))
  names(v) <- unique(meta[, i])
  v
})
names(color.ls) <- colnames(colann)
png(paste0(pdir, 'heatmap_cluster_proportion.png'), width = 2000, height = 1800, res = 100)
pheatmap(prop, border_color = NA,
         annotation_col = colann,
         annotation_colors = color.ls,
         fontsize = 100)

dev.off()
