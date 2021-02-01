library(here)
setwd(here())

tb = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/Marker_Genes.csv',header=T,sep=',',as.is=T)
gvec <- unlist(tb)

gvec <- gvec[!is.na(gvec)]
str(gvec)
gvec <- unique(gvec)


## GBMonly
pdir <- 'atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/plot/'
brain <- readRDS('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
mat <- brain@assays$RNA@counts
gvec.tmp <- intersect(rownames(mat), gvec)
clu <- as.character(brain@active.ident)
names(clu) <- colnames(mat)
agg <- sapply(0:(length(unique(clu))-1), function(i){
  print(i)
  id <- which(clu==i)
  print(str(id))
  tmp <- rowMeans(as.matrix(mat[gvec.tmp, id]))
})
colnames(agg) <- paste0('cluster', 0:(length(unique(clu))-1))
library(pheatmap)
pdf(paste0(pdir, 'markerGene_cluster_mean_hm.pdf'), width = 7, height = 14)
pheatmap(agg,
         border_color = NA)
dev.off()
  

## human atlas
pdir <- 'atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/plot/'
brain <- readRDS('atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/res/humanAtlas_harmony.rds')
mat <- brain@assays$RNA@counts
gvec.tmp <- intersect(rownames(mat), gvec)
clu <- as.character(brain@active.ident)
names(clu) <- colnames(mat)
agg <- sapply(0:(length(unique(clu))-1), function(i){
  print(i)
  id <- which(clu==i)
  print(str(id))
  tmp <- rowMeans(as.matrix(mat[gvec.tmp, id]))
})
colnames(agg) <- paste0('cluster', 0:(length(unique(clu))-1))
library(pheatmap)
pdf(paste0(pdir, 'markerGene_cluster_mean_hm.pdf'), width = 7, height = 14)
pheatmap(agg,
         border_color = NA)
dev.off()



## mouse atlas
pdir <- 'atlasGBM/mouseatlas/integrate/harmony/seuratGene/plot/'
brain <- readRDS('atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/humanAtlas_harmony.rds')
mat <- brain@assays$RNA@counts
gvec.tmp <- intersect(rownames(mat), gvec)
clu <- as.character(brain@active.ident)
names(clu) <- colnames(mat)

agg <- sapply(0:(length(unique(clu))-1), function(i){
  print(i)
  id <- which(clu==i)
  print(str(id))
  tmp <- rowMeans(as.matrix(mat[gvec.tmp, id]))
})
colnames(agg) <- paste0('cluster', 0:(length(unique(clu))-1))
library(pheatmap)
pdf(paste0(pdir, 'markerGene_cluster_mean_hm.pdf'), width = 7, height = 14)
pheatmap(agg,
         border_color = NA)
dev.off()

