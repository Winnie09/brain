setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result')
key = as.character(commandArgs(trailingOnly = T)[[1]])
# key = 'atlasOnly_homolog_sub'
# key = 'atlas_GBM_homolog_sub'
aggregateFunc <- function(m){
  clus = clu[match(colnames(m), names(clu))]
  me <- sapply(unique(clus),function(i) rowMeans(m[,clus==i,drop=F]))
  colnames(me) <- unique(clus)
  me = me[,order(as.numeric(colnames(me)),decreasing=F)]
  colnames(me) = paste0('clu',colnames(me))
  me
}
library(Seurat)
d = readRDS(paste0('./',key,'/integrated.rds'))
count = as.matrix(d@assays$RNA@counts)
rm(d)
clu = readRDS(paste0('./',key,'/cluster/cluster.rds'))

a = gsub('La:','La_', colnames(count))
ds = gsub(':.*','',a)

meta = readRDS(paste0('./',key,'/meta.rds'))
count = count[, meta$cell]
dir.create(paste0('./',key,'/cluster_mean/species/'), recursive = T)
for (i in unique(meta$species)){
  print(i)
  mat = count[, meta$species == i]
  agmat = aggregateFunc(mat)
  saveRDS(agmat, paste0('./',key,'/cluster_mean/species/', i, '.rds'))
}

dir.create(paste0('./',key,'/cluster_mean/study/'), recursive = T)
for (i in unique(ds)){
  print(i)
  mat = count[, ds == i]
  agmat = aggregateFunc(mat)
  saveRDS(agmat, paste0('./',key,'/cluster_mean/study/', i, '.rds'))
}
