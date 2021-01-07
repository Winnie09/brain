result_path = as.character(commandArgs(trailingOnly = T)[1])
key = as.character(commandArgs(trailingOnly = T)[[2]])
if(!exists('result_path')) print('please provide a working directory one level above the result path containing umap, meta, cluster')
setwd(result_path)
# result_path ='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/'
# key = 'sub'
# key = 'full'
# result_path = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/'
# key = 'referenceVec'
# key = 'referenceVec_allenDEG'

library(Seurat)
library(ggplot2)
library(cowplot)
u = readRDS(paste0('./result/',key,'/umap/umap.rds'))
meta = readRDS(paste0('./result/',key,'/meta.rds'))
clu = readRDS(paste0('./result/',key,'/cluster/cluster.rds'))  

ct = as.character(meta$celltype_super)  
id = intersect(which(!is.na(meta$celltype)), which(ct=='neuronal'))
ct[id] <- as.character(meta$celltype)[id]  ##### detailed neuronal

ct[grepl('GBM',meta$study)] <- meta[grepl('GBM',meta$study),'study']
names(ct) = rownames(u)

sample = rep('atlas',nrow(meta))
sample[grep('GBM',meta$study)] = as.character(meta[grep('GBM',meta$study),'study'])
clu = clu[rownames(u)]
allenct <- sapply(1:max(clu),function(i) {
  names(which.max(table(ct[clu==i & !grepl('GBM',ct)])))
})
gbmenrich <- sapply(1:max(clu),function(i) {
  round(mean(grepl('GBM',ct[clu==i]))/mean(grepl('GBM',ct)),4)
})
percent <- sapply(unique(sample),function(i) {
  sapply(1:max(clu),function(j) {
    round(sum(clu==j & sample==i)/sum(sample==i),4)
  })
})
res <- data.frame(cluster=1:max(clu),atlas_celltype=allenct,gbm_enrichment=gbmenrich,percent,stringsAsFactors = F)
colnames(res)[4:ncol(res)] <- colnames(percent)
dir.create(paste0('./originCt/result/',key),showWarnings = F,recursive = T)
write.csv(res,paste0('./originCt/result/',key,'/originCt_prop.csv'),row.names = F)
saveRDS(res,paste0('./originCt/result/',key,'/originCt_prop.rds'))

  
