library(Seurat)
library(Matrix)
library(here)
setwd(here())
k <- readRDS('atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/res/atlasGBM_harmony.rds')
source('/home-4/zji4@jhu.edu/scratch/raisin/software/raisin/raisin.R')
d <- k$RNA@data
meta <- k@meta.data
clu <- Idents(k)
for (rid in sort(unique(clu))) {
  t1cell <- names(clu)[which(clu==rid)]
  t2cell <- names(clu)[which(clu!=rid)]
  ctype <- paste0(meta[colnames(d),'study'],':',as.character(clu[colnames(d)]==rid))
  ut <- unique(ctype)
  design <- data.frame(sample=ut,feature=sub('.*:','',ut),individual=sub(':.*','',ut),stringsAsFactors = F)
  tab <- table(design[,3])
  sel <- names(tab)[tab==2]
  design <- design[design[,3] %in% sel,]
  cid <- which(ctype %in% design[,1])
  if (sum(tab > 1) == 2) {
    fit <- RAISINfit(d[,cid],ctype[cid],design,testtype='paired',ncores=10)
    res <- RAISINtest(fit)[,-3]
    res <- cbind(res,within=rowMeans(fit$mean[row.names(res),design[design[,2]=='TRUE',1]]),outside=rowMeans(fit$mean[row.names(res),design[design[,2]=='FALSE',1]]))
    dir.create('atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/diff/')
    saveRDS(res,file=paste0('atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/diff/',rid,'.rds'))
  }
}
