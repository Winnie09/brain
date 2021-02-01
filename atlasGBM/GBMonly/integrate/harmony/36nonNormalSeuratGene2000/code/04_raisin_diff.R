library(Seurat)
library(Matrix)
library(here)
setwd(here())
k <- readRDS('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
source('/home-4/zji4@jhu.edu/scratch/raisin/software/raisin/raisin.R')
d <- k$RNA@data
meta <- k@meta.data
clu <- Idents(k)
for (rid in as.character(sort(unique(clu)))) {
  t1cell <- names(clu)[which(clu==rid)]
  t2cell <- names(clu)[which(clu!=rid)]
  ctype <- paste0(meta[colnames(d),'Sample.ID'],':',as.character(clu[colnames(d)]==rid))
  ut <- unique(ctype)
  design <- data.frame(sample=ut,feature=sub('.*:','',ut),individual=sub(':.*','',ut),stringsAsFactors = F)
  tab <- table(design[,3])
  sel <- names(tab)[tab==2]
  design <- design[design[,3] %in% sel,]
  cid <- which(ctype %in% design[,1])
  if (sum(tab==2) > 1) {
    fit <- RAISINfit(d[,cid],ctype[cid],design,testtype='paired',ncores=10)
    res <- RAISINtest(fit)[,-3]
    res <- cbind(res,within=rowMeans(fit$mean[row.names(res),design[design[,2]=='TRUE',1]]),outside=rowMeans(fit$mean[row.names(res),design[design[,2]=='FALSE',1]]))
    dir.create('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/diff/')
    saveRDS(res,file=paste0('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/diff/',rid,'.rds'))
    write.csv(res,file=paste0('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/diff/',rid,'.csv'), row.names = TRUE)
  }
}

