library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/humanGBM/diffgene/res/'

mat =  readRDS('./atlasGBM/humanGBM/data/combine_mat.rds')
meta =  readRDS('./atlasGBM/humanGBM/data/meta_allcell.rds')
donor <-  paste0(meta$study, ';',meta$gender,';', meta$age,';', meta$location, ';',meta$time, ';',meta$donor)

dmean <- sapply(unique(donor),function(i) colSums(t(mat[,donor==i])))
dmean <- dmean[, !grepl('2014_Pollen_NatBiotech', colnames(dmean))]

study <- sub(';.*', '', unique(donor))
study <- study[!grepl('2014_Pollen_NatBiotech', study)]


d=dmean
d <- sweep(d,2,colSums(d)/1e6,'/')
d <- log2(d + 1)
library(preprocessCore)
dn <- dimnames(d)
d <- normalize.quantiles(d)
dimnames(d) <- dn

library(limma)
res <- topTable(eBayes(lmFit(d,design=cbind(1,grepl('^GBM',study)))),n=nrow(d),coef=2)
res <- res[order(res[,'adj.P.Val'], -abs(res[,'logFC'])),]
saveRDS(res,file= paste0(rdir, 'tumornormal_diffgene_limma_res.rds'))

diffgene <- row.names(res)[res[,'adj.P.Val'] < 0.05]  
saveRDS(diffgene,file= paste0(rdir, 'tumornormal_diffgene_limma.rds'))

 tb = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/Marker_Genes.csv',header=T,sep=',',as.is=T)
 v = as.vector(unlist(tb))
pd <- data.frame(diffgene = diffgene, fdr = res[diffgene, 'adj.P.Val'], marker = sapply(diffgene, function(i) ifelse(i %in% v, 'Marker', 'unknown')))
library(ggplot2)
library(RColorBrewer)
p1 <- ggplot(data = pd, aes(x = marker, y = fdr, color = marker)) + 
  geom_violin(scale = 'width', alpha = 0.5) +
  geom_jitter(alpha = 0.5, size = 0.5) +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2')

p2 <- ggplot(data = pd[pd[,2]<0.01,], aes(x = marker, y = fdr, color = marker)) + 
  geom_violin(scale = 'width', alpha = 0.5) +
  geom_jitter(alpha = 0.5, size = 0.5) +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2')


p3 <- ggplot(data = pd[pd[,2]<0.001,], aes(x = marker, y = fdr, color = marker)) + 
  geom_violin(scale = 'width', alpha = 0.5) +
  geom_jitter(alpha = 0.5, size = 0.5) +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2')

p4 <- ggplot(data = pd[pd[,2]<1e-5,], aes(x = marker, y = fdr, color = marker)) + 
  geom_violin(scale = 'width', alpha = 0.5) +
  geom_jitter(alpha = 0.5, size = 0.5) +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2')
pdf(paste0(rdir, 'tumornormal_diffgene_limma_compared_with_markerGenes.rds'), width = 8, height = 6)
gridExtra::grid.arrange(p1,p2,p3,p4, nrow = 2)
dev.off()
