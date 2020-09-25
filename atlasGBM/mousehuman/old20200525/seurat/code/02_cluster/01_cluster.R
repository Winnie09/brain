library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/result/')
# key = as.character(commandArgs(trailingOnly = T)[[1]])
# key = 'atlasOnly_homolog_sub'
key = 'atlas_GBM_homolog_sub'
pca = readRDS(paste0('./',key,'/pca/pca.rds'))
pcvar = apply(pca,2,var)
pdf(paste0('./',key,'/pca/pcvar.pdf'),height=4,width=4)
plot(pcvar,pch=20,xlab='num.PC',ylab='variance')
dev.off()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/cluster.R')
cluster(paste0('./',key,'/'),hclu=F,useMNN=F,useUMAP=F,usePCA=T,num.PC=10,clunum=50,k=2)

