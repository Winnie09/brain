library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
# key = as.character(commandArgs(trailingOnly = T)[[1]])
key = 'atlasOnly_homolog_sub'
pca = readRDS(paste0('./brain/atlas/result/',key,'/pca/pca.rds'))
pcvar = apply(pca,2,var)
pdf(paste0('./brain/atlas/result/',key,'/pca/pcvar.pdf'),height=4,width=4)
plot(pcvar,pch=20,xlab='num.PC',ylab='variance')
dev.off()
source('./raisin/proc/cluster.R')
cluster(paste0('./brain/atlas/result/',key,'/'),hclu=F,useMNN=F,useUMAP=F,usePCA=T,num.PC=15,clunum=50,k=2)

