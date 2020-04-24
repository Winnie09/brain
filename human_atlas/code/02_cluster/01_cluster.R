key = as.character(commandArgs(trailingOnly = T)[[1]])
# key = 'referenceVec'
# key = 'referenceVec_allenDEG'
library(Seurat)
library(ggplot2)
library(cowplot)
setwd(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/',key))
pca = readRDS(paste0('./pca/pca.rds'))
pcvar = apply(pca,2,var)
pdf('./pca/pcvar.pdf',height=4,width=4)
plot(pcvar,pch=20,xlab='num.PC',ylab='variance')
dev.off()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/cluster.R')
cluster(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/',key),hclu=F,useMNN=F,useUMAP=F,usePCA=T,num.PC=10)

