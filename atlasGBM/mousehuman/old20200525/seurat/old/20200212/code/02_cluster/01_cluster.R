library(Seurat)
library(ggplot2)
library(cowplot)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
pca = readRDS('./brain/atlas/res/pca/pca.rds')
pcvar = apply(pca,2,var)
pdf('./brain/atlas/res/pca/pcvar.pdf',height=4,width=4)
plot(pcvar,pch=20,xlab='num.PC',ylab='variance')
dev.off()
source('./raisin/proc/cluster.R')
cluster('./brain/atlas/res/',hclu=F,useMNN=F,useUMAP=F,usePCA=T,num.PC=10)

