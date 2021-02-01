setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
# key = as.character(commandArgs(trailingOnly = T)[[1]])
key = 'atlas_GBM_homolog_sub'
pca = readRDS(paste0('./brain/atlas/result/',key,'/pca/pca.rds'))
pcvar = apply(pca,2,var)
pdf(paste0('./brain/atlas/result/',key,'/pca/pcvar.pdf'),height=4,width=4)
plot(pcvar,pch=20,xlab='num.PC',ylab='variance')
dev.off()
source('./raisin/proc/cluster.R')
savedir = paste0('./brain/atlas/result/',key,'/')
d <- readRDS(paste0(savedir,'/pca/pca.rds')) 
d <- d[,1:15]
set.seed(12345)
res = kmeans(d, 30)$cluster
dir.create(paste0(savedir,'cluster/'))
saveRDS(res, paste0(savedir,'/cluster/cluster.rds'))

