library(copykat)
library(Matrix)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/count.rds')
samp <- sub('_.*','',colnames(d))
us <- unique(samp)
id <- as.numeric(commandArgs(trailingOnly = T))
d <- d[,samp==us[id]]
d <- d[rowSums(d) > 0,]
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/copykat/plot/')
copykat <- copykat(rawmat=as.matrix(d),n.cores=10,sam.name=us[id])
saveRDS(copykat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/copykat/res/',us[id],'.rds'))

