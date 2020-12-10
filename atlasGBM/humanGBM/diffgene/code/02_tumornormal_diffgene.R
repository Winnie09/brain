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
res <- row.names(res)[res[,'adj.P.Val'] < 0.05]  
saveRDS(res,file= paste0(rdir, 'tumornormal_diffgene_limma.rds'))

