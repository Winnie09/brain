dlist = list()
dlist[['2016_La:EmbryoMoleculeCounts']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds')
dlist[['2016_La:ESMoleculeCounts']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/ESMoleculeCounts.rds')
dlist[['2018_Fan_CellResearch']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')
dlist[['2018_Lake_NatBiotech:CerebellarHem']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/CerebellarHem.rds')
dlist[['2018_Lake_NatBiotech:FrontalCortex']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/FrontalCortex.rds')
dlist[['2018_Lake_NatBiotech:VisualCortex']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/VisualCortex.rds')
dlist[['2018_Li_Science:adult']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/adult.rds')
dlist[['2018_Li_Science:fetal']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds')
dlist[['allen']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds')

gbmgn <- row.names(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix/count.rds'))
gn <- table(unlist(sapply(dlist,row.names)))
gn <- names(gn)[gn==length(dlist)]
gn <- intersect(gn,gbmgn)

dlist <- sapply(dlist,function(i) i[gn,])
dlist <- do.call(cbind,dlist)

meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec_allenDEG/meta.rds')
meta <- meta[match(colnames(dlist),meta[,'cell']),]
time <- meta$time_super
time[time=='unknown'] <- 'adult'
meta$study <- sub('cell_.*','',meta$cell)
samp <- paste0(time,'_',meta$study)

d <- sapply(unique(samp),function(ss) {
  rowSums(dlist[,samp==ss])
})

d <- sweep(d,2,colSums(d)/1e6,'/')
d <- log2(d + 1)
library(preprocessCore)
dn <- dimnames(d)
d <- normalize.quantiles(d)
dimnames(d) <- dn

dcn <- sub('_.*','',colnames(d))
library(limma)

res <- topTable(eBayes(lmFit(d,design=cbind(1,dcn=='adult'))),n=nrow(d),coef=2)
gl <- row.names(res)[res[,'adj.P.Val'] < 0.05]  
saveRDS(gl,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/diffgene/res/time.rds')  


