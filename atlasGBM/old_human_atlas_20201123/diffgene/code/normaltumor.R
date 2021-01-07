## load atlas
dlist = list()
dlist[['2016_La:EmbryoMoleculeCounts']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds'))
dlist[['2016_La:ESMoleculeCounts']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/ESMoleculeCounts.rds'))
dlist[['2018_Fan_CellResearch']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds'))
dlist[['2018_Lake_NatBiotech:CerebellarHem']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/CerebellarHem.rds'))
dlist[['2018_Lake_NatBiotech:FrontalCortex']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/FrontalCortex.rds'))
dlist[['2018_Lake_NatBiotech:VisualCortex']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/VisualCortex.rds'))
dlist[['2018_Li_Science:adult']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/adult.rds'))
dlist[['2018_Li_Science:fetal']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds'))
dlist[['allen']] <- rowSums(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds'))

d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix/count.rds')
p <- sub('_.*','',colnames(d))
dmean <- sapply(unique(p),function(i) rowSums(d[,p==i]))
rm('d')

gn <- table(unlist(sapply(dlist,names)))
gn <- names(gn)[gn==length(dlist)]
gn <- intersect(gn,row.names(dmean))
d <- sapply(dlist,function(i) i[gn])
d <- cbind(d,dmean[gn,])
d <- sweep(d,2,colSums(d)/1e6,'/')
d <- log2(d + 1)
library(preprocessCore)
dn <- dimnames(d)
d <- normalize.quantiles(d)
dimnames(d) <- dn

library(limma)
res <- topTable(eBayes(lmFit(d,design=cbind(1,grepl('^GBM',colnames(d))))),n=nrow(d),coef=2)
res <- row.names(res)[res[,'adj.P.Val'] < 0.05]  
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/diffgene/res/normaltumor.rds')
