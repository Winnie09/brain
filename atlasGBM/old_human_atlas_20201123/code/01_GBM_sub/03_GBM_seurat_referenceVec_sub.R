# ## select genes
# gl <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/res/gl/GBM_allen.rds')
# gl <- row.names(gl)[gl[,'fdr'] > 0.1]
# ctgl <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/res/gl/allen_ct.rds')
# ctgl <- unique(unlist(sapply(ctgl,function(i) {
#   row.names(i)[i[,'fdr'] < 0.05 & abs(i[,'log2FoldChange']) > 1]
# })))
# gl <- intersect(gl,ctgl)

library(Seurat)
countfunc <- function(d) {
  log2(t(t(d)/(colSums(d) / 10000) + 1))
}

## load atlas
dlist = list()
dlist[['2016_La:EmbryoMoleculeCounts']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds'))
dlist[['2016_La:ESMoleculeCounts']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/ESMoleculeCounts.rds'))
dlist[['2016_La:iPSMoleculeCounts']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/iPSMoleculeCounts.rds'))

dlist[['2018_Fan_CellResearch']] <- log2(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')+1)

dlist[['2018_Lake_NatBiotech:CerebellarHem']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/CerebellarHem.rds'))
dlist[['2018_Lake_NatBiotech:FrontalCortex']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/FrontalCortex.rds'))
dlist[['2018_Lake_NatBiotech:VisualCortex']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/VisualCortex.rds'))

dlist[['2018_Li_Science:adult']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/adult.rds')
dlist[['2018_Li_Science:fetal']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds'))

dlist[['allen']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds'))

gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix/count.rds')

gn <- table(unlist(sapply(dlist,row.names)))
gn <- names(gn)[gn==length(dlist)]
gn <- intersect(gn,row.names(gbm))
gbm <- gbm[gn,]
p = sub('_.*','',colnames(gbm))

for (i in names(dlist)) {
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]][gn,],project=i))
}

for (sp in unique(p)) {
  id <- which(p==sp)
  id <- sample(id,min(1000,length(id)))
  dlist[[sp]] <- FindVariableFeatures(CreateSeuratObject(countfunc(gbm[,id]),project=sp))
}

reference.vec = c('2018_Li_Science:fetal','2018_Li_Science:adult','2018_Lake_NatBiotech:VisualCortex','2018_Lake_NatBiotech:FrontalCortex','2018_Lake_NatBiotech:CerebellarHem','2018_Fan_CellResearch','2016_La:iPSMoleculeCounts','2016_La:ESMoleculeCounts','2016_La:EmbryoMoleculeCounts')
d.anchors <- FindIntegrationAnchors(object.list = dlist, reference = grep('GBM',names(dlist),invert = T))
d.integrated <- IntegrateData(anchorset = d.anchors)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/',showWarnings = F, recursive = T)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/GBM_integrated_sub.rds')
