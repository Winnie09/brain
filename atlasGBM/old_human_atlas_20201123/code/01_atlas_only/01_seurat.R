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

gn <- table(unlist(sapply(dlist,row.names)))
gn <- names(gn)[gn==length(dlist)]

for (i in names(dlist)) {
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]][gn,],project=i))
}

d.anchors <- FindIntegrationAnchors(object.list = dlist)
d.integrated <- IntegrateData(anchorset = d.anchors)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/integrated.rds')

