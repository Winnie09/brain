library(Seurat)
countfunc <- function(d) {
  log2(t(t(d)/(colSums(d) / 10000) + 1))
}

### load atlas
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


dlist[['GBM']] <- FindVariableFeatures(CreateSeuratObject(countfunc(gbm),project='GBM'))


### load reference
ref.integrated <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/integrated.rds')

### query GBM
d.query <- dlist[grep('GBM',names(dlist))]
d.anchors <- FindTransferAnchors(reference = ref.integrated, query = d.query, dims = 1:30)
predictions <- TransferData(anchorset = d.anchors, refdata = ref.integrated$celltype, dims = 1:30)
d.query <- AddMetaData(d.query, metadata = predictions)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/queryGBM/',showWarnings = F, recursive = T)
saveRDS(d.query,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/queryGBM/GBM_integrated.rds')
