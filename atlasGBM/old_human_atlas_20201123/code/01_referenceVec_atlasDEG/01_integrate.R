
########################
## select genes
c <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/diffgene/res/celltype.rds')
t <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/diffgene/res/time.rds')
n <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/diffgene/res/normaltumor.rds')
gl <- setdiff(union(c,t),n)

### 
library(Seurat)
countfunc <- function(d) {
  log2(t(t(d)/(colSums(d) / 10000) + 1))
}

## load atlas
dlist = list()
dlist[['2016_La:EmbryoMoleculeCounts']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds')[])
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
print('intersect genes of atlas and GBM\n')
str(gn)
gn <- intersect(gn,gl) ## use DEG
print('intersect genes with DEGs\n')
str(gn)

gbm <- gbm[gn,]
p = sub('_.*','',colnames(gbm))

for (i in names(dlist)) {
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]][gn,],project=i))
}

for (sp in unique(p)) {
  dlist[[sp]] <- FindVariableFeatures(CreateSeuratObject(countfunc(gbm[,p==sp]),project=sp))
}

reference.vec = c('2018_Li_Science:fetal','2018_Li_Science:adult','2018_Lake_NatBiotech:VisualCortex','2018_Lake_NatBiotech:FrontalCortex','2018_Lake_NatBiotech:CerebellarHem','2018_Fan_CellResearch','2016_La:iPSMoleculeCounts','2016_La:ESMoleculeCounts','2016_La:EmbryoMoleculeCounts')
d.anchors <- FindIntegrationAnchors(object.list = dlist, reference = grep('GBM',names(dlist),invert = T))
d.integrated <- IntegrateData(anchorset = d.anchors)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec_atlasDEG/',showWarnings = F, recursive = T)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec_atlasDEG/integrated.rds')
