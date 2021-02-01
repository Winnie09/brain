library(Seurat)
countfunc <- function(d) {
  log2(t(t(d)/(colSums(d) / 10000) + 1))
}
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
## load atlas
dlist = list()
dlist[['2018_Fan_CellResearch']] <- log2(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')+1)
dlist[['2018_Lake_NatBiotech:CerebellarHem']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/CerebellarHem.rds'))
dlist[['2018_Lake_NatBiotech:FrontalCortex']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/FrontalCortex.rds'))
dlist[['2018_Lake_NatBiotech:VisualCortex']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/VisualCortex.rds'))
dlist[['2018_Li_Science:adult']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/adult.rds')
dlist[['2018_Li_Science:fetal']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds'))
dlist[['allen']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds'))
dlist[['2018_Rosenberg']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/proc/expr.rds'))
row.names(dlist[['2018_Rosenberg']]) <- toupper(row.names(dlist[['2018_Rosenberg']]))
dlist[['2018_Rosenberg']] <- dlist[['2018_Rosenberg']][!duplicated(row.names(dlist[['2018_Rosenberg']])),]
dlist[['2018_Zeisel']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/mat.rds'))
row.names(dlist[['2018_Zeisel']]) <- toupper(row.names(dlist[['2018_Zeisel']]))
dlist[['2018_Zeisel']] <- dlist[['2018_Zeisel']][!duplicated(row.names(dlist[['2018_Zeisel']])),]

dlist[['2017_Nowakowski_Science']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2017_Nowakowski_Science/proc/log2tpm.rds')  ## add
dlist[['2017_Nowakowski_Science']] <- dlist[['2017_Nowakowski_Science']][!duplicated(row.names(dlist[['2017_Nowakowski_Science']])),]
dlist[['2016_La:EmbryoMoleculeCounts']] <- log2TPM_from_fludigm_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds'))
dlist[['2016_La:ESMoleculeCounts']] <- log2TPM_from_fludigm_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/ESMoleculeCounts.rds'))
dlist[['2016_La:MouseEmbryoMoleculeCounts.rds']] <- log2TPM_from_fludigm_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/MouseEmbryoMoleculeCounts.rds'))
dlist[['2016_La:MouseAdultDAMoleculeCounts.rds']] <- log2TPM_from_fludigm_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/MouseAdultDAMoleculeCounts.rds'))
dlist[['2016_Tasic_NatNeuro']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_Tasic_NatNeuro/proc/log2tpm.rds')
dlist[['2016_Tasic_NatNeuro']] <- dlist[['2016_Tasic_NatNeuro']][!duplicated(row.names(dlist[['2016_Tasic_NatNeuro']])),]
dlist[['2015_Zeisel_Science']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Zeisel_Science/proc/log2tpm.rds') ## add
dlist[['2015_Zeisel_Science']] <- dlist[['2015_Zeisel_Science']][!duplicated(row.names(dlist[['2015_Zeisel_Science']])),]

dlist[['allen']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds'))

for (i in names(dlist)) {
  set.seed(12345)
  if (ncol(dlist[[i]]) > 20000) dlist[[i]] <- dlist[[i]][,sample(1:ncol(dlist[[i]]),20000)]
}

# gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix/count.rds')
gn <- table(unlist(sapply(dlist,row.names)))
gn <- names(gn)[gn==length(dlist)]
# gn <- intersect(gn,row.names(gbm))
# gbm <- gbm[gn,]
# p = sub('_.*','',colnames(gbm))
for (i in names(dlist)) {
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]][gn,],project=i))
}
# for (sp in unique(p)) {
#   dlist[[sp]] <- FindVariableFeatures(CreateSeuratObject(countfunc(gbm[,p==sp]),project=sp))
# }
# reference.vec = c('2018_Rosenberg','2018_Zeisel','2018_Li_Science:fetal','2018_Li_Science:adult','2018_Lake_NatBiotech:VisualCortex','2018_Lake_NatBiotech:FrontalCortex','2018_Lake_NatBiotech:CerebellarHem','2018_Fan_CellResearch','2016_La:ESMoleculeCounts','2016_La:EmbryoMoleculeCounts')
# d.anchors <- FindIntegrationAnchors(object.list = dlist, reference = grep('GBM',names(dlist),invert = T))
# d.integrated <- IntegrateData(anchorset = d.anchors)
d.anchors <- FindIntegrationAnchors(object.list = dlist)
d.integrated <- IntegrateData(anchorset = d.anchors)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/atlasOnly_sub/',showWarnings = F, recursive = T)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/atlasOnly_sub/integrated.rds')
