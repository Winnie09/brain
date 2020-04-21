library(Seurat)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')

## load data
dlist = list()
dlist[['2018_Fan_CellResearch']] <- log2(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')+1)
dlist[['2018_Fan_CellResearch']] <- dlist[['2018_Fan_CellResearch']][!duplicated(row.names(dlist[['2018_Fan_CellResearch']])),]


dlist[['2016_La:EmbryoMoleculeCounts']] <- log2TPM_from_fludigm_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds'))
row.names(dlist[['2016_La_EmbryoMoleculeCounts']]) <- toupper(row.names(dlist[['2016_La_EmbryoMoleculeCounts']]))  ## mouse
dlist[['2016_La_EmbryoMoleculeCounts']] <- dlist[['2016_La_EmbryoMoleculeCounts']][!duplicated(row.names(dlist[['2016_La_EmbryoMoleculeCounts']])),]

gn <- table(unlist(sapply(dlist,row.names)))
gn <- names(gn)[gn==length(dlist)]

for (i in names(dlist)) {
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]][gn,],project=i))
}

d.anchors <- FindIntegrationAnchors(object.list = dlist)
d.integrated <- IntegrateData(anchorset = d.anchors)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/01_integrate/',showWarnings = F, recursive = T)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/01_integrate/integrated.rds')



