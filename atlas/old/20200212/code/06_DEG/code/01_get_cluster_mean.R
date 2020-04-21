countfunc <- function(d) {
  log2(t(t(d)/(colSums(d) / 10000) + 1))
}

aggregateFunc <- function(m){
  clus = clu[match(colnames(m), names(clu))]
  me <- sapply(unique(clus),function(i) rowMeans(m[,clus==i,drop=F]))
  colnames(me) <- unique(clus)
  me = me[,sort(colnames(me),decreasing=F)]
}

clu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/cluster/cluster.rds')
## load atlas
dlist = list()
dlist[['2016_La:EmbryoMoleculeCounts']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds')))
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['2016_La:ESMoleculeCounts']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/ESMoleculeCounts.rds')))
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

#dlist[['2016_La:iPSMoleculeCounts']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/iPSMoleculeCounts.rds'))
dlist = list()
dlist[['2018_Fan_CellResearch']] <- aggregateFunc(log2(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')+1))
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['2018_Lake_NatBiotech:CerebellarHem']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/CerebellarHem.rds')))
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['2018_Lake_NatBiotech:FrontalCortex']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/FrontalCortex.rds')))
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['2018_Lake_NatBiotech:VisualCortex']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/VisualCortex.rds')))
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['2018_Li_Science:adult']] <- aggregateFunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/adult.rds'))
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['2018_Li_Science:fetal']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds')))
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['2018_Rosenberg']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/proc/expr.rds')))
row.names(dlist[['2018_Rosenberg']]) <- toupper(row.names(dlist[['2018_Rosenberg']]))
dlist[['2018_Rosenberg']] <- dlist[['2018_Rosenberg']][!duplicated(row.names(dlist[['2018_Rosenberg']])),]
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['2018_Zeisel']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/mat.rds')))
row.names(dlist[['2018_Zeisel']]) <- toupper(row.names(dlist[['2018_Zeisel']]))
dlist[['2018_Zeisel']] <- dlist[['2018_Zeisel']][!duplicated(row.names(dlist[['2018_Zeisel']])),]
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

dlist = list()
dlist[['allen']] <- aggregateFunc(countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds'))
)
saveRDS(dlist, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/sub/DEG/cluster_mean/',names(dlist)))

