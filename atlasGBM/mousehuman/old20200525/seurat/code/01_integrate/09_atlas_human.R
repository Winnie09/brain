

#############
library(Seurat)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
## load atlas
dlist = list()
print('1')
a <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/adult.rds'), platform = '10x')
b <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds'),  platform = '10x')
rn <- intersect(rownames(a), rownames(b))
dlist[['2018_Li_Science']] <- cbind(a[rn,],b[rn,])
print('4')
a <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/CerebellarHem.rds'), platform = '10x')
b <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/FrontalCortex.rds'), platform = '10x')
c <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/VisualCortex.rds'), platform = '10x')
rn <- intersect(rownames(a), rownames(b))
rn <- intersect(rn, rownames(c))
dlist[['2018_Lake_NatBiotech']] <- cbind(a[rn,], b[rn,], c[rn,])
print('5')
dlist[['2018_Fan_CellResearch']] <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds'), normalizeByLibsize=FALSE, log.transform.only = T) 
print('6')
dlist[['allen']] <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds'))
print('7')
dlist[['2017_Nowakowski_Science']] <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2017_Nowakowski_Science/proc/log2tpm.rds'), normalizeByLibsize=F)  ## add
print('9')
a <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds'))
print('9')
b <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/ESMoleculeCounts.rds'))
rn <- intersect(rownames(a), rownames(b))
dlist[['2016_La']] <- cbind(a[rn, ], b[rn, ])
print('11')
dlist[['2018_zhong_nature']] <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_zhong_nature/proc/expr.rds'), normalizeByLibsize = FALSE, log.transform.only = TRUE)

## select intersect genes
gn <- sapply(dlist,row.names)
gn <- table(unlist(gn))
gn <- names(gn)[gn==length(dlist)]
for (i in 1:length(dlist)) dlist[[i]] <- dlist[[i]][gn,]

for (i in names(dlist)) {
  print('findvariablefeatures...')
  print(i)
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]],project=i))
}
print('18')
d.anchors <- FindIntegrationAnchors(object.list = dlist)
print('19')
d.integrated <- IntegrateData(anchorset = d.anchors)
print('20')

a <- unlist(sapply(dlist,colnames))
v <- rep('0',length(a))
names(v) <- a
for (i in names(dlist)) v[colnames(dlist[[i]])] <- i
d.integrated@active.ident <- as.factor(v)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/result/atlas_human/',showWarnings = F, recursive = T)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/result/atlas_human/integrated.rds')

