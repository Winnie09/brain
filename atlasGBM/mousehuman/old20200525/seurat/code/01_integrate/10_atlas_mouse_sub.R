library(Seurat)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
## load atlas
dlist = list()
## load atlas
dlist = list()
print('1')
dlist[['2018_Zeisel']] <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/sub10000mat.rds'), platform = '10x', maxCell = 1e4)
print('2')
dlist[['2018_Rosenberg']] <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/proc/expr.rds'), platform = '10x', maxCell = 1e4) ## 10x count
print('3')
dlist[['2016_Tasic_NatNeuro']] <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_Tasic_NatNeuro/proc/log2tpm.rds'), normalizeByLibsize = FALSE, maxCell = 1e4)
print('9')
c <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/MouseEmbryoMoleculeCounts.rds'), maxCell = 1e4)
print('11')
d <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/MouseAdultDAMoleculeCounts.rds'), maxCell = 1e4)
rn <- intersect(rownames(c), rownames(d))
dlist[['2016_La']] <- cbind(c[rn, ], d[rn, ])
print('10')
dlist[['2015_Zeisel_Science']] <- preprocess(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Zeisel_Science/proc/log2tpm.rds'), normalizeByLibsize = FALSE, maxCell = 1e4)


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

dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/atlas_mouse_sub/',showWarnings = F, recursive = T)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/atlas_mouse_sub/integrated.rds')
