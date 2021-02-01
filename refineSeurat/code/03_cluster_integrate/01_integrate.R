library(Seurat)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')

mat1 <-  log2(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')+1)
row.names(mat1) <- toupper(row.names(mat1)) ## mouse
mat1 <- mat1[!duplicated(row.names(mat1)),]
mat1 = mat1[rowSums(mat1>0)>0.01,]
gs = findVariableGene(mat1,1000)
set.seed(12345)
clu1 = kmeans(t(mat1[gs,]), 1000)$cluster

mat2 <- log2TPM_from_fludigm_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds'))
row.names(mat2) <- toupper(row.names(mat2))  ## mouse
mat2 <- mat2[!duplicated(row.names(mat2)),]
mat2 = mat2[rowSums(mat2>0)>0.01,]
gs = findVariableGene(mat2,1000)
set.seed(12345)
clu2 = kmeans(t(mat2[gs,]), 1000)$cluster


dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate',recursive = T)
saveRDS(clu1, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/clu1.rds')
saveRDS(clu2, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/clu2.rds')

data1 = sapply(1:max(clu1), function(c) rowMeans(mat1[,names(clu1[clu1==c]),drop=F]) )
data2 = sapply(1:max(clu2), function(c) rowMeans(mat2[,names(clu2[clu2==c]),drop=F]) )
colnames(data1) = paste0('2018_Fan_CellResearch:clu',1:max(clu1))
colnames(data2) = paste0('2016_La:clu',1:max(clu2))

dlist = list()
dlist[['2018_Fan_CellResearch']] = data1
dlist[['2016_La']] = data2

gn <- table(unlist(sapply(dlist,row.names)))
gn <- names(gn)[gn==length(dlist)]
dlist[[1]] = dlist[[1]][gn,]
dlist[[2]] = dlist[[2]][gn,]
mat = do.call(cbind,dlist)
saveRDS(mat,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/mat_beforeSeurat.rds')

saveRDS(cbind(mat1[rownames(mat),],mat2[rownames(mat),]), '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/original_matrix.rds')

for (i in names(dlist)) {
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]][gn,],project=i))
}
saveRDS(dlist,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/dlist.rds')

d.anchors <- FindIntegrationAnchors(object.list = dlist)
d.integrated <- IntegrateData(anchorset = d.anchors)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/',showWarnings = F, recursive = T)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/integrated.rds')

