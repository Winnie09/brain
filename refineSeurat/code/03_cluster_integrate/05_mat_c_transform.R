library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
source('./resource/function.R')
d.integrated = readRDS('./brain/refineSeurat/result/03_cluster_integrate/integrated.rds')
mat1 = as.matrix(d.integrated@assays$integrated@data)
mat2 = as.matrix(d.integrated@assays$RNA@data)
c = mat2[rownames(mat1),] - mat1

mat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/original_matrix.rds')
clu1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/clu1.rds')
clu2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/clu2.rds')

cell = colnames(mat)
study = sapply(cell, function(i) sub(':.*','',i))
clu = sapply(cell, function(i) {
  if (i %in% names(clu1)){
    clu1[which(names(clu1) == i)]
  } else {
    clu2[which(names(clu2) == i)]
  }
} )
df = data.frame(cell = cell, study = study, clu = clu, stringsAsFactors = F)
saveRDS(df, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/clu_meta.rds')
c2 <- sapply(colnames(mat),function(i){
  study = df[which(df$cell == i), 'study']
  clu = df[which(df$cell == i), 'clu']
  c[, colnames(c) %in% paste0(study,':clu',clu)]
  })
mat_transform = mat[rownames(c2),] - c2
saveRDS(mat_transform,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/mat_c_transform.rds')
