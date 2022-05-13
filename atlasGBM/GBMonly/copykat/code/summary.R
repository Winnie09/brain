setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/copykat/res')
cancer <- sapply(list.files(),function(i) {
  a <- readRDS(i)
  rownames(a)[a[,2]=='aneuploid']
})

normal <- sapply(list.files(),function(i) {
  a <- readRDS(i)
  rownames(a)[a[,2]=='diploid']
})

sapply(cancer,length)/(sapply(normal,length)+sapply(cancer,length))

clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')

table(clu[unlist(cancer)])/table(clu)

