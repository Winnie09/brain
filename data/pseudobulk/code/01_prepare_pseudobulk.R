library(Seurat)
library(Matrix)
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/pseudobulk/data/'

## ===========================================================================
## mouse atlas each cluster is a pseudobulk: output: gene by pseudobulk matrix
## ===========================================================================
mouseatlas <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/humanAtlas_harmony.rds')
ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/celltype_annotation.rds')
expr = mouseatlas@assays$RNA@counts
rm(mouseatlas)

identical(names(ct), colnames(expr))
summary(expr[,1])
summary(expr[1,])
range(expr)

mouseatlaspb <- sapply(unique(ct),function(i) {
 rowMeans(expr[,ct==i])
})
saveRDS(mouseatlaspb, paste0(rdir, 'mouseatlas_pb.rds'))

## ===================================================================================
## human allen 10x each cell type as a pseudobulk: output: gene by pseudobulk matrix
## ===================================================================================
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
atlas <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/proc/count.rds'), platform='10x', log2normalize = TRUE, species='human')
colnames(atlas) <- gsub('allen:', 'allen_human_10x:', colnames(atlas))
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/data/meta_allcell.rds')
meta[,1] <- sub('allen_human_10x;', '', meta[,1])
rownames(meta) <- meta[,1]
meta <- meta[colnames(atlas), ]

humanatlaspb <- sapply(unique(meta$celltype),function(i) {
 rowMeans(atlas[,meta$celltype==i])
})
saveRDS(humanatlaspb, paste0(rdir, 'humanatlas_pb.rds'))

## ========================================
## GBM each tumor samples is a pseudo bulk
## ========================================
gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')

expr <- gbm$RNA@data
pat <- gbm@meta.data$orig.ident
clu <- Idents(gbm)
gbmpb <- sapply(unique(pat),function(i) {
  tmp <- expr[,pat==i]
  rowMeans(sapply(unique(clu[colnames(tmp)]),function(j) {
    rowMeans(tmp[,clu[colnames(tmp)]==j,drop=F]) ##   
  }))
})
colnames(gbmpb) <- as.character(unique(pat))
saveRDS(gbmpb, paste0(rdir, 'gbm_pb.rds'))

## =====================
## human cortex atlas
## ===================
library(Seurat)
tb = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/humancortexatlas_celltype_annotation.csv', as.is = T)


atlas <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/humanAtlas_harmony.rds')


expr = atlas@assays$RNA@counts
actide = as.numeric(atlas@active.ident)-1
n = names(atlas@active.ident)
ct = tb[match(actide, tb[,1]),3]
names(ct) = n
saveRDS(ct, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/celltype.rds')

identical(names(ct), colnames(expr))
summary(expr[,1])
summary(expr[1,])
range(expr)

atlaspb <- sapply(unique(ct),function(i) {
 rowMeans(expr[,ct==i])
})
saveRDS(atlaspb, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/pseudobulk/data/humancortexatlas_pb.rds')




