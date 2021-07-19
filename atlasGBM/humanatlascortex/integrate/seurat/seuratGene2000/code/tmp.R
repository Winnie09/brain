library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/'  ##### 
mat <- readRDS(paste0(rdir, 'humanAtlas_harmony.rds'))
mat <- RunUMAP(mat, reduction = "pca", dims = 1:30,return.model=TRUE)
mat <- FindNeighbors(mat, reduction = "pca", dims = 1:30)
mat <- FindClusters(mat, resolution = 0.8)

saveRDS(mat,paste0(rdir, 'humanAtlas_harmony.rds'))
umap <- mat@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))

