## /home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/code/
library(Seurat)
library(future)
library(future.apply)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/'  ##### 
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
mat =  readRDS('./atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')
d <- readRDS(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/proc/tumor/doublet/doublet.rds'))  
#normcelllist <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/res/nonmalignant_cells_high_confident.rds')
nfeatures <- 2000
for (n in names(mat)) {
  tmp <- mat[[n]]
  tmp <- tmp[,!colnames(tmp) %in% d[[n]]]
  tmp <- tmp[!grepl('^MT-',rownames(tmp)),]
  tmp <- tmp[rowMeans(tmp > 0) >= 0.01,]
  tmp <- CreateSeuratObject(counts = tmp)
  tmp <- FindVariableFeatures(tmp,selection.method = "vst", nfeatures = nfeatures)
  mat[[n]] <- tmp
}

integrated.features <- SelectIntegrationFeatures(mat, nfeatures = nfeatures)
print(sort(integrated.features))
print(length(integrated.features))
saveRDS(integrated.features, paste0(rdir, nfeatures,'features.rds'))

for (n in names(mat)) {
  x <- mat[[n]]
  x <- ScaleData(x, features = integrated.features, verbose = FALSE)
  mat[[n]] <- RunPCA(x, features = integrated.features, verbose = FALSE)
}

anchors <- FindIntegrationAnchors(object.list = mat, reduction = "rpca")
mat <- IntegrateData(anchorset = anchors)

DefaultAssay(mat) <- "integrated"
mat <- ScaleData(mat, verbose = FALSE)
mat <- RunPCA(mat, npcs = 30, verbose = FALSE)
mat <- RunUMAP(mat, reduction = "pca", dims = 1:30,return.model=TRUE)
mat <- FindNeighbors(mat, reduction = "pca", dims = 1:30)
mat <- FindClusters(mat, resolution = 0.8)
saveRDS(mat,paste0(rdir, 'humanAtlas_harmony.rds'))
umap <- mat@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))


