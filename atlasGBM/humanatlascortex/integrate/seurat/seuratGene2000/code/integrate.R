library(Seurat)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/'  ##### 
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
mat =  readRDS('./atlasGBM/humanatlascortex/data/dlist.rds')

nfeatures <- 2000
for (n in names(mat)) {
  tmp <- mat[[n]]
  tmp <- tmp[, colSums(tmp > 0) > 500]
  tmp <- tmp[!grepl('^MT-',rownames(tmp)),]
  tmp <- tmp[rowMeans(tmp > 0) >= 0.01,]
  tmp <- t(t(tmp)/apply(tmp,2,max))*10
  tmp <- CreateSeuratObject(counts = tmp)
  tmp <- FindVariableFeatures(tmp,selection.method = "vst", nfeatures = nfeatures)
  mat[[n]] <- tmp
}

integrated.features <- SelectIntegrationFeatures(mat, nfeatures = nfeatures)
print(sort(integrated.features))
print(length(integrated.features))
saveRDS(integrated.features, paste0(rdir, nfeatures,'features.rds'))

anc <- FindIntegrationAnchors(object.list = mat, anchor.features = integrated.features)
mat <- IntegrateData(anchorset = anc)

DefaultAssay(mat) <- "integrated"
mat <- ScaleData(mat, verbose = FALSE)
mat <- RunPCA(mat, npcs = 30, verbose = FALSE)
mat <- RunUMAP(mat, reduction = "pca", dims = 1:30,return.model=TRUE)
mat <- FindNeighbors(mat, reduction = "pca", dims = 1:30)
mat <- FindClusters(mat, resolution = 0.8)

saveRDS(mat,paste0(rdir, 'humanAtlas_harmony.rds'))
umap <- mat@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))

