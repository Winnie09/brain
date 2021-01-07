library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/humanatlas/integrate/seurat/noAllenSeuratGene/res/'##### 
pdir <- './atlasGBM/humanatlas/integrate/seurat/noAllenSeuratGene/plot/'##### 
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
mat =  readRDS('./atlasGBM/humanatlas/data/combine_mat.rds')
meta =  readRDS('./atlasGBM/humanatlas/data/meta_allcell.rds')
rownames(meta) <- meta$cell
selectcell <- sub('.*;', '', meta$cell[!grepl('2014_Pollen_NatBiotech', meta$study) & !grepl('allen', meta$study)])
mat <- mat[,selectcell]
mat <- mat[, colSums(mat > 0) > 500]
mat <- mat[rowMeans(mat > 0) >= 0.01,]
mat <- mat[!grepl('^MT-',rownames(mat)),]
meta <- meta[colnames(mat), ]
nfeatures <- 2000
brain <- CreateSeuratObject(counts = mat, project = "atlasGBM") ## input is log-normalized => no more Normalized().
brain@meta.data <- cbind(brain@meta.data, meta)
brain@meta.data$mydonor <- paste0(meta$study, ';',meta$gender,';', meta$age,';', meta$location, ';',meta$time, ';',meta$donor)
print(length(unique(brain@meta.data$mydonor)))
print(sort(unique(brain@meta.data$mydonor)))

brain.ls <- SplitObject(brain,split.by = 'mydonor')
print(length(brain.ls))
for(i in 1:length(brain.ls)){
    brain.ls[[i]] <- FindVariableFeatures(brain.ls[[i]],selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    print(i)
}

anchors <- FindIntegrationAnchors(object.list = brain.ls, k.filter=k.filter)
integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)

brain <- RunUMAP(integrated, reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony.rds'))
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))


DimPlot(brain, reduction = 'umap',label = T)
ggsave(paste0(pdir, 'umap_dimplot.pdf'),height = 10,width =12)

umap.dat <- data.frame(brain@reductions$umap@cell.embeddings)
umap.dat$seurat_clusters <- Idents(brain)
umap.dat$cell <- rownames(umap.dat)
saveRDS(umap.dat,paste0(rdir,'umap_with_clusters.rds'))

print(sessionInfo())

