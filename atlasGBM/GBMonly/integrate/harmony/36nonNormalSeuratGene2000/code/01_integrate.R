## /home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/code/
library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/'  ##### 
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
mat =  readRDS('./atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_mat.rds')
meta =  readRDS('./atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
rownames(meta) <- meta$cell
meta <- meta[colnames(mat), ]
nfeatures <- 2000
mat <- mat[!grepl('^MT-',rownames(mat)),]
mat <- mat[rowMeans(mat > 0) >= 0.01,]

brain <- CreateSeuratObject(counts = mat, project = "atlasGBM") ## input is log-normalized => no more Normalized().
brain@meta.data <- cbind(brain@meta.data, meta)
brain.ls <- SplitObject(brain, split.by = 'Sample.ID')
print(length(brain.ls))
for(i in 1:length(brain.ls)){
    brain.ls[[i]] <- FindVariableFeatures(brain.ls[[i]],selection.method = "vst", nfeatures = nfeatures)
    print(i)
}

integrated.features <- SelectIntegrationFeatures(brain.ls, nfeatures = nfeatures)
print(sort(integrated.features))
print(length(integrated.features))
saveRDS(integrated.features, paste0(rdir, nfeatures,'features.rds'))

brain <- ScaleData(brain,features = integrated.features,verbose = FALSE)
brain <- RunPCA(brain,npcs = 30,features = integrated.features, verbose = FALSE)
brain <- RunHarmony(brain, c("Sample.ID"))
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony_tmp.rds'))

brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony.rds'))
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))


