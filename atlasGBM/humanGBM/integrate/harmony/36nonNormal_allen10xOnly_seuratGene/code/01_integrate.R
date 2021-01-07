library(Seurat)
library(harmony)
library(data.table)
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
rdir = 'atlasGBM/humanGBM/integrate/harmony/36nonNormal_allen10xOnly_seuratGene/res/'  ##### 
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
nfeatures <- 2000
mat <- readRDS('atlasGBM/humanGBM/data/36nonNormalGBM_humanatlas_combined_log2norm_mat.rds')
meta <- readRDS('atlasGBM/humanGBM/data/36nonNormalGBM_humanatlas_combined_meta.rds')
rownames(meta) <- meta$cell
meta <- meta[colnames(mat), ]

meta <- meta[!grepl('allen_human_fluidigm', meta$study),]
mat = mat[, meta$cell]
mat <- mat[!grepl('^MT', rownames(mat)), ]
mat <- mat[rowMeans(mat > 0) > 0.01, ]

meta$study <- sapply(1:nrow(meta), function(i){
  if (meta$study[i] == 'GBM'){
    return(meta$Sample.ID[i])
  } else {
    return(meta$study[i])
  }
})
table(meta$study)

brain <- CreateSeuratObject(counts = mat, project = "atlasGBM") ## input is log-normalized => no more Normalized().
brain@meta.data <- data.frame(brain@meta.data, meta, normalTumor = ifelse(grepl('GBM', meta$study), 'GBM', 'atlas'))
brain.ls <- SplitObject(brain, split.by = 'study')

for(i in 1:length(brain.ls)){
    brain.ls[[i]] <- FindVariableFeatures(brain.ls[[i]],selection.method = "vst", nfeatures = nfeatures)
    print(i)
}
integrated.features <- SelectIntegrationFeatures(brain.ls, nfeatures = nfeatures)
print(sort(integrated.features))
print(length(integrated.features))
saveRDS(integrated.features,paste0(rdir, nfeatures,'features.rds'))

brain <- ScaleData(brain,features = integrated.features,verbose = FALSE)
brain <- RunPCA(brain,npcs = 30,features = integrated.features, verbose = FALSE)
brain <- RunHarmony(brain, c("normalTumor", "study"))
saveRDS(brain,paste0(rdir, 'atlasGBM_harmony_tmp.rds'))

brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20, n.epochs=500) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
saveRDS(brain,paste0(rdir, 'atlasGBM_harmony.rds'))
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))

print(sessionInfo())



