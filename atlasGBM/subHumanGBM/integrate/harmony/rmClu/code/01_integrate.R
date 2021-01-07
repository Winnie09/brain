library(Seurat)
library(harmony)
library(data.table)
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/brain')
ddir <- 'atlasGBM/humanGBM/integrate/harmony/36nonNormal_allen10xOnly_seuratGene/res/'
rdir <- 'atlasGBM/subHumanGBM/integrate/harmony/rmClu/result/'
library(Seurat)
seu <- readRDS(paste0(ddir, 'atlasGBM_harmony.rds'))
str(seu)
clu.sl <- c(14,15,3,18,29,31)

active.ident <- as.character(seu@active.ident)
names(active.ident) <- names(seu@active.ident)
str(active.ident)
# seu@active.ident <- factor(sapply(active.ident, function(i) ifelse(i %in% clu.sl, i, NA)))
# DimPlot(seu)

tmp <-sapply(active.ident, function(i) ifelse(!i %in% clu.sl, i, NA))
names(tmp) <- names(active.ident)
cell.sl <- names(tmp[!is.na(tmp)])

mat <- seu@assays$RNA@counts
meta <- seu@meta.data[, c(-1:-3)]

mat <- mat[, cell.sl]
meta <- meta[cell.sl, ]

################
nfeatures = 2000
brain <- CreateSeuratObject(counts = mat, project = "atlasGBM") ## input is log-normalized => no more Normalized().
brain@meta.data <- meta
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
saveRDS(brain,paste0(rdir, 'sub_atlasGBM_harmony_tmp.rds'))

brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20, n.epochs=500) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
saveRDS(brain,paste0(rdir, 'sub_atlasGBM_harmony.rds'))
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'sub_umap_embeddings.rds'))

print(sessionInfo())




