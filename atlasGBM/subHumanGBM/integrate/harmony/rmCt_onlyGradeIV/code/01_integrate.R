library(Seurat)
library(harmony)
library(data.table)
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/brain')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
rdir <- 'atlasGBM/subHumanGBM/integrate/harmony/rmCt_onlyGradeIV/result/'
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
nfeatures <- 2000

atlas.gbm <- readRDS('atlasGBM/humanGBM/integrate/harmony/36nonNormal_allen10xOnly_seuratGene/res/atlasGBM_harmony_with_ct.rds')
colnames(atlas.gbm@meta.data)[ncol(atlas.gbm@meta.data)] <- 'celltype_anno'
ct <- atlas.gbm@meta.data[,'celltype_anno']
id <- intersect(which(!ct %in% c('endothelial', 'microglia', 'oligodendrocyte')),
                which(atlas.gbm@meta.data[,'Tumor.Grade'] == 'IV'))
mat <- atlas.gbm@assays$RNA@counts[, id]

brain <- CreateSeuratObject(counts = mat, project = "atlasGBM") ## input is log-normalized => no more Normalized().
brain@meta.data <- cbind(brain@meta.data[,1:3], atlas.gbm@meta.data[id, c(-1:-3)])

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

brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
saveRDS(brain,paste0(rdir, 'atlasGBM_harmony.rds'))
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))

print(sessionInfo())



